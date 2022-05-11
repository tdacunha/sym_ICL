package main

import (
	"time"
	"fmt"
	"log"
	"os"
	"runtime"
	"strconv"
	"github.com/phil-mansfield/lmc_ges_tracking/lib"
)

var (
	// What level in the file are the HR particles stored at?
	HRLevel = 1
	NFiles = 8
)

func main() {
	haloIndex := -1
    if len(os.Args) > 3 || len(os.Args) <= 1 {
        panic("You must supply a config file and an (optional) index to " +
			"restrict analysis to.")
    } else if len(os.Args) == 3 {
		var err error
		haloIndex, err = strconv.Atoi(os.Args[2])
		if err != nil { panic(err.Error()) }
	}

	inputName := os.Args[1]
	
	cfg := lib.ParseConfig(inputName)
	
	for i := range cfg.BaseDir {
		if haloIndex != -1 && i != haloIndex { continue }
		TagHalo(cfg, i)
		runtime.GC()
	}
}

func TagHalo(cfg *lib.Config, cfgi int) {
	log.Printf("Tagging halo %d (%s)", cfgi, cfg.BaseDir[cfgi])
	lib.MemoryUsage()

	baseDir, snapFmt := cfg.BaseDir[cfgi], cfg.SnapFormat[cfgi]
	
	mergers := lib.ReadMergers(lib.MergerFileName(baseDir))

	maxSnap := len(mergers.Rvir[0]) - 1
	mpeak := CalcMpeak(mergers)
	
	tags := lib.NewTags(mergers.Haloes)
	
	nWorkers := lib.SetThreads(-1)
	fileName := fmt.Sprintf(snapFmt, maxSnap, 0)
	header := lib.OpenGadget2Zoom(fileName, []string{"x", "v", "id32"})
	
	idxList := lib.NewCompactList(int32(header.NTot[HRLevel]))
	snapList := lib.NewCompactList(int32(header.NTot[HRLevel]))
	workers := NewWorkerArray(nWorkers, header.NTot[HRLevel], header.L)
	
	var (
		dt1, dt2, dt3, dt4 time.Duration
	)

	for snap := 0; snap <= maxSnap; snap++ {
		log.Printf("Snapshot %3d", snap)
		if snap % 10 == 0 {
			lib.MemoryUsage()
			log.Printf(`
    Reading:   %.2f
    Tagging:   %.2f
    Inserting: %.2f
    Combining: %.2f
`, float64(dt1)/1e9, float64(dt2)/1e9, float64(dt3)/1e9, float64(dt4)/1e9,
			)
		}
		
		hx, hr := ExtractHaloes(mergers, snap)

		for b := 0; b < cfg.Blocks[cfgi]; b++ {
			fileName := fmt.Sprintf(snapFmt, snap, b)
			t0 := time.Now()
			px, pid := ReadLevel(fileName, HRLevel)

			t1 := time.Now()
			TagParticles(workers, px, pid, hx, hr, mpeak)
			t2 := time.Now()
			lib.InsertOwnersInLists(workers, snap, idxList, snapList, mpeak)
			t3 := time.Now()
			tags.AddChangedParticles(idxList, snapList, snap)
			t4 := time.Now()

			dt1 += t1.Sub(t0)
			dt2 += t2.Sub(t1)
			dt3 += t3.Sub(t2)
			dt4 += t4.Sub(t3)
		}
		runtime.GC()
	}

	lib.OrderTags(tags)
	lookup := lib.NewTagLookup(header.NTot[HRLevel])
	lookup.AddTags(tags)
	lib.WriteTags(baseDir, NFiles, tags, lookup)

	EstimateSpaceSavings(tags)
}

func ExtractHaloes(m *lib.Mergers, snap int) (hx [][3]float32, hr []float32) {
	hx, hr = make([][3]float32, len(m.X)), make([]float32, len(m.X))
	for ih := range m.X {
		hx[ih], hr[ih] = m.X[ih][snap], m.Rvir[ih][snap]
	}
	return hx, hr
}

func NewWorkerArray(nWorkers, nTot int, L float64) []*lib.TagWorker {
	workers := make([]*lib.TagWorker, nWorkers)
	lib.ThreadSplitArray(nTot, nWorkers, func(worker, start, end, step int) {
		workers[worker] = lib.NewTagWorker(float32(L), start, end, step)
	})
	return workers
}

func TagParticles(
	w []*lib.TagWorker, px [][3]float32, pid []int32,
	hx [][3]float32, hr, mpeak []float32,
) {
	lib.ThreadSplitArray(len(px), len(w), func(worker, start, end, skip int) {
		w[worker].ResetBounds(start, end, skip)
		w[worker].LoadParticles(px, pid)
		w[worker].FindParticleOwners(hx, hr, mpeak)
	})
}

func ReadLevel(fileName string, level int) (xp [][3]float32, idp []int32) {
	f := lib.OpenGadget2Zoom(fileName, []string{"x", "v", "id32"})

	xp = make([][3]float32, f.N[level])
	idp = make([]int32, f.N[level])

	f.Read("x", level, xp)
	f.Read("id32", level, idp)

	return xp, idp
}

func CalcMpeak(m *lib.Mergers) []float32 {
	mpeak := make([]float32, len(m.Mvir))
	for i := range mpeak {
		mpeak[i] = lib.Mpeak(m.Mvir[i])
	}
	return mpeak
}

func EstimateSpaceSavings(tags *lib.Tags) {
	nTot, n0 := 0, 0
	n0Infall := 0
	for i := range tags.Snap {
		nTot += len(tags.Snap[i])
		n0 += int(tags.N0[i])
		
		for j := range tags.Snap[i] {
			n0Infall += 236 - int(tags.Snap[i][j])
		}
	}

	f0 := float64(n0) / float64(nTot)
	fInfall := float64(n0Infall) / (236*float64(n0))
	log.Printf("f0 = %.3f, fin = %.3f, ftot = %.3f",
		f0, fInfall, 0.5*f0*fInfall)
}
