package main

import (
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
	maxSnap := len(mergers.Rvir)
	mpeak := CalcMpeak(mergers)
	
	tags := lib.NewTags(mergers.Haloes)
	
	nWorkers := lib.SetThreads(-1)
	fileName := fmt.Sprintf(snapFmt, maxSnap, 0)
	header := lib.OpenGadget2Zoom(fileName, []string{"x", "v", "id32"})
	
	idxList := lib.NewCompactList(int32(header.NTot[HRLevel]))
	snapList := lib.NewCompactList(int32(header.NTot[HRLevel]))
	workers := NewWorkerArray(nWorkers, header.NTot[HRLevel], header.L)
	
	for snap := 0; snap <= maxSnap; snap++ {
		log.Printf("Snapshot %3d", snap)
		lib.MemoryUsage()
		
		hx, hr := ExtractHaloes(mergers, snap)
		
		for b := 0; b < cfg.Blocks[cfgi]; b++ {
			log.Printf("    Block %d", b)

			fileName := fmt.Sprintf(snapFmt, snap, b)
			px, pid := ReadLevel(fileName, HRLevel)
			
			TagParticles(workers, px, hx, hr, mpeak)
			lib.InsertOwnersInLists(workers, snap, idxList, snapList, mpeak)
			tags.AddChangedParticles(pid, idxList, snapList, snap)
		}

		runtime.GC()
	}

	lib.OrderTags(tags)
	lib.WriteTags(baseDir, NFiles, tags)
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
	lib.ThreadSplitArray(nTot, nTot, func(worker, start, end, step int) {
		workers[worker] = lib.NewTagWorker(float32(L), start, end, step)
	})
	return workers
}

func TagParticles(
	w []*lib.TagWorker, px, hx [][3]float32, hr, mpeak []float32,
) {
	
	lib.ThreadSplitArray(len(w), len(w), func(worker, _, _, _ int) {
		w[worker].LoadParticles(px)
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
		lib.Mpeak(m.Mvir[i])
	}
	return mpeak
}
