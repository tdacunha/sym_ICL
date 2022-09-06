package main

import (
	"fmt"
	"runtime"
	"os"
	"encoding/binary"
	"strconv"
	"log"
	
	"github.com/phil-mansfield/symphony_pipeline/lib"
)


// HaloTrack contains the evolution history of a single halo.
func main() {	
	haloIndex := -1
    if len(os.Args) > 3 || len(os.Args) <= 1 {
        panic("You must provide a config file.")
    } else if len(os.Args) == 3 {
		var err error
		haloIndex, err = strconv.Atoi(os.Args[2])
		if err != nil { panic(err.Error()) }
	}

	inputName := os.Args[1]
	cfg := lib.ParseConfig(inputName)

	for i := range cfg.BaseDir {
		if haloIndex != -1 && haloIndex !=  i { continue }

		log.Printf("Analysing %s", cfg.BaseDir[i])
		runtime.GC()

		h := ReadFullTree(cfg.BaseDir[i])
		log.Println("Creating lookup table")
		h.IDTable = NewLookupTable(h.ID)
		log.Println("Creating tracks")
		t := CalcTracks(h, cfg.MatchID[i], cfg.MaxSnap[i])
		log.Println("Writing tracks")
		WriteTracks(cfg.BaseDir[i], t)
		MemoryLog()
	}
	log.Println("Done")
}


func ReadFullTree(baseDir string) *Haloes {
	fileNames := lib.TreeFileNames(baseDir)
	files := make([]*os.File, len(fileNames))
	for i := range files {
		var err error
		files[i], err = os.Open(fileNames[i])
		if err != nil { panic(err.Error()) }
	}

	n := lib.TreeNTot(files)
	h := &Haloes{
		ID: make([]int32, n), DescID: make([]int32, n),
		UPID: make([]int32, n), DFID: make([]int32, n),
		Snap: make([]int32, n), Mvir: make([]float32, n),
		IDTable: nil, 
	}

	lib.ReadTreeVarFullInt(files, "ID", h.ID)
	lib.ReadTreeVarFullInt(files, "DescID", h.DescID)
	lib.ReadTreeVarFullInt(files, "UPID", h.UPID)
	lib.ReadTreeVarFullInt(files, "DFID", h.DFID)
	lib.ReadTreeVarFullInt(files, "Snap", h.Snap)
	lib.ReadTreeVarFullFloat(files, "Mvir", h.Mvir)

	return h
}

// Haloes is a collection of raw columns from the tree.dat files.
type Haloes struct {
	ID, DescID, UPID, DFID, Snap []int32
	Mvir []float32
	IDTable *LookupTable
}

type LookupTable struct {
	Idx []int32
}

func NewLookupTable(id []int32) *LookupTable {
	max := int32(-1)
	for i := range id { 
		if id[i] > max { max = id[i] }
	}
	tab := &LookupTable{ make([]int32, max+1) }
	for i := range tab.Idx { tab.Idx[i] = -2 }
	for i := range id {
		tab.Idx[id[i]] = int32(i)
	}
	return tab
}

func (tab *LookupTable) Find(id int32) int32 {
	return tab.Idx[id]
}

type Tracks struct {
	N, MWIdx int
	Starts, Ends []int32
	IsReal, IsDisappear, IsCentralSub []bool
	HostIdx, HostSnap, TrackIdx [][]int32
	IsReverseMerger, IsReverseSub, IsValidHost [][]bool
	Mpeak []float32
	MostMassivePreHostTrack []int32
	PreprocessSnap []int32
}

func WriteTracks(dir string, t *Tracks) {
	f, err := os.Create(lib.BranchesFileName(dir))
	if err != nil { panic(err.Error()) }
	defer f.Close()
	
	err = binary.Write(f, lib.ByteOrder, int32(t.N))
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, lib.ByteOrder, int32(t.MWIdx))
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, lib.ByteOrder, int32(0))
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, lib.ByteOrder, t.Ends)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, lib.ByteOrder, t.IsReal)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, lib.ByteOrder, t.IsDisappear)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, lib.ByteOrder, t.IsCentralSub)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, lib.ByteOrder, t.MostMassivePreHostTrack)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, lib.ByteOrder, t.PreprocessSnap)
}

func CalcTracks(h *Haloes, centralID int32, maxSnap int32) *Tracks {
	t := &Tracks{ }

	t.Starts, t.Ends = StartsEnds(h)
	t.N = len(t.Starts)
	t.MWIdx = FindCentral(h, t, centralID)
	t.IsReal = IsReal(h, t)
	t.IsDisappear = IsDisappear(h, t, maxSnap)
	t.IsCentralSub = IsCentralSub(h, t)
	t.HostIdx, t.HostSnap = FindAllHosts(h, t)	
	t.TrackIdx = FindTrackIndices(h, t)
	t.IsReverseMerger = IsReverseMerger(h, t)
	t.IsReverseSub = IsReverseSub(h, t)
	t.IsValidHost = IsValidHost(h, t)
	t.Mpeak = Mpeak(h, t)
	t.MostMassivePreHostTrack = MostMasssivePreHostTrack(h, t)
	t.PreprocessSnap = PreprocessSnap(h, t)

	return t
}

func MemoryLog() {
	ms := &runtime.MemStats{ }
	runtime.ReadMemStats(ms)
	log.Printf("Allocated: %.1f In Use: %.1f Idle: %.1f\n",
		float64(ms.Alloc) / 1e9, float64(ms.HeapInuse) / 1e9,
		float64(ms.HeapIdle) / 1e9)
}

func StartsEnds(h *Haloes) (starts, ends []int32) {
	starts, ends = []int32{ 0 }, []int32{ }
	for i := 0; i < len(h.DFID) - 1; i++ {
		if h.DescID[i+1] != h.ID[i] {
			starts = append(starts, int32(i+1))
			ends = append(ends, int32(i+1))
		}
	}
	ends = append(ends, int32(len(h.DFID)))

	return starts, ends
}

func FindCentral(h *Haloes, t *Tracks, centralID int32) int {
	fmt.Println(len(h.ID), len(t.Starts), centralID)
	for i := range t.Starts {
		if h.ID[t.Starts[i]] == centralID {
			return i
		}
	}
	panic(fmt.Sprintf("Halo with ID %d not found", centralID))
}

func IsReal(h *Haloes, t *Tracks) []bool {
	isReal := make([]bool, t.N)
	for i := range isReal {
		isReal[i] = h.UPID[t.Ends[i] - 1] == -1
	}
	return isReal
}

func IsDisappear(h *Haloes, t *Tracks, maxSnap int32) []bool {
	out := make([]bool, t.N)
	
	for i := range out {
		out[i] = h.DescID[t.Starts[i]] == -1 &&
			h.Snap[t.Starts[i]] != maxSnap
	}

	return out
}

func IsCentralSub(h *Haloes, t *Tracks) []bool {
	minID, maxID := h.DFID[t.Starts[t.MWIdx]], h.DFID[t.Ends[t.MWIdx] - 1]

	out := make([]bool, t.N)
	for i := range out {
		start, end := t.Starts[i], t.Ends[i]
		for j := start; j < end; j++ {
			if h.UPID[j] != -1 {
				k := h.IDTable.Find(h.UPID[j])
				out[i] = out[i] || h.DFID[k] >= minID && h.DFID[k] <= maxID
			}
		}
	}

	return out
}

func FindAllHosts(h *Haloes, t *Tracks) ([][]int32, [][]int32) {
	outIdx, outSnap := make([][]int32, t.N), make([][]int32, t.N)

	for i := range t.Starts {
		for j := t.Starts[i]; j < t.Ends[i]; j++ {
			if h.UPID[j] != -1 {
				k := h.IDTable.Find(h.UPID[j])
				outIdx[i] = append(outIdx[i], k)
				outSnap[i] = append(outSnap[i], h.Snap[j])
			}
		}
	}

	return outIdx, outSnap
}

func FindTrackIndices(h *Haloes, t *Tracks) [][]int32 {
	out := make([][]int32, t.N)
	idxTable := h.Snap // Not used after this

	for i := range t.Starts {
		for j := t.Starts[i]; j < t.Ends[i]; j++ {
			idxTable[j] = int32(i)
		}
	}

	for i := range t.HostIdx {
		out[i] = make([]int32, len(t.HostIdx[i]), len(t.HostIdx[i]))

		for j := range out[i] {
			out[i][j] = idxTable[t.HostIdx[i][j]]
		}
	}

	return out
}

func IsReverseMerger(h *Haloes, t *Tracks) [][]bool {
	out := make([][]bool, len(t.HostIdx))
	
	for i := range out {
		out[i] = make([]bool, len(t.HostIdx[i]), len(t.HostIdx[i]))

		minID, maxID := h.DFID[t.Starts[i]], h.DFID[t.Ends[i] - 1]
		for j := range t.HostIdx[i] {
			// The last halo of the host (first in DF order)
			start := t.Starts[t.TrackIdx[i][j]]
			// If the descendent is -1, that means this wasn't a merger
			// and we're okay
			if h.DescID[start] != -1 {
				// If the host's descendant has a DFID within this halo's DFID
				// range, that means the "host" merged with the descendant.
				k := h.IDTable.Find(h.DescID[start])
				out[i][j] = h.DFID[k] >= minID && h.DFID[k] <= maxID
			}
		}
	}

	return out
}

func IsReverseSub(h *Haloes, t *Tracks) [][]bool {
	out := make([][]bool, len(t.HostIdx))
	
	for i := range out {
		out[i] = make([]bool, len(t.HostIdx[i]), len(t.HostIdx[i]))

		minID, maxID := h.DFID[t.Starts[i]], h.DFID[t.Ends[i] - 1]
		for j := range t.HostIdx[i] {
			// Last halo of the host (first in DF order)
			start := t.Starts[t.TrackIdx[i][j]]

			// Loop from last host halo to the timestep after the host-sub
			// pairing was found
			for k := start; k < t.HostIdx[i][j]; k++ {
				// If this thing is a host halo, everything is fine
				if h.UPID[k] != -1 {
					// Check if you're a subhalo
					l := h.IDTable.Find(h.UPID[k])
					if h.DFID[l] >= minID && h.DFID[l] <= maxID {
						out[i][j] = true
						break
					}
				}
			}
		}
	}

	return out
}

func FindFirstIndex(h *Haloes, i int32) int32 {
	j := i
	for ; j > 0; j-- {
		if h.DescID[j] != h.ID[j-1] { return j }
	}
	return 0
}

func FindLastIndex(h *Haloes, i int32) int32 {
	j := int(i)
	for ; j < len(h.ID) - 1; j++ {
		if h.DescID[j+1] != h.ID[j] { return int32(j) }
	}
	
	return int32(len(h.ID) - 1)
}

func IsValidHost(h *Haloes, t *Tracks) [][]bool {
	out := make([][]bool, t.N)
	for i := range out {
		out[i] = make([]bool, len(t.TrackIdx[i]), len(t.TrackIdx[i]))
		for j := range out[i] {
			k := t.TrackIdx[i][j]
			out[i][j] = t.IsReal[k] && !t.IsDisappear[k] &&
				!t.IsReverseMerger[i][j] && !t.IsReverseSub[i][j]
			
		}
	}

	return out
}

func Mpeak(h *Haloes, t *Tracks) []float32 {
	out := make([]float32, t.N)
	for i := range out {
		for j := t.Starts[i]; j < t.Ends[i]; j++ {
			if h.Mvir[j] > out[i] { out[i] = h.Mvir[j] }
		}
	}

	return out
}

func MostMasssivePreHostTrack(h *Haloes, t *Tracks) []int32 {
	out := make([]int32, t.N)
	for i := range out { out[i] = -1 }
	for i := range out {
		for j := range t.TrackIdx[i] {
			if !t.IsValidHost[i][j] { continue }
			k := t.TrackIdx[i][j]
			if k == int32(t.MWIdx) {
				out[i] = -1
				continue
			}
			if (out[i] == -1 || t.Mpeak[k] > t.Mpeak[out[i]]) &&
				t.Mpeak[k] > t.Mpeak[i] {
				out[i] = k
			}
		}
	}

	return out
}

func PreprocessSnap(h *Haloes, t *Tracks) []int32 {
	out := make([]int32, t.N)
	for i := range out { out[i] = -1 }
	for i := range out {
		for j := range t.HostSnap[i] {
			if !t.IsValidHost[i][j] && t.TrackIdx[i][j] != int32(t.MWIdx) {
				continue }
			if out[i] == -1 || out[i] > t.HostSnap[i][j] {
				out[i] = t.HostSnap[i][j]
			}
		}
	}
	return out
}

func IDOrder(id []int32) []int32 {
	fid := make([]float64, len(id))
	for i := range id {
		fid[i] = float64(id[i])
	}

	idx := lib.QuickSortIndex(fid)
	idx32 := make([]int32, len(idx))
	for i := range idx32 { idx32[i] = int32(idx[i]) }
	
	return idx32
}
