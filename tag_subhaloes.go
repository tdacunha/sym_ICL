package main

import (
	"fmt"
	"sort"
	"runtime"
	"os"
	"encoding/binary"
	
	"github.com/phil-mansfield/lmc_ges_tracking/lib"
	"github.com/phil-mansfield/guppy/lib/catio"
)

var (
	MaxSnap = int32(235)
	TreeFileName="/scratch/users/enadler/Halo416/rockstar/trees/tree_0_0_0.dat"
	OutputFileName="/scratch/users/phil1/lmc_ges_tracking/Halo416/sub_info.dat"
)

// Haloes is a collection of raw columns from the tree.dat files.
type Haloes struct {
	ID, DescID, UPID, DFID, Snap []int32
	Mvir []float64
	IDTable *LookupTable
}

type LookupTable struct {
	Order, ID []int32
}

func NewLookupTable(id []int32) *LookupTable {
	order := IDOrder(id)
	return &LookupTable{ order, id }
}

func (t *LookupTable) Find(id int32) int32 {
	j := sort.Search(len(t.ID), func (j int) bool {
		return t.ID[t.Order[j]] >= id
	})
	return t.Order[j]
}

// HaloTrack contains the evolution history of a single halo.
func main() {	
	fmt.Println()
	
	cfg := catio.DefaultConfig
	cfg.SkipLines = 45

	//rd := catio.TextFile(TreeFileName, cfg)
	//
	// Read the columns
	//cols := rd.ReadInts([]int{ 1, 3, 6, 28, 31 })
	//id, descid, upid, dfid, snap := cols[0], cols[1], cols[2], cols[3],cols[4]
	//mvir := rd.ReadFloat64s([]int{ 10 })[0]
	//h := &Haloes{ ID: id, DescID: descid, UPID: upid,
	// DFID: dfid, Snap: snap, Mvir: mvir }
	h := ReadBinaryHaloes("tmp_data/tree.bin")
	SortHaloes(h)
	runtime.GC()
	t := CalcTracks(h)

	_ = t
}

func ReadBinaryHaloes(fname string) *Haloes {
	f, err := os.Open(fname)
	if err != nil { err.Error() }
	defer f.Close()

	n := int32(0)
	order := binary.LittleEndian
	err = binary.Read(f, order, &n)
	if err != nil { panic(err.Error()) }

	h := &Haloes{ }
	h.ID, h.DescID, h.UPID = make([]int32,n), make([]int32,n), make([]int32,n)
	h.DFID, h.Snap, h.Mvir = make([]int32,n), make([]int32,n), make([]float64,n)

	err = binary.Read(f, order, h.ID)
	if err != nil { panic(err.Error()) }
	err = binary.Read(f, order, h.DescID)
	if err != nil { panic(err.Error()) }
	err = binary.Read(f, order, h.UPID)
	if err != nil { panic(err.Error()) }
	err = binary.Read(f, order, h.DFID)
	if err != nil { panic(err.Error()) }
	err = binary.Read(f, order, h.Snap)
	if err != nil { panic(err.Error()) }
	err = binary.Read(f, order, h.Mvir)
	if err != nil { panic(err.Error()) }
	
	return h
}

func ReorderInts(x, order []int32) []int32 {
	out := make([]int32, len(order))
	for i := range out {
		out[i] = x[order[i]]
	}
	return out
}

func ReorderFloat64s(x []float64, order []int32) []float64 {
	out := make([]float64, len(order))
	for i := range out {
		out[i] = x[order[i]]
	}
	return out
}

func SortHaloes(h *Haloes) {
	order := IDOrder(h.DFID)
	h.ID = ReorderInts(h.ID, order)
	h.DescID = ReorderInts(h.DescID, order)
	h.UPID = ReorderInts(h.UPID, order)
	h.DFID = ReorderInts(h.DFID, order)
	h.Snap = ReorderInts(h.Snap, order)
	h.Mvir = ReorderFloat64s(h.Mvir, order)
	h.IDTable = NewLookupTable(h.ID)
}

type Tracks struct {
	N int
	Starts, Ends []int32
	IsReal, IsDisappear, IsMWSub []bool
	HostIdx, TrackIdx [][]int32
	IsReverseMerger, IsReverseSub, IsValidHost [][]bool
	Mpeak []float64
	MostMassivePreHostTrack []int32
}

func CalcTracks(h *Haloes) *Tracks {
	t := &Tracks{ }
	
	t.Starts, t.Ends = StartsEnds(h)
	t.N = len(t.Starts)
	
	t.IsReal = IsReal(h, t)
	t.IsDisappear = IsDisappear(h, t)
	t.IsMWSub = IsMWSub(h, t)
		
	t.HostIdx = FindAllHosts(h, t)	
	t.TrackIdx = FindTrackIndices(h, t)
	
	t.IsReverseMerger = IsReverseMerger(h, t)
	t.IsReverseSub = IsReverseSub(h, t)
	t.IsValidHost = IsValidHost(h, t)

	t.Mpeak = Mpeak(h, t)
	t.MostMassivePreHostTrack = MostMasssivePreHostTrack(h, t)

	return t
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

func IsReal(h *Haloes, t *Tracks) []bool {
	isReal := make([]bool, t.N)
	for i := range isReal {
		isReal[i] = h.UPID[t.Ends[i] - 1] == -1
	}
	return isReal
}

func IsDisappear(h *Haloes, t *Tracks) []bool {
	out := make([]bool, t.N)
	
	for i := range out {
		out[i] = h.DescID[t.Starts[i]] == -1 &&
			h.Snap[t.Starts[i]] != MaxSnap
	}

	return out
}

func IsMWSub(h *Haloes, t *Tracks) []bool {
	minID, maxID := h.DFID[t.Starts[0]], h.DFID[t.Ends[0] - 1]
	
	out := make([]bool, t.N)
	for i := range out {
		start, end := t.Starts[i], t.Ends[i]
		for j := start; j < end; j++ {
			if h.UPID[j] != -1 {
				k := h.IDTable.Find(h.UPID[j])
				out[i] = h.DFID[k] >= minID && h.DFID[k] <= maxID
			}
		}
	}

	return out
}

func FindAllHosts(h *Haloes, t *Tracks) [][]int32 {
	out := make([][]int32, t.N)

	for i := range t.Starts {
		for j := t.Starts[i]; j < t.Ends[i]; j++ {
			if h.UPID[j] != -1 {
				k := h.IDTable.Find(h.UPID[j])
				out[i] = append(out[i], k)
				if i == 1 {
					fmt.Println(j, h.UPID[j], k)
				}
			}
		}
	}

	return out
}

func FindTrackIndices(h *Haloes, t *Tracks) [][]int32 {
	out := make([][]int32, t.N)

	for i := range t.HostIdx {
		out[i] = make([]int32, len(t.HostIdx[i]), len(t.HostIdx[i]))

		for j := range out[i] {
			start := FindFirstIndex(h, t.HostIdx[i][j])
			out[i][j] =  int32(sort.Search(t.N, func (k int) bool {
				return t.Starts[k] >= start
			}))
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
			start := t.Starts[t.TrackIdx[i][j]]
			if h.DescID[start] != -1 {
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
			start := t.Starts[t.TrackIdx[i][j]]

			for k := start; k < t.HostIdx[i][j]; k++ {
				if h.UPID[k] != -1 {
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
		out[i] = make([]bool, len(t.TrackIdx[i]), len(t.TrackIdx))
		for j := range out[i] {
			k := t.TrackIdx[i][j]
			out[i][j] = t.IsReal[k] && !t.IsDisappear[k] &&
				!t.IsReverseMerger[i][j] && !t.IsReverseSub[i][j]
			
		}
	}

	return out
}

func Mpeak(h *Haloes, t *Tracks) []float64 {
	out := make([]float64, t.N)
	for i := range out {
		for j := t.Starts[i]; j <= t.Ends[j]; j++ {
			if h.Mvir[j] > out[i] { out[i] = h.Mvir[j] }
		}
	}

	return out
}

func MostMasssivePreHostTrack(h *Haloes, t *Tracks) []int32 {
	out := make([]int32, t.N)
	for i := range out {
		out[i] = -1
		for j := range t.TrackIdx[i] {
			if t.TrackIdx[i][j] == 0 { continue }
			k := t.TrackIdx[i][j]
			if out[i] == -1 || t.Mpeak[k] > t.Mpeak[out[i]] {
				out[i] = k
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
