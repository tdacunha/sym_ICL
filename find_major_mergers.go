package main

import (
	"fmt"
	"sort"
	"runtime"
	"io"
	"os"
	"path"
	"encoding/binary"
	"io/ioutil"
	"strings"
	"strconv"
	"log"
	
	"github.com/phil-mansfield/lmc_ges_tracking/lib"
	"github.com/phil-mansfield/guppy/lib/catio"
)

var (
	MaxSnap = int32(235)
	NMerger = 10
)

// Haloes is a collection of raw columns from the tree.dat files.
type Haloes struct {
	ID, DescID, UPID, DFID, Snap []int32
	Mvir, Vmax []float32
	X, V [][3]float32
	IDTable *LookupTable
}

type LookupTable struct {
	Order, ID []int32
}


// HaloTrack contains the evolution history of a single halo.
func main() {	
    if len(os.Args) != 2 {
        panic(fmt.Sprintf("You must supply a file with tree names, " +
            "output names, and z=0 MW IDs"))
    }

	inputName := os.Args[1]
	treeDirs, outFileNames, mwIDs := ParseInputFile(inputName)

	log.Printf("Running on %d haloes", len(mwIDs))

	for i := range treeDirs {
		h := ReadTree(treeDirs[i])
		MemoryLog()

		log.Printf("Analyzing halo-%d's mergers", mwIDs[i])
		runtime.GC()

		SortHaloes(h)
		runtime.GC()

		t := CalcTracks(h, mwIDs[i])
		mergers := MajorMergers(t, NMerger)
		MemoryLog()
		log.Printf("Writing %d mergers to %s", NMerger, outFileNames[i])
		WriteMergers(outFileNames[i], h, t, mergers)
		runtime.GC()
	}
}

func MemoryLog() {
	ms := &runtime.MemStats{ }
	runtime.ReadMemStats(ms)
	log.Printf("Allocated: %.1f In Use: %.1f Idle: %.1f\n",
		float64(ms.Alloc) / 1e9, float64(ms.HeapInuse) / 1e9,
		float64(ms.HeapIdle) / 1e9)
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

func CountHeaderLines(fileName string) (n int , empty bool) {
	f, err := os.Open(fileName)
	if err != nil { panic(err.Error()) }
	defer f.Close()

	stat, err := f.Stat()
	if err != nil { panic(err.Error()) }
	nMax := stat.Size()

	if nMax > 50000 { nMax = 50000 }
	
	buf := make([]byte, nMax)
	var nRead int
	nRead, err = io.ReadFull(f, buf)
	if err != nil { panic(err.Error()) }
	if int64(nRead) < nMax { return 0, true }

	headerLines := 0
	line := 1
	for i := range buf {
		switch buf[i] {
		case '\n':
			line++
		case '#':
			headerLines = line
		}
	}

	return headerLines, false
}

func ReadTree(treeDir string) *Haloes {
	files := TreeFileNames(treeDir)
	hs := []*Haloes{ }
	for i := range files {
		runtime.GC()
		h := ReadTreeFile(files[i])
		if h != nil { hs = append(hs, h) }
	}

	return JoinHaloes(hs)
}

func TreeFileNames(dir string) []string {
	files, err := ioutil.ReadDir(dir)
	if err != nil { panic(err.Error()) }
	
	out := []string{ }
	for i := range files {
		name := files[i].Name()
		if len(name) >= 8 &&
			name[:4] == "tree" &&
			name[len(name)-4:] == ".dat" {
			out = append(out, path.Join(dir, name))
		}
	}

	return out
}

func ReadTreeFile(file string) *Haloes {
	cfg := catio.DefaultConfig

	var empty bool
	if cfg.SkipLines, empty = CountHeaderLines(file); empty {
		return nil
	}


	log.Printf("Parsing %s", file)
	rd := catio.TextFile(file, cfg)
	
	// Read the columns
	cols := rd.ReadInts([]int{ 1, 3, 6, 28, 31 })
	id, descid, upid, dfid, snap := cols[0], cols[1],
		cols[2], cols[3],cols[4]
	fcols := rd.ReadFloat32s([]int{ 10, 16, 17, 18, 19, 20, 21, 22})
	mvir, vmax, x, y, z, vx, vy, vz := fcols[0], fcols[1],
		fcols[2], fcols[3], fcols[4], fcols[5], fcols[6], fcols[7]
	MemoryLog()
	return &Haloes{ ID: ToInt32(id), DescID: ToInt32(descid),
		UPID: ToInt32(upid), DFID: ToInt32(dfid), Snap: ToInt32(snap),
		Mvir: mvir, Vmax: vmax, X: ToVector(x, y, z), V: ToVector(vx, vy, vz) }
}

func JoinHaloes(hs []*Haloes) *Haloes {
	h := &Haloes{ }
	log.Println("IDs")
	MemoryLog()
	JoinID(hs, h)
	runtime.GC()
	log.Println("DescID")
	MemoryLog()
	JoinDescID(hs, h)
	runtime.GC()
	log.Println("UPID")
	MemoryLog()
	JoinUPID(hs, h)
	runtime.GC()
	log.Println("DFID")
	MemoryLog()
	JoinDFID(hs, h)
	runtime.GC()
	log.Println("Snap")
	MemoryLog()
	JoinSnap(hs, h)
	runtime.GC()
	log.Println("Mvir")
	MemoryLog()
	JoinMvir(hs, h)
	runtime.GC()
	log.Println("VMax")
	MemoryLog()
	JoinVmax(hs, h)
	runtime.GC()
	log.Println("X")
	MemoryLog()
	JoinX(hs, h)
	runtime.GC()
	log.Println("V")
	MemoryLog()
	JoinV(hs, h)
	return h
}

func JoinID(hs []*Haloes, h *Haloes) {
	xs := [][]int32{ }
	for i := range hs { xs = append(xs, hs[i].ID) }
	h.ID = JoinInt32(xs)
	for i := range hs { hs[i].ID = nil }
}

func JoinDescID(hs []*Haloes, h *Haloes) {
	xs := [][]int32{ }
	for i := range hs { xs = append(xs, hs[i].DescID) }
	h.DescID = JoinInt32(xs)
	for i := range hs { hs[i].DescID = nil }
}

func JoinUPID(hs []*Haloes, h *Haloes) {
	xs := [][]int32{ }
	for i := range hs { xs = append(xs, hs[i].UPID) }
	h.UPID = JoinInt32(xs)
	for i := range hs { hs[i].UPID = nil }
}

func JoinDFID(hs []*Haloes, h *Haloes) {
	xs := [][]int32{ }
	for i := range hs { xs = append(xs, hs[i].DFID) }
	h.DFID = JoinInt32(xs)
	for i := range hs { hs[i].DFID = nil }
}

func JoinSnap(hs []*Haloes, h *Haloes) {
	xs := [][]int32{ }
	for i := range hs { xs = append(xs, hs[i].Snap) }
	h.Snap = JoinInt32(xs)
	for i := range hs { hs[i].Snap = nil }
}

func JoinMvir(hs []*Haloes, h *Haloes) {
	xs := [][]float32{ }
	for i := range hs { xs = append(xs, hs[i].Mvir) }
	h.Mvir = JoinFloat32(xs)
	for i := range hs { hs[i].Mvir = nil }
}

func JoinVmax(hs []*Haloes, h *Haloes) {
	xs := [][]float32{ }
	for i := range hs { xs = append(xs, hs[i].Vmax) }
	h.Vmax = JoinFloat32(xs)
	for i := range hs { hs[i].Vmax = nil }
}

func JoinX(hs []*Haloes, h *Haloes) {
	xs := [][][3]float32{ }
	for i := range hs { xs = append(xs, hs[i].X) }
	h.X = JoinVec32(xs)
	for i := range hs { hs[i].X = nil }
}

func JoinV(hs []*Haloes, h *Haloes) {
	xs := [][][3]float32{ }
	for i := range hs { xs = append(xs, hs[i].V) }
	h.V = JoinVec32(xs)
	for i := range hs { hs[i].V = nil }
}


func JoinInt32(xs [][]int32) []int32 {
	n := 0
	for i := range xs { n += len(xs[i]) }
	
	out := make([]int32, n)
	j := 0
	for i := range xs {
		for k := range xs[i] {
			out[j] = xs[i][k]
			j++
		}
	}

	return out
}

func JoinFloat32(xs [][]float32) []float32 {
	n := 0
	for i := range xs { n += len(xs[i]) }
	
	out := make([]float32, n)
	j := 0
	for i := range xs {
		for k := range xs[i] {
			out[j] = xs[i][k]
			j++
		}
	}

	return out
}

func JoinVec32(xs [][][3]float32) [][3]float32 {
	n := 0
	for i := range xs { n += len(xs[i]) }
	
	out := make([][3]float32, n)
	j := 0
	for i := range xs {
		for k := range xs[i] {
			out[j] = xs[i][k]
			j++
		}
	}

	return out
}

func ToInt32(x []int) []int32 {
	out := make([]int32, len(x))
	for i := range x {
		out[i] = int32(x[i])
	}

	return out
}

func ToVector(x, y, z []float32) [][3]float32 {
	out := make([][3]float32, len(x))
	for i := range out {
		out[i] = [3]float32{ x[i], y[i], z[i] }
	}
	return out
}

func ReorderInts(x, order []int32) []int32 {
	out := make([]int32, len(order))
	for i := range out {
		out[i] = x[order[i]]
	}
	return out
}

func ReorderFloats(x []float32, order []int32) []float32 {
	out := make([]float32, len(order))
	for i := range out {
		out[i] = x[order[i]]
	}
	return out
}

func ReorderVectors(x [][3]float32, order []int32) [][3]float32 {
	out := make([][3]float32, len(order))
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
	h.Mvir = ReorderFloats(h.Mvir, order)
	h.Vmax = ReorderFloats(h.Vmax, order)
	h.X = ReorderVectors(h.X, order)
	h.V = ReorderVectors(h.V, order)
	h.IDTable = NewLookupTable(h.ID)
}

type Tracks struct {
	N, MWIdx int
	Starts, Ends []int32
	IsReal, IsMWSub []bool
	MpeakInfall []float32
}

func CalcTracks(h *Haloes, mwID int32) *Tracks {
	t := &Tracks{ }
	
	t.Starts, t.Ends = StartsEnds(h)
	t.N = len(t.Starts)
	
	t.IsReal = IsReal(h, t)
	t.IsMWSub = IsMWSub(h, t, mwID)
	t.MpeakInfall = MpeakInfall(h, t)
	
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

func IsMWSub(h *Haloes, t *Tracks, mwID int32) []bool {
	mwIdx := 0
	for ; mwIdx < t.N && h.ID[t.Starts[mwIdx]] != mwID; mwIdx++ {
	}

	if mwIdx == t.N {
		panic(fmt.Sprintf("Could not find a z=0 halo with ID %d", mwID))
	}
	
	t.MWIdx = mwIdx
	minID, maxID := h.DFID[t.Starts[mwIdx]], h.DFID[t.Ends[mwIdx] - 1]
	
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

func MpeakInfall(h *Haloes, t *Tracks) []float32 {
	out := make([]float32, t.N)
	for i := range out {
		for j := t.Ends[i]-1; j >= t.Starts[i] && h.UPID[j] == -1; j-- {
			if h.Mvir[j] > out[i] { out[i] = h.Mvir[j] }
		}
	}

	return out
}

func MajorMergers(t *Tracks, n int) []int {
	idx := []int{ }
	mpeak := []float64{ }

	for i := 0; i < t.N; i++ {
		if t.IsMWSub[i] {
			idx = append(idx, i)
			mpeak = append(mpeak, float64(t.MpeakInfall[i]))
		}
	}

	order := lib.QuickSortIndex(mpeak)

	mergers := []int{ }
	for i := len(order) - 1; i > len(order) - 1 - n; i-- {
		mergers = append(mergers, idx[order[i]])
	}

	return mergers
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

func WriteMergers(fname string, h *Haloes, t *Tracks, mergers []int) {
	f, err := os.Create(fname)
	if err != nil { panic(err.Error()) }
	defer f.Close()

	order := binary.LittleEndian
	err = binary.Write(f, order, int64(len(mergers)))
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, int64(MaxSnap))
	if err != nil { panic(err.Error()) }
	
	WriteHalo(f, h, t, t.MWIdx)
	for i := range mergers {
		WriteHalo(f, h, t, mergers[i])
	}
}

func WriteHalo(f *os.File, h *Haloes, t *Tracks, i0 int) {
	id := make([]int32, MaxSnap+1)
	mvir := make([]float32, MaxSnap+1)
	vmax := make([]float32, MaxSnap+1)
	x := make([][3]float32, MaxSnap+1)
	v := make([][3]float32, MaxSnap+1)

	for i := t.Starts[i0]; i < t.Ends[i0]; i++ {
		s := h.Snap[i]
		id[s] = h.ID[i]
		mvir[s] = h.Mvir[i]
		vmax[s] = h.Vmax[i]
		x[s] = h.X[i]
		v[s] = h.V[i]
	}

	order := binary.LittleEndian
	err := binary.Write(f, order, id)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, mvir)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, vmax)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, x)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, v)
	if err != nil { panic(err.Error()) }
}

func ParseInputFile(fname string) ([]string, []string, []int32) {
	b, err := ioutil.ReadFile(fname)
	if err != nil {
		panic(fmt.Sprintf("Could not open input file: %s", err.Error()))
	}
	s := string(b)

	trees, outs, mwIDs := []string{ }, []string{ }, []int32{ }
	
	lines := strings.Split(s, "\n")
	for i := range lines {
		line := strings.Trim(lines[i], " ")
		if len(line) == 0 { continue }

		tok := strings.Split(lines[i], " ")
		cols := []string{ }
		for i := range tok {
			if len(tok[i]) > 0 {
				cols = append(cols, tok[i])
			}
		}
		
		if len(cols) != 3 {
			panic(fmt.Sprintf("Line %d of %s is '%s', but you need there " +
				"to be three columns.", i+1, fname, line))
		}
		
		trees = append(trees, cols[0])
		outs = append(outs, cols[1])
		mwID, err := strconv.Atoi(cols[2])
		if err != nil {
			panic(fmt.Sprintf("Could not parse the ID on line %d of " +
				"%s: %s", i+1, fname, cols[2]))
		}
		mwIDs = append(mwIDs, int32(mwID))
	}
	
	return trees, outs, mwIDs
}
