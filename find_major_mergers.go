package main

import (
	"fmt"
	"runtime"
	"os"
	"encoding/binary"
	"io/ioutil"
	"strings"
	"strconv"
	"log"
	"path"
	
	"github.com/phil-mansfield/lmc_ges_tracking/lib"
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
	treeDirs, mwIDs := ParseInputFile(inputName)

	log.Printf("Running on %d haloes", len(mwIDs))

	for i := range treeDirs {
		log.Println(treeDirs[i])
		MemoryLog()
		branchesName := lib.BranchesName(treeDirs[i])
		fileNames := lib.TreeFileNames(treeDirs[i])
		files := make([]*os.File, len(fileNames))
		for i := range files {
			var err error
			files[i], err = os.Open(fileNames[i])
			if err != nil { panic(err.Error()) }
		}
		b := lib.ReadBranches(branchesName)

		n := lib.TreeNTot(files)
		mvir := make([]float32, n)
		lib.ReadTreeVarFullFloat(files, "Mvir", mvir)
		
		mpeak := Mpeak(b, mvir)
		mIdx := FindLargestMergers(b, mpeak)
		runtime.GC()
		snaps := GetSnaps(files, b, n, mIdx)
		runtime.GC()

		mFile, err := os.Create(path.Join(treeDirs[i], "mergers.dat"))
		if err != nil { panic(err.Error()) }

		WriteHeader(MaxSnap, mIdx, mFile)
		WriteFloat(files, b, n, mIdx, snaps, "Mvir", mFile)
		runtime.GC()
		WriteFloat(files, b, n, mIdx, snaps, "Vmax", mFile)
		runtime.GC()
		WriteInt(files, b, n, mIdx, snaps, "ID", mFile)
		runtime.GC()
		WriteVector(files, b, n, mIdx, snaps, "X", mFile)
		runtime.GC()
		WriteVector(files, b, n, mIdx, snaps, "V", mFile)
		runtime.GC()

		mFile.Close()
	}
}

func WriteHeader(snaps int32, mIdx []int32, out *os.File) {
	err := binary.Write(out, lib.ByteOrder, snaps+1)
	if err != nil { panic(err.Error()) }
	err = binary.Write(out, lib.ByteOrder, int32(len(mIdx)))
	if err != nil { panic(err.Error()) }
	err = binary.Write(out, lib.ByteOrder, mIdx)
	if err != nil { panic(err.Error()) }
}

func WriteFloat(files []*os.File, b *lib.Branches, n int, 
	mIdx []int32, snaps [][]int32, varName string, out *os.File) {
	
	x := make([]float32, n)
	lib.ReadTreeVarFullFloat(files, varName, x)

	for k, i := range mIdx {
		xi := x[b.Starts[i]: b.Ends[i]]
		xFlat := make([]float32, MaxSnap+1)
		for j := range xFlat { xFlat[j] = -1 }
		for j := range xi {
			xFlat[snaps[k][j]] = xi[j]
		}
		binary.Write(out, lib.ByteOrder, xFlat)
	}
}

func WriteInt(files []*os.File, b *lib.Branches, n int, 
	mIdx []int32, snaps [][]int32, varName string, out *os.File) {

	x := make([]int32, n)
	lib.ReadTreeVarFullInt(files, varName, x)

	for k, i := range mIdx {
		xi := x[b.Starts[i]: b.Ends[i]]
		xFlat := make([]int32, MaxSnap+1)
		for j := range xFlat { xFlat[j] = -1 }
		for j := range xi {
			xFlat[snaps[k][j]] = xi[j]
		}
		binary.Write(out, lib.ByteOrder, xFlat)
	}
}

func WriteVector(files []*os.File, b *lib.Branches, n int, 
	mIdx []int32, snaps [][]int32, varName string, out *os.File) {
	x := make([][3]float32, n)
	lib.ReadTreeVarFullVector(files, varName, x)

	for k, i := range mIdx {
		xi := x[b.Starts[i]: b.Ends[i]]
		xFlat := make([][3]float32, MaxSnap+1)
		for j := range xFlat { xFlat[j] = [3]float32{ -1, -1, -1 } }
		for j := range xi {
			xFlat[snaps[k][j]] = xi[j]
		}
		binary.Write(out, lib.ByteOrder, xFlat)
	}
}

func GetSnaps(files []*os.File, b *lib.Branches,
	n int, mIdx []int32) [][]int32 {

	snap := make([]int32, n)
	snaps := make([][]int32, len(mIdx))
	lib.ReadTreeVarFullInt(files, "Snap", snap)
	for i, j := range mIdx {
		snapi := snap[b.Starts[j]: b.Ends[j]]
		snaps[i] = make([]int32, len(snapi))
		for k := range snapi { snaps[i][k] = snapi[k] }
	}
	return snaps
}

func Mpeak(b *lib.Branches, mvir []float32) []float32 {
	mpeak := make([]float32, b.N)
	for i := range b.Starts {
		if !b.IsReal[i] { continue }
		for j := b.Starts[i]; j < b.Ends[i]; j++ {
			if mvir[j] > mpeak[i] { mpeak[i] = mvir[j] }
		}
	}
	return mpeak
}

func FindLargestMergers(b *lib.Branches, mpeak []float32) []int32 {
	mpeak64 := make([]float64, b.N)
	for i := range mpeak64 { mpeak64[i] = float64(mpeak[i]) }
	order := lib.QuickSortIndex(mpeak64)
	
	idx := []int32{ b.CentralIdx }
	for i := len(mpeak) - 1; i >= 0; i-- {
		if b.IsMWSub[order[i]] {
			idx = append(idx, int32(order[i]))
			if len(idx) == NMerger + 1 { break }
		}
	}

	return idx
}

func MemoryLog() {
	ms := &runtime.MemStats{ }
	runtime.ReadMemStats(ms)
	log.Printf("Allocated: %.1f In Use: %.1f Idle: %.1f\n",
		float64(ms.Alloc) / 1e9, float64(ms.HeapInuse) / 1e9,
		float64(ms.HeapIdle) / 1e9)
}

func ParseInputFile(fname string) ([]string, []int32) {
	b, err := ioutil.ReadFile(fname)
	if err != nil {
		panic(fmt.Sprintf("Could not open input file: %s", err.Error()))
	}
	s := string(b)

	trees, mwIDs := []string{ }, []int32{ }
	
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
		
		if len(cols) != 2 {
			panic(fmt.Sprintf("Line %d of %s is '%s', but you need there " +
				"to be three columns.", i+1, fname, line))
		}
		
		trees = append(trees, cols[0])
		mwID, err := strconv.Atoi(cols[1])
		if err != nil {
			panic(fmt.Sprintf("Could not parse the ID on line %d of " +
				"%s: %s", i+1, fname, cols[1]))
		}
		mwIDs = append(mwIDs, int32(mwID))
	}
	
	return trees, mwIDs
}
