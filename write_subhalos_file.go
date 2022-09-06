package main

import (
	"runtime"
	"os"
	"encoding/binary"
	"strconv"
	"log"
	
	"github.com/phil-mansfield/symphony_pipeline/lib"
)

const (
	NPeakMin = 300
)

type LookupTable struct {
	Order, ID []int32
}

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
		
		log.Printf("Running halo %d (%s)", i, cfg.BaseDir[i])
		lib.MemoryUsage()
		if haloIndex != -1 && i != haloIndex { continue }

		branchesName := lib.BranchesFileName(cfg.BaseDir[i])
		fileNames := lib.TreeFileNames(cfg.BaseDir[i])
		
		files := make([]*os.File, len(fileNames))
		for j := range files {
			var err error
			files[j], err = os.Open(fileNames[j])
			if err != nil { panic(err.Error()) }
		}
		b := lib.ReadBranches(branchesName)

		n := lib.TreeNTot(files)
		mvir := make([]float32, n)
		lib.ReadTreeVarFullFloat(files, "Mvir", mvir)

		log.Printf("Mpeak cutoff: %g", cfg.Mp[i]*NPeakMin)
		
		mpeak := Mpeak(b, mvir)

		mIdx := FindLargestMergers(b, mpeak, float32(cfg.Mp[i]*NPeakMin))
		runtime.GC()

		log.Printf("%d haloes tracked.", len(mIdx))
		
		snaps := GetSnaps(files, b, n, mIdx)
		runtime.GC()

		mFile, err := os.Create(lib.MergerFileName(cfg.BaseDir[i]))
		if err != nil { panic(err.Error()) }

		maxSnap := cfg.MaxSnap[i]
		
		WriteHeader(maxSnap, mIdx, mFile)
		WriteFloat(files, b, n, mIdx, snaps, maxSnap, "Mvir", mFile)
		runtime.GC()
		WriteFloat(files, b, n, mIdx, snaps, maxSnap, "Vmax", mFile)
		runtime.GC()
		WriteFloat(files, b, n, mIdx, snaps, maxSnap, "RVmax", mFile)
		runtime.GC()
		WriteInt(files, b, n, mIdx, snaps, maxSnap, "ID", mFile)
		runtime.GC()
		WriteVector(files, b, n, mIdx, snaps, maxSnap, "X", mFile)
		runtime.GC()
		WriteVector(files, b, n, mIdx, snaps, maxSnap, "V", mFile)
		runtime.GC()

		mFile.Close()		
	}
	log.Println("Finishing write_subhalos_file")
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
	mIdx []int32, snaps [][]int32, maxSnap int32,
	varName string, out *os.File) {
	
	x := make([]float32, n)
	lib.ReadTreeVarFullFloat(files, varName, x)

	for k, i := range mIdx {
		xi := x[b.Starts[i]: b.Ends[i]]
		xFlat := make([]float32, maxSnap+1)
		for j := range xFlat { xFlat[j] = -1 }
		for j := range xi {
			xFlat[snaps[k][j]] = xi[j]
		}
		binary.Write(out, lib.ByteOrder, xFlat)
	}
}

func WriteInt(files []*os.File, b *lib.Branches, n int, 
	mIdx []int32, snaps [][]int32, maxSnap int32,
	varName string, out *os.File) {

	x := make([]int32, n)
	lib.ReadTreeVarFullInt(files, varName, x)

	for k, i := range mIdx {
		xi := x[b.Starts[i]: b.Ends[i]]
		xFlat := make([]int32, maxSnap+1)
		for j := range xFlat { xFlat[j] = -1 }
		for j := range xi {
			xFlat[snaps[k][j]] = xi[j]
		}
		binary.Write(out, lib.ByteOrder, xFlat)
	}
}

func WriteVector(files []*os.File, b *lib.Branches, n int, 
	mIdx []int32, snaps [][]int32, maxSnap int32,
	varName string, out *os.File) {
	x := make([][3]float32, n)
	lib.ReadTreeVarFullVector(files, varName, x)

	for k, i := range mIdx {
		xi := x[b.Starts[i]: b.Ends[i]]
		
		xFlat := make([][3]float32, maxSnap+1)
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

func FindLargestMergers(
	b *lib.Branches, mpeak []float32, mpeakMin float32,
) []int32 {
	
	mpeak64 := make([]float64, b.N)
	for i := range mpeak64 { mpeak64[i] = float64(mpeak[i]) }
	order := lib.QuickSortIndex(mpeak64)

	idx := []int32{ b.CentralIdx }
	for i := len(mpeak) - 1; i >= 0; i-- {
		if mpeak[order[i]] < mpeakMin { continue }
		if b.IsMWSub[order[i]] {
			idx = append(idx, int32(order[i]))
		}
	}

	return idx
}
