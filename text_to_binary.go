package main

import (
	"fmt"
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
	DFIDCol = 28
	// ID, DescID, UPID, Phantom, Snap, NextProg
	IntCols = []int{ 1, 3, 6, 8, 31, 32 }
	// Mvir, Rs, Vmax, M200b, M200c, M500c, Xoff, SpinBullock, b/a, c/a, T/U
	FloatCols = []int{ 10, 12, 16, 39, 40, 41, 43, 45, 46, 47, 56 }
	// X(3), V(3), J(3), A(3)
	VectorCols = []int{17, 18, 19, 20, 21, 22, 23, 24, 25, 48, 49, 50 }
)

// HaloTrack contains the evolution history of a single halo.
func main() {	
    if len(os.Args) != 2 {
        panic(fmt.Sprintf("You must supply a file with tree names and " +
            "output directories"))
    }

	inputName := os.Args[1]
	treeDirs, outDirs := ParseInputFile(inputName)

	log.Printf("Running on %d haloes", len(treeDirs))

	for i := range treeDirs {
		ConvertTree(treeDirs[i], outDirs[i])
	}
}

func MemoryLog() {
	ms := &runtime.MemStats{ }
	runtime.ReadMemStats(ms)
	log.Printf("Allocated: %.1f GB In Use: %.1f GB Idle: %.1f GB\n",
		float64(ms.Alloc) / 1e9, float64(ms.HeapInuse) / 1e9,
		float64(ms.HeapIdle) / 1e9)
}

func CountHeaderLines(fileName string) (nLines, nTrees int) {
	f, err := os.Open(fileName)
	if err != nil { panic(err.Error()) }
	defer f.Close()

	stat, err := f.Stat()
	if err != nil { panic(err.Error()) }
	nMax := stat.Size()

	if nMax > 50000 { nMax = 50000 }
	
	buf := make([]byte, nMax)
	_, err = io.ReadFull(f, buf)
	if err != nil { panic(err.Error()) }

	nLines = 0
	line := 1
	lineStart := 0
	nTrees = -1
	for i := range buf {
		switch buf[i] {
		case '\n', 0:
			line++

			if nTrees == -1 {
				lineEnd := i
				line := buf[lineStart: lineEnd]

				if line[0] != '#' {
					nTrees, err = strconv.Atoi(string(line))
					if err != nil { nTrees = -1 }
				}
				
				lineStart = lineEnd
			}

			if buf[i] == 0 { break }
			
		case '#':
			nLines = line
		}
	}

	return nLines, nTrees
}

func ConvertTree(treeDir, outDir string) {
	files := TreeFileNames(treeDir)
	binFiles := BinFileNames(outDir, files)
	for i := range files {
		runtime.GC()
		ConvertTreeFile(files[i], binFiles[i])
	}
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

func BinFileNames(outDir string, treeFileNames []string) []string {
	out := []string{ }
	for i := range treeFileNames {
		name := path.Base(treeFileNames[i])
		nameBase := strings.Split(name, ".")[0]
		outName := fmt.Sprintf("%s.df.bin", nameBase)
		out = append(out, path.Join(outDir, outName))
	}
	return out
}

func ConvertTreeFile(treeName, outName string) {
	cfg := catio.DefaultConfig

	var nTrees int
	if cfg.SkipLines, nTrees = CountHeaderLines(treeName); nTrees == 0 {
		return
	}

	log.Printf("Parsing %s", treeName)
	rd := catio.TextFile(treeName, cfg)

	f, err := os.Create(outName)
	if err != nil { panic(err.Error()) }
	defer f.Close()
	
	order := ConvertInts(rd, f)
	runtime.GC()
	ConvertFloats(rd, f, order)
	runtime.GC()
	ConvertVectors(rd, f, order)
	runtime.GC()
}

func WriteHeader(f *os.File, n int) {
	WriteInt(f, n)
	WriteInt(f, 1 + len(IntCols) + len(FloatCols) + len(VectorCols))
	WriteInt(f, DFIDCol)
	for i := range IntCols { WriteInt(f, IntCols[i]) }
	for i := range FloatCols { WriteInt(f, FloatCols[i]) }
	for i := range VectorCols { WriteInt(f, VectorCols[i]) }
}

func ConvertInts(rd catio.Reader, f *os.File) []int32 {
	// Read the columns
	colIdxs := []int{ DFIDCol }
	colIdxs = append(colIdxs, IntCols...)
	cols := rd.ReadInts(colIdxs)
	dfid := ToInt32(cols[0])
	order := IDOrder(dfid)

	hd := CreateTreeHeader(len(order))
	lib.WriteTreeHeader(f, hd)

	buf := dfid
	byteOrder := binary.LittleEndian
	for i := range cols {
		ReorderToInt32(cols[i], order, buf)
		err := binary.Write(f, byteOrder, buf)
		if err != nil { panic(err.Error()) }
	}

	return order
}

func CreateTreeHeader(n int) *lib.TreeHeader {
	cols := []int32{ int32(DFIDCol) }
	for i := range IntCols {
		cols = append(cols, int32(IntCols[i]))
	}
	for i := range FloatCols {
		cols = append(cols, int32(FloatCols[i]))
	}
	for i := 0; i < len(VectorCols); i += 3 {
		cols = append(cols, int32(VectorCols[i]))
	}
	return &lib.TreeHeader{
		lib.FixedHeader{
			int32(n),
			int32(1 + len(IntCols)),
			int32(len(FloatCols)),
			int32(len(VectorCols) / 3),
		}, cols,
	}
}

func ConvertFloats(rd catio.Reader, f *os.File, order []int32) {
	cols := rd.ReadFloat32s(FloatCols)

	buf := make([]float32, len(cols[0]))
	for i := range cols {
		ReorderToFloat32(cols[i], order, buf)
		err := binary.Write(f, binary.LittleEndian, buf)
		if err != nil { panic(err.Error()) }
	}
}

func ConvertVectors(rd catio.Reader, f *os.File, order []int32) {
	cols := rd.ReadFloat32s(VectorCols)

	buf := make([][3]float32, len(cols[0]))
	for i := 0; i < len(cols); i += 3 {
		ReorderToVector32([3][]float32{cols[i],cols[i+1],cols[i+2]}, order, buf)
		err := binary.Write(f, binary.LittleEndian, buf)
		if err != nil { panic(err.Error()) }
	}
}

func ToInt32(x []int) []int32 {
	out := make([]int32, len(x))
	for i := range x {
		out[i] = int32(x[i])
	}

	return out
}

func ReorderToInt32(x []int, order, out []int32) {
	for i := range out {
		out[i] = int32(x[order[i]])
	}
}

func ReorderToFloat32(x []float32, order []int32, out []float32){
	for i := range out {
		out[i] = x[order[i]]
	}
}

func ReorderToVector32(x [3][]float32, order []int32, out [][3]float32) {
	for i := range out {
		for dim := 0; dim < 3; dim++ {
			out[i][dim] = x[dim][order[i]]
		}
	}
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

func WriteInt(f *os.File, n int) {
	err := binary.Write(f, binary.LittleEndian, int32(n))
	if err != nil { panic(err.Error()) }
}

func ParseInputFile(fname string) (trees, outs []string) {
	b, err := ioutil.ReadFile(fname)
	if err != nil {
		panic(fmt.Sprintf("Could not open input file: %s", err.Error()))
	}
	s := string(b)

	trees, outs = []string{ }, []string{ }
	
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
		outs = append(outs, cols[1])
	}
	
	return trees, outs
}
