package lib

import (
	"encoding/binary"
	"fmt"
	"os"
)

var (
	ByteOrder = binary.LittleEndian
	ColLookup = map[string]int32{
		"DFID": 28,
		"ID": 1,
		"DescID": 3,
		"UPID": 6,
		"Phantom": 8,
		"Snap": 31,
		"NextProg": 32,
		"Mvir": 10,
		"Rs": 12,
		"Vmax": 16,
		"M200b": 39,
		"M200c": 40,
		"M500c": 41,
		"Xoff": 43,
		"SpinBullock": 45,
		"BToA": 46,
		"CToA": 47,
		"VirialRatio": 56,
		"RVmax": 60,
		"X": 17,
		"V": 20,
		"J": 23,
		"A": 48,
	}
)

type FixedHeader struct {
	N, NInt, NFloat, NVector int32
}

type TreeHeader struct {
	FixedHeader
	Cols []int32
}

func WriteTreeHeader(f *os.File, hd *TreeHeader) {
	err := binary.Write(f, ByteOrder, hd.FixedHeader)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, ByteOrder, hd.Cols)
	if err != nil { panic(err.Error()) }
}

func ReadTreeHeader(f *os.File) *TreeHeader {
	hd := &TreeHeader{ }
	_, err := f.Seek(0, 0)
	if err != nil { panic(err.Error()) }
	
	err = binary.Read(f, ByteOrder, &hd.FixedHeader)
	if err != nil { panic(err.Error()) }
	hd.Cols = make([]int32, hd.NInt + hd.NFloat + hd.NVector)
	err = binary.Read(f, ByteOrder, hd.Cols)
	if err != nil { panic(err.Error()) }

	return hd
}

func ReadTreeVar(f *os.File, varName string, buf interface{}) {
	hd := ReadTreeHeader(f)

	varCol, ok := ColLookup[varName]
	if !ok {
		panic(fmt.Sprintf("Unrecognized variable name, \"%s\"", varName))
	}

	varColIdx := -1
	for i := range hd.Cols {
		if hd.Cols[i] == varCol {
			varColIdx = i
			break
		}
	}
	if varColIdx == -1 {
		panic(fmt.Sprintf("\"%s\" is a valid variable name, but it " +
			"isn't stored in this tree file.", varName))
	}

	offset := 4*int64(hd.N)*int64(varColIdx)
	if varColIdx > int(hd.NInt + hd.NFloat) {
		offset += 8*int64(hd.N)*int64(varColIdx - int(hd.NInt + hd.NFloat))
	}

	_, err := f.Seek(offset, 1)
	if err != nil { panic(err.Error()) }
	err = binary.Read(f, ByteOrder, buf)
	if err != nil { panic(err.Error()) }
}

func TreeNTot(files []*os.File) int {
	nTot := 0
	for i := range files {
		hd := ReadTreeHeader(files[i])
		nTot += int(hd.N)
	}
	return nTot
}

func ReadTreeVarFullInt(files []*os.File, varName string, buf []int32) {
	start := 0
	for i := range files {
		hd := ReadTreeHeader(files[i])
		end := start + int(hd.N)
		subBuf := buf[start: end]
		ReadTreeVar(files[i], varName, subBuf)
		start = end
	}
}

func ReadTreeVarFullFloat(files []*os.File, varName string, buf []float32) {
	start := 0
	for i := range files {
		hd := ReadTreeHeader(files[i])
		end := start + int(hd.N)
		subBuf := buf[start: end]
		ReadTreeVar(files[i], varName, subBuf)
		start = end
	}
}

func ReadTreeVarFullVector(files []*os.File, varName string, buf [][3]float32) {
	start := 0
	for i := range files {
		hd := ReadTreeHeader(files[i])
		end := start + int(hd.N)
		subBuf := buf[start: end]
		ReadTreeVar(files[i], varName, subBuf)
		start = end
	}
}

type Branches struct {
	N, CentralIdx int32
	Starts, Ends []int32
	IsReal, IsDisappear, IsMWSub []bool
	PreprocessIdx []int32
}

func ReadBranches(name string) *Branches {
	f, err := os.Open(name)
	if err != nil { panic(err.Error()) }
	defer f.Close()

	b := &Branches{ }

	err = binary.Read(f, ByteOrder, &b.N)
	if err != nil { panic(err.Error()) }
	err = binary.Read(f, ByteOrder, &b.CentralIdx)
	if err != nil { panic(err.Error()) }

	edges := make([]int32, b.N+1)
	b.IsReal = make([]bool, b.N)
	b.IsDisappear = make([]bool, b.N)
	b.IsMWSub = make([]bool, b.N)
	b.PreprocessIdx = make([]int32, b.N)
	
	err = binary.Read(f, ByteOrder, edges)
	if err != nil { panic(err.Error()) }
	b.Starts = edges[:len(edges) - 1]
	b.Ends = edges[1:]
	err = binary.Read(f, ByteOrder, b.IsReal)
	if err != nil { panic(err.Error()) }
	err = binary.Read(f, ByteOrder, b.IsDisappear)
	if err != nil { panic(err.Error()) }
	err = binary.Read(f, ByteOrder, b.IsMWSub)
	if err != nil { panic(err.Error()) }
	err = binary.Read(f, ByteOrder, b.PreprocessIdx)
	if err != nil { panic(err.Error()) }
	
	return b
}
