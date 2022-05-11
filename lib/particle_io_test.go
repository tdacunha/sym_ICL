package lib

import (
	"testing"
)

func TestReadTags(t *testing.T) {
	tags := &Tags{
		N0: []int32{ 1, 1, 0, 3, 0 },
		ID: [][]int32{
			{1, 2},
			{3},
			{ },
			{4, 5, 6, 7, 8},
			{9},
		},
		Snap: [][]int16{
			{0, 1},
			{0},
			{ },
			{0, 1, 2, 3, 4},
			{0},
		},
		Flag: [][]uint8{
			{0, 0},
			{1},
			{ },
			{2, 2, 2, 2, 2},
			{3},
		},
	}

	lookup := &TagLookup{
		Halo: make([]int16, 100),
		Index: make([]int32, 100),
	}

	for i := range lookup.Halo {
		lookup.Halo[i] = int16(i % 5)
		lookup.Index[i] = int32(i % 21)
	}

	WriteTags("test_dir", 3, tags, lookup)
	hd := ReadParticleHeader("test_dir")
	rTags := ReadTags("test_dir", hd)

	if !Int32Eq(tags.N0, rTags.N0) ||
		!Int32Eq2D(tags.ID, rTags.ID) ||
		!Int16Eq2D(tags.Snap, rTags.Snap) ||
		!Uint8Eq2D(tags.Flag, rTags.Flag) {
		t.Errorf("Wrote tags: %v, but read tags %v", tags, rTags)
	}

	rLookup := ReadTagLookup("test_dir")
	if !Int32Eq(lookup.Index, rLookup.Index) ||
		!Int16Eq(lookup.Halo, rLookup.Halo) {
		t.Errorf("Wrote lookup: %v, but read lookup %v", lookup, rLookup)
	}
}

func TestReadFloat(t *testing.T) {
	snap := 999
	
	x := [][]float32{
		{1.0, 2.0},
		{3.0},
		{ },
		{4.0, 5.0, 6.0, 7.0, 8.0},
		{9.0},
	}

	WriteFloat("test_dir", 3, "test_float", snap, x)
	hd := ReadParticleHeader("test_dir")
	rx:= ReadFloat("test_dir", "test_float", snap, hd)

	if !Float32Eq2D(x, rx, 1e-3) {
		t.Errorf("Wrote floats %v, but read floats %v", x, rx)
	}
}

func TestReadVector(t *testing.T) {
	snap := 999
	
	x := [][][3]float32{
		[][3]float32{{1,1,1}, {2,2,2}},
		[][3]float32{{3,3,3}},
		[][3]float32{ },
		[][3]float32{{4,4,4}, {5,5,5}, {6,6,6}, {7,7,7}, {8,8,8}},
		[][3]float32{{9,9,9}},
	}

	WriteVector("test_dir", 3, "test_vector", snap, x)
	hd := ReadParticleHeader("test_dir")
	rx:= ReadVector("test_dir", "test_vector", snap, hd)

	if !Vec32Eq2D(x, rx, 1e-3) {
		t.Errorf("Wrote vectors %v, but read vectors %v", x, rx)
	}
}

func Float32Eq2D(x, y [][]float32, eps float32) bool {
	if len(x) != len(y) { return false }
	for i := range x {
		if !Float32Eq(x[i], y[i], eps) { return false }
	}
	return true
}

func Vec32Eq2D(x, y [][][3]float32, eps float32) bool {
	if len(x) != len(y) { return false }
	for i := range x {
		if !Vec32Eq(x[i], y[i], eps) { return false }
	}
	return true
}


func Int32Eq2D(x, y [][]int32) bool {
	if len(x) != len(y) { return false }
	for i := range x {
		if !Int32Eq(x[i], y[i]) { return false }
	}
	return true
}

func Int16Eq2D(x, y [][]int16) bool {
	if len(x) != len(y) { return false }
	for i := range x {
		if !Int16Eq(x[i], y[i]) { return false }
	}
	return true
}

func Uint8Eq2D(x, y [][]uint8) bool {
	if len(x) != len(y) { return false }
	for i := range x {
		if !Uint8Eq(x[i], y[i]) { return false }
	}
	return true
}

func EpsEq(x, y, eps float32) bool {
	return x + eps >= y && x - eps <= y
}

func Float32Eq(x, y []float32, eps float32) bool {
	if len(x) != len(y) { return false }
	for i := range x {
		if !EpsEq(x[i], y[i], eps) { return false }
	}
	return true
}

func Vec32Eq(x, y [][3]float32, eps float32) bool {
	if len(x) != len(y) { return false }
	for i := range x {
		for k := 0; k < 3; k++ {
			if !EpsEq(x[i][k], y[i][k], eps) { return false }
		}
	}
	return true
}

func Int16Eq(x, y []int16) bool {
	if len(x) != len(y) { return false }
	for i := range x {
		if x[i] != y[i] { return false }
	}
	return true
}

func Uint8Eq(x, y []uint8) bool {
	if len(x) != len(y) { return false }
	for i := range x {
		if x[i] != y[i] { return false }
	}
	return true
}
