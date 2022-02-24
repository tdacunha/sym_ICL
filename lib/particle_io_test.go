package lib

import (
	"testing"
)

func TestReadTags(t *testing.T) {
	tags := &Tags{
		Idx: [][]int32{
			{1, 2},
			{3},
			{4, 5, 6, 7, 8},
			{9},
		},
		Snap: [][]int16{
			{0, 1},
			{0},
			{0, 1, 2, 3, 4},
			{0},
		},
		Flag: [][]uint8{
			{0, 0},
			{1},
			{2, 2, 2, 2, 2},
			{3},
		},
	}

	WriteTags("test_dir", 3, tags)
	hd := ReadParticleHeader("test_dir")
	rTags := ReadTags("test_dir", hd)

	if !Int32Eq2D(tags.Idx, rTags.Idx) ||
		!Int16Eq2D(tags.Snap, rTags.Snap) ||
		!Uint8Eq2D(tags.Flag, rTags.Flag) {
		t.Errorf("Wrote tags: %v, but read tags %v", tags, rTags)
	}
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
