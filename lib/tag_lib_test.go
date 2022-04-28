package lib

import (
	"testing"
	"math"
)

func TestTrimFloat(t *testing.T) {
	tags := &Tags{
		Snap: [][]int16{
			{},
			{5, 5, 5, 5, 5},
			{1, 3, 4, 5, 6},
			{10, 1, 10, 1, 1, 10, 1, 10},
		},
	}

	tests := []struct{
		iHalo, snap int
		in, out []float32
	} {
		{0, 1, []float32{}, []float32{}},
			{1, 5, []float32{1, 2, 3, 4, 5}, []float32{1, 2, 3, 4, 5}},
			{1, 1, []float32{1, 2, 3, 4, 5}, []float32{}},
			{2, 4, []float32{0, 1, 2, 3, 4}, []float32{0, 1, 2}},
			{3, 5, []float32{8, 7, 6, 5, 4, 3, 2, 1}, []float32{7, 5, 4, 2}},
	}

	for i := range tests {

		out := make([]float32, len(tests[i].in))
		copy(out, tests[i].in)
		out = tags.TrimFloat(tests[i].iHalo, out, tests[i].snap)

		if !Float32Eq(out, tests[i].out, 0.0) {
			t.Errorf("%d) Expected %.1f (halo = %d, snap = %d) to trim to " + 
				"%.1f, got %.1f", i+1, tests[i].in, tests[i].iHalo,
				tests[i].snap, tests[i].out, out)
		}
	}
}

func TestExpandFloat(t *testing.T) {
	tags := &Tags{
		Snap: [][]int16{
			{},
			{5, 5, 5, 5, 5},
			{1, 3, 4, 5, 6},
			{10, 1, 10, 1, 1, 10, 1, 10},
		},
	}

	fill := float32(math.Inf(0))

	tests := []struct{
		iHalo, snap int
		out, in []float32
	} {
		{0, 1, []float32{}, []float32{}},
			{1, 5, []float32{1, 2, 3, 4, 5}, []float32{1, 2, 3, 4, 5}},
			{1, 1, []float32{fill, fill, fill, fill, fill}, []float32{}},
			{2, 4, []float32{0, 1, 2, fill, fill}, []float32{0, 1, 2}},
			{3, 5, []float32{fill, 7, fill, 5, 4, fill, 2, fill}, 
				[]float32{7, 5, 4, 2}},
	}

	for i := range tests {

		out := make([]float32, len(tests[i].in))
		copy(out, tests[i].in)
		out = tags.ExpandFloat(tests[i].iHalo, out, tests[i].snap, fill)

		if !Float32Eq(out, tests[i].out, 0.0) {
			t.Errorf("%d) Expected %.1f (halo = %d, snap = %d) to expand to " + 
				"%.1f, got %.1f", i+1, tests[i].in, tests[i].iHalo,
				tests[i].snap, tests[i].out, out)
		}
	}
}

func TestTrimVector(t *testing.T) {
	tags := &Tags{
		Snap: [][]int16{
			{},
			{5, 5, 5, 5, 5},
			{1, 3, 4, 5, 6},
			{10, 1, 10, 1, 1, 10, 1, 10},
		},
	}

	tests := []struct{
		iHalo, snap int
		in, out [][3]float32
	} {
		{0, 1, [][3]float32{}, [][3]float32{}},
			{1, 5, [][3]float32{{1,1,1}, {2,2,2}, {3,3,3}, {4,4,4}, {5,5,5}}, [][3]float32{{1,1,1}, {2,2,2}, {3,3,3}, {4,4,4}, {5,5,5}}},
			{1, 1, [][3]float32{{1,1,1}, {2,2,2}, {3,3,3}, {4,4,4}, {5,5,5}}, [][3]float32{}},
			{2, 4, [][3]float32{{0,0,0}, {1,10,100}, {2,20,200}, {3,30,300}, {4,40,400}}, [][3]float32{{0,0,0}, {1,10,100}, {2,20,200}}},
			{3, 5, [][3]float32{{8,8,8}, {7,7,7}, {6,6,6}, {5,5,5}, {4,4,4}, {3,3,3}, {2,2,2}, {1,1,1}}, [][3]float32{{7,7,7}, {5,5,5}, {4,4,4}, {2,2,2}}},
	}

	for i := range tests {

		out := make([][3]float32, len(tests[i].in))
		copy(out, tests[i].in)
		out = tags.TrimVector(tests[i].iHalo, out, tests[i].snap)

		if !Vec32Eq(out, tests[i].out, 0.0) {
			t.Errorf("%d) Expected %.1f (halo = %d, snap = %d) to trim to " + 
				"%.1f, got %.1f", i+1, tests[i].in, tests[i].iHalo,
				tests[i].snap, tests[i].out, out)
		}
	}
}

func TestExpandVector(t *testing.T) {
	tags := &Tags{
		Snap: [][]int16{
			{},
			{5, 5, 5, 5, 5},
			{1, 3, 4, 5, 6},
			{10, 1, 10, 1, 1, 10, 1, 10},
		},
	}

	fill := [3]float32{  float32(math.Inf(0)),
		float32(math.Inf(0)), float32(math.Inf(0)) }

	tests := []struct{
		iHalo, snap int
		out, in [][3]float32
	} {
		{0, 1, [][3]float32{}, [][3]float32{}},
			{1, 0, [][3]float32{fill, fill, fill, fill, fill}, [][3]float32{{1,1,1}, {2,2,2}, {3,3,3}, {4,4,4}, {5,5,5}}},
			{1, 5, [][3]float32{{1,1,1}, {2,2,2}, {3,3,3}, {4,4,4}, {5,5,5}}, [][3]float32{{1,1,1}, {2,2,2}, {3,3,3}, {4,4,4}, {5,5,5}}},
			{2, 4, [][3]float32{{0,0,0}, {1,10,100}, {2,20,200}, fill, fill}, [][3]float32{{0,0,0}, {1,10,100}, {2,20,200}}},
			{3, 5, [][3]float32{fill, {7,7,7}, fill, {5,5,5}, {4,4,4}, fill, {2,2,2}, fill}, [][3]float32{{7,7,7}, {5,5,5}, {4,4,4}, {2,2,2}}},
			}

	for i := range tests {

		out := make([][3]float32, len(tests[i].in))
		copy(out, tests[i].in)
		out = tags.ExpandVector(tests[i].iHalo, out, tests[i].snap, fill)

		if !Vec32Eq(out, tests[i].out, 0.0) {
			t.Errorf("%d) Expected %.1f (halo = %d, snap = %d) to expand to " + 
				"%.1f, got %.1f", i+1, tests[i].in, tests[i].iHalo,
				tests[i].snap, tests[i].out, out)
		}
	}
}


func TestOrderTags(t *testing.T) {
	tags := &Tags{
		ID: [][]int32{
			{},
			{10, 20, 30, 40, 50},
			{10, 20, 30, 40, 50},
			{10, 20, 30, 40, 50, 60, 70, 80},
		},
		Snap: [][]int16{
			{},
			{1, 2, 3, 4, 5},
			{1, 2, 3, 4, 5},
			{1, 2, 3, 4, 5, 6, 7, 8},
		},
		Flag: [][]uint8{
			{},
			{0, 0, 0, 1, 1},
			{1, 1, 1, 0, 0},
			{1, 1, 0, 1, 0, 0, 1, 0},
		},
	}
	exp := &Tags{
		N0: []int32{ 0, 3, 2, 4 },
		ID: [][]int32{
			{},
			{10, 20, 30, 40, 50},
			{40, 50, 30, 10, 20},
			{50, 60, 30, 80, 10, 20, 70, 40},
		},
		Snap: [][]int16{
			{},
			{1, 2, 3, 4, 5},
			{4, 5, 3, 1, 2},
			{5, 6, 3, 8, 1, 2, 7, 4},
		},
		Flag: [][]uint8{
			{},
			{0, 0, 0, 1, 1},
			{0, 0, 1, 1, 1},
			{0, 0, 0, 0, 1, 1, 1, 1},
		},
	}

	tags.OrderTags()
	
	if !Int32Eq(tags.N0, exp.N0) {
		t.Errorf("Expected reordered tags.N0 = %d, got %d", exp.N0, tags.N0)
	}

	for i := range tags.N0 {
		if !Int32Eq(tags.ID[i], exp.ID[i]) {
			t.Errorf("Expected tags.ID[%d] = %d, got %d",
				i, exp.ID[i], tags.ID[i])
		}
		if !Int16Eq(tags.Snap[i], exp.Snap[i]) {
			t.Errorf("Expected tags.Snap[%d] = %d, got %d",
				i, exp.Snap[i], tags.Snap[i])
		}
		if !Uint8Eq(tags.Flag[i], exp.Flag[i]) {
			t.Errorf("Expected tags.Flag[%d] = %d, got %d",
				i, exp.Flag[i], tags.Flag[i])
		}
	}
}
