package lib

import (
	"testing"
)

func TestFinder(t *testing.T) {
	L := float32(200)
	x := [][3]float32{
		{ 50,  50,  50},
		{150, 150, 150},
		{100, 100, 100},
		{100, 100, 101},
		{100, 100, 102},
		{100, 100, 103},
		{100, 100, 105},
		{100, 100, 106},
	}

	f := NewFinder(L, x)

	idx0 := []int32{ 2, 3, 4, 5 }
	idx1 := []int32{ 0, 2, 3, 4, 5, 6, 7, 1 }
	
	idx := f.Find([3]float32{100, 100, 100}, 4)
	if !Int32Eq(idx, idx0) {
		t.Errorf("expected Find(1) to give %d, but got %d", idx0, idx)
	}
	
	idx = f.Find([3]float32{100, 100, 100}, 99)
	if !Int32Eq(idx, idx1) {
		t.Errorf("expected Find(1) to give %d, but got %d", idx1, idx)
	}

	for i := range x {
		for k := 0; k < 3; k++ {
			x[i][k] = (x[i][k] - 50)*2
		}
	}
	x[0] = [3]float32{1, 1, 1}
	x[1] = [3]float32{199, 199, 199}
	x = append(x, [3]float32{198, 198, 198})

	f.Reuse(x)
	idx0 = []int32{ 2, 3, 4, 5 }
	idx1 = []int32{ 2, 3, 4, 5, 6, 7 }

	idx = f.Find([3]float32{100, 100, 100}, 8)
	if !Int32Eq(idx, idx0) {
		t.Errorf("expected Find(1) to give %d, but got %d", idx0, idx)
	}
	
	idx = f.Find([3]float32{100, 100, 100}, 99)
	if !Int32Eq(idx, idx1) {
		t.Errorf("expected Find(1) to give %d, but got %d", idx1, idx)
	}
}
