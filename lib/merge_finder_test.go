package lib

import (
	"testing"
	"math/rand"
)

func TestMergePair(t *testing.T) {
	tests := []struct {
		X, Y, Res []int32
	} {
		{[]int32{}, []int32{}, []int32{}},
			{[]int32{1}, []int32{}, []int32{1}},
			{[]int32{}, []int32{1}, []int32{1}},
			{[]int32{}, []int32{1, 2, 3}, []int32{1, 2, 3}},
			{[]int32{1, 2, 3}, []int32{}, []int32{1, 2, 3}},
			{[]int32{1, 2, 3}, []int32{4, 5, 6}, []int32{1, 2, 3, 4, 5, 6}},
			{[]int32{4, 5, 6}, []int32{1, 2, 3}, []int32{1, 2, 3, 4, 5, 6}},
			{[]int32{1, 4, 6}, []int32{2, 3, 5}, []int32{1, 2, 3, 4, 5, 6}},
			{[]int32{2, 2, 3}, []int32{1, 2, 3, 4}, []int32{1, 2, 2, 2, 3, 3, 4}},
			{[]int32{1, 2, 3, 4}, []int32{2, 2, 3}, []int32{1, 2, 2, 2, 3, 3, 4}},
		}
	
	for i := range tests {
		res := MergePair(tests[i].X, tests[i].Y)
		if !Int32Eq(res, tests[i].Res) {
			t.Errorf("%d) Expected PairMerge(%d, %d) = %d, got %d.",
				i+1, tests[i].X, tests[i].Y, tests[i].Res, res)
		}
	}
}

func Int32Eq(x, y []int32) bool {
	if len(x) != len(y) { return false }
	for i := range x {
		if x[i] != y[i] { return false }
	}
	return true
}

func BenchmarkQuickSort3e6(b *testing.B) {
	n := 3000000
	x := make([]int32, n)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		b.StopTimer()

		for i := range x {
			x[i] = int32(rand.Intn(int(1.1*float64(n))))
		}
		
		b.StartTimer()
		
		QuickSortInt32(x)
	}
}

func BenchmarkQuickSort1e6(b *testing.B) {
	n := 1000000
	x := make([]int32, n)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		b.StopTimer()

		for i := range x {
			x[i] = int32(rand.Intn(int(1.1*float64(n))))
		}
		
		b.StartTimer()
		
		QuickSortInt32(x)
	}
}

func BenchmarkQuickSort1e4(b *testing.B) {
	n := 10000
	x := make([]int32, n)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		b.StopTimer()

		for i := range x {
			x[i] = int32(rand.Intn(int(1.1*float64(n))))
		}
		
		b.StartTimer()
		
		QuickSortInt32(x)
	}
}

func BenchmarkQuickSort1e2(b *testing.B) {
	n := 100
	x := make([]int32, n)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		b.StopTimer()

		for i := range x {
			x[i] = int32(rand.Intn(int(1.1*float64(n))))
		}
		
		b.StartTimer()
		
		QuickSortInt32(x)
	}
}

func BenchmarkMap1e6(b *testing.B) {
	n := 1000000
	x := make([]int32, n)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		b.StopTimer()

		m := map[int32]bool { }
		for i := 0; i < n; i++ {
			m[int32(rand.Intn(int(1.1*float64(n))))] = true
		}
		for i := range x {
			x[i] = int32(rand.Intn(int(1.1*float64(n))))
		}
		
		b.StartTimer()
		
		for i := range x {
			i32 := int32(i)
			if _, ok := m[i32]; !ok {
				m[i32] = true
			}
		}
	}
}


func BenchmarkMap3e6(b *testing.B) {
	n := 3000000
	x := make([]int32, n)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		b.StopTimer()

		m := map[int32]bool { }
		for i := 0; i < n; i++ {
			m[int32(rand.Intn(int(1.1*float64(n))))] = true
		}
		for i := range x {
			x[i] = int32(rand.Intn(int(1.1*float64(n))))
		}
		
		b.StartTimer()
		
		for i := range x {
			i32 := int32(i)
			if _, ok := m[i32]; !ok {
				m[i32] = true
			}
		}
	}
}


func BenchmarkMap1e4(b *testing.B) {
	n := 100
	x := make([]int32, n)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		b.StopTimer()

		m := map[int32]bool { }
		for i := 0; i < n; i++ {
			m[int32(rand.Intn(int(1.1*float64(n))))] = true
		}
		for i := range x {
			x[i] = int32(rand.Intn(int(1.1*float64(n))))
		}
		
		b.StartTimer()
		
		for i := range x {
			i32 := int32(i)
			if _, ok := m[i32]; !ok {
				m[i32] = true
			}
		}
	}
}

func BenchmarkMap1e2(b *testing.B) {
	n := 10000
	x := make([]int32, n)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		b.StopTimer()

		m := map[int32]bool { }
		for i := 0; i < n; i++ {
			m[int32(rand.Intn(int(1.1*float64(n))))] = true
		}
		for i := range x {
			x[i] = int32(rand.Intn(int(1.1*float64(n))))
		}
		
		b.StartTimer()
		
		for i := range x {
			i32 := int32(i)
			if _, ok := m[i32]; !ok {
				m[i32] = true
			}
		}
	}
}
