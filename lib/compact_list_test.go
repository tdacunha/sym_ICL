package lib

import (
	"testing"
)

func TestSimpleListFuncs(t *testing.T) {
	start := []int32{listEnd, 0, listEnd, 2, 5 }
	next := []int32{ 1, 3, 4, -1, -1, -1 }
	data := []int32{ 10, 11, 12, 13, 14, 15, 16 }
	l := &CompactList{ start: start, next: next, data: data }

	arrays := [][]int32{
		{},
		{10, 11, 13}, 
		{},
		{12, 14},
		{15},
	}

	heads := []int32{ listEnd, 10, listEnd, 12, 15 }

	arr := []int32{ }
	for i := range start {
		if head, _ := l.Head(int32(i)); head != heads[i] {
			t.Errorf("Expected Head(%d) = %d, got %d", i, heads[i], head)
		}

		arr := l.GetArray(int32(i), arr)
		if !Int32Eq(arr, arrays[i]) {
			t.Errorf("Expected GetArray(%d) = %d, got %d", i, arrays[i], arr)
		}
	}
}

func TestPush(t *testing.T) {
	pushes := []struct {
		id, data int32
	} {
		{0, 10},
			{0, 11},
			{0, 12},
			{2, 13},
			{3, 14},
			{2, 15},
	}
	
	arrays := [][]int32{
		{12, 11, 10},
		{},
		{15, 13},
		{14},
		{},
	}

	l := NewCompactList(5)
	
	for i := range pushes {
		l.Push(pushes[i].id, pushes[i].data)
	}
	
	buf := []int32{ }
	for i := range arrays {
		buf = l.GetArray(int32(i), buf)
		if !Int32Eq(buf, arrays[i]) {
			t.Errorf("Expected GetArray(%d) = %d, got %d.", i, arrays[i], buf)
		}
	}
}
