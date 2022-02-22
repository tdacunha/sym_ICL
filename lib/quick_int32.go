package lib

// sort3 sorts three values from largest to smallest.
func sort3Int32(x, y, z int32) (max, mid, min int32) {
	if x > y {
		if x > z {
			if y > z {
				return x, y, z
			} else {
				return x, z, y
			}
		} else {
			return z, x, y
		}
	} else {
		if y > z {
			if x > z {
				return y, x, z
			} else {
				return y, z, x
			}
		} else {
			return z, y, x
		}
	}
}

// QuickSort sorts an array in place via quicksort (and returns the result for
// convenience.)
//
// QuickSort is significantly faster than the standard library's quicksort on
// arrays of floats.
func QuickSortInt32(xs []int32) []int32 {
	if len(xs) < manualLen {
		return ShellSortInt32(xs)
	} else {
		pivIdx := partitionInt32(xs)
		QuickSortInt32(xs[0:pivIdx])
		QuickSortInt32(xs[pivIdx:len(xs)])
		return xs
	}
}

// partition rearranges the elements of a slice, xs, into two contiguous
// groups, such that every element of the first group is smaller than every
// element of the second. partition then returns the length of the first group.
func partitionInt32(xs []int32) int {
	n, n2 := len(xs), len(xs)/2
	// Take three values. The median will be the pivot, the other two will
	// be sentinel values so that we can avoid bounds checks.
	max, mid, min := sort3Int32(xs[0], xs[n2], xs[n-1])
	xs[0], xs[n2], xs[n-1] = min, mid, max
	xs[1], xs[n2] = xs[n2], xs[1]

	lo, hi := 1, n-1
	for {
		lo++
		for xs[lo] < mid {
			lo++
		}
		hi--
		for xs[hi] > mid {
			hi--
		}
		if hi < lo {
			break
		}
		xs[lo], xs[hi] = xs[hi], xs[lo]
	}

	// Swap the pivot into the middle
	xs[1], xs[hi] = xs[hi], xs[1]

	return hi
}
