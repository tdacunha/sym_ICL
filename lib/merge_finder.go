package lib

// MergePair merges two sorted arrays into a single sorted array.
func MergePair(x, y []int32) []int32 {
	out := make([]int32, len(x)+len(y))

	ix, iy, j := 0, 0, 0
	for {
		// Handle the cases where we've reached the end of one array or
		// the other
		if ix == len(x) {
			for ; iy < len(y); iy++ {
				out[j] = y[iy]
				j++
			}
			break
		} else if iy == len(y) {
			for ; ix < len(x); ix++ {
				out[j] = x[ix]
				j++
			}
			break
		}

		if x[ix] >= y[iy] {
			out[j] = y[iy]
			iy++
		} else {
			out[j] = x[ix]
			ix++
		}
		j++
	}

	return out
}

// MergeAll merges an arbitray number of sorted arrays together into a single
// sorted array.
func MergeAll(xs [][]int32) []int32 {
	// The most common case has one of the arrays being much longer than the
	// others. This case handles that.

	iMax := 0
	for i := range xs {
		if len(xs[iMax]) > len(xs[i]) {
			iMax = i
		}
	}

	// This is slow if you have a very long xs array. COme back to this later
	// if needed.
	var m []int32
	for i := range xs {
		if i == iMax { continue }
		m = MergePair(m, xs[i])
	}

	return MergePair(m, xs[iMax])
}

// UniqueIDs returns the IDs in x that are not contained within ref. Both ref
// and x must be sorted.
func UniqueIDs(ref, x []int32) []int32 {
	out := []int32{ }
	ix, iref := 0, 0
	for {
		if ix == len(x) { return out }
		if iref == len(ref) {
			return append(out, x[ix:]...)
		}

		if x[ix] < ref[iref] {
			out = append(out, x[ix])
			ix++
		} else if x[ix] == ref[iref] {
			ix++
			iref++
		} else if x[ix] > ref[iref] {
			
		}
	}

	return out
}
