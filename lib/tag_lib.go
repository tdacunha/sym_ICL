package lib

// PeakIndex returns the index of the largest value in x.
func PeakIndex(x []float32) int {
	iMax := 0
	for i := range x {
		if x[i] > 0 && x[i] > x[iMax] {
			iMax = i
		}
	}
	return iMax
}

// Mpeak returns M_peak given an M_vir history.
func Mpeak(mvir []float32) float32 {
	return mvir[PeakIndex(mvir)]
}

// IsNewOwner returns true if the halo at index i is the new owner of a
// particle, given the indices of its previous owners. mpeak gives the peak
// mass of all the tracked haloes.
func IsNewOwner(mpeak []float32, prevOwners []int32, i int32) bool {
	last := len(prevOwners) - 1
	return len(prevOwners) <= 0 || mpeak[prevOwners[last]] < mpeak[i]
}

// BetterOwner returns which halo index, i or j, is better owner of a given
// particle for a fixed snapshot.
func BetterOwner(mpeak []float32, i, j int32) int32 {
	if mpeak[i] > mpeak[j] {
		return  i
	} else {
		return j
	}
}

type TagWorker struct {
	start, end, skip int32
	x [][3]float32
	id []int32
	boxSize float32
	finder *Finder
	
	// Owners is the index of the halo that each particle belongs to.
	Owners []int32 
}

// NewTagWorker initalizes a new buffer for a box with a given width and a given
// start and skip. Use the same SplitArray() call that you'd use to do whatever
// work you're doing to get this right.
func NewTagWorker(L float32, start, end, skip int) *TagWorker {
	return &TagWorker{
		start: int32(start),
		end: int32(end),
		skip: int32(skip),
		x: [][3]float32{ },
		boxSize: L,
		finder: nil,
	}
}

func (w *TagWorker) ResetBounds(start, end, skip int) {
	w.start = int32(start)
	w.end = int32(end)
	w.skip = int32(skip)
}

// OriginalIndex maps an index in the worker array to the corresponding index
// in the particle array.
func (w *TagWorker) OriginalIndex(i int32) int32 {
	return w.start + i*w.skip
}

// LoadParticles loads the particles x into w's internal buffers and resets
// the Owners buffer.
func (w *TagWorker) LoadParticles(x [][3]float32, id []int32) {
	w.x = w.x[:0]
	w.id = w.id[:0]
	for i := w.start; i < w.end; i += w.skip {
		w.x = append(w.x, x[i])
		w.id = append(w.id, id[i])
	}
	if len(w.x) == 0 {
		w.Owners = w.Owners[:0]
		return
	}
	
	// TODO: make an option that reuses finder space.
	if w.finder == nil {
		L := 3*x[0][0] // L doesn't matter, just needs to be big
		w.finder = NewFinder(L, w.x)
	} else {
		w.finder.Reuse(w.x)
	}

	// Initialize 
	if cap(w.Owners) >= len(w.x) {
		w.Owners = w.Owners[:len(w.x)]
	} else {
		w.Owners = w.Owners[:cap(w.Owners)]

		nNew := len(w.x) - len(w.Owners)
		w.Owners = append(w.Owners, make([]int32, nNew)...)
	}

	for i := range w.Owners { w.Owners[i] = -1 }
}

// FindParticleOwners sets a worker's Owners field given a set of halo
// positions, virial radii, and mpeak values.
func (w *TagWorker) FindParticleOwners(x [][3]float32, r, mpeak []float32) {
	// The Finder object doesn't get reset for empty files.
	if len(w.x) == 0 { return }

	for i := int32(0); i < int32(len(x)); i++ {
		if r[i] <= 0 { continue } // Halo doesn't exist at this snapshot
		idx := w.finder.Find(x[i], r[i])
		for _, j := range idx {
			if w.Owners[j] == -1 {
				w.Owners[j] = i
			} else {
				w.Owners[j] = BetterOwner(mpeak, i, w.Owners[j])
			}
		}
	}
}


// InsertOwnersInLists mergest owner indices/snapshots with the existing
// per-particle lists of these values. snap is the current snapshot and mpeak
// is a per-halo array of M_peak values.
func InsertOwnersInLists(workers []*TagWorker, snap int,
	idxList, snapList *CompactList, mpeak []float32) {

	snap32 := int32(snap)
	
	for _, w := range workers {
		for j := int32(0); j < int32(len(w.Owners)); j++ {
			if w.Owners[j] == -1 { continue }

			idIdx := w.id[j] - 1 // Gadget Ids are 1-indexed

			currOwner := w.Owners[j]
			prevOwner, ok := idxList.Head(idIdx)
			if !ok || (prevOwner != currOwner &&
				currOwner == BetterOwner(mpeak, currOwner, prevOwner)) {

				idxList.Push(idIdx, currOwner)
				snapList.Push(idIdx, snap32)
			}
		}
	}
}




func NewTags(nHalo int) *Tags {
	return &Tags {
		make([]int32, nHalo),
		make([][]int32, nHalo),
		make([][]int16, nHalo),
		make([][]uint8, nHalo),
	}
}

// AddChangedParticles adds particles from the two lists which have changed
// ownership in the current snapshot. The two lists have lengths equal to the 
// 
func (buf *Tags) AddChangedParticles(idxList, snapList *CompactList, snap int) {

	snap32, snap16 := int32(snap), int16(snap)

	for i := range idxList.start {
		idx := idxList.start[i]
		if idx == listEnd { continue }

		snapi := snapList.data[idx]
		haloi := idxList.data[idx]
		flagi := uint8(0)
		
		if idxList.next[idx] != listEnd { flagi = uint8(1) }

		if snapi == snap32 {
			// Gadget IDs are 1-indexed
			buf.ID[haloi] = append(buf.ID[haloi], int32(i+1))
			buf.Snap[haloi] = append(buf.Snap[haloi], snap16)
			buf.Flag[haloi] = append(buf.Flag[haloi], flagi)
		}
	}
}


// OrderTags orders the elements of Idx, Snap, and Flag so that all the 0-tag
// particles are first and all the 1-tag particles are the 
func (tags *Tags) OrderTags() {
	tags.N0 = make([]int32, len(tags.Flag))
	for i := range tags.Flag {
		id, snap, flag := tags.ID[i], tags.Snap[i], tags.Flag[i]

		for j := range flag {
			if flag[j] == 0 { tags.N0[i]++ }
		}

		k0, k1 := 0, int(tags.N0[i])
		n0, n1 := k1, len(tags.Flag[i])

		for { // Loops over all swaps
			for ; k0 < n0 && flag[k0] == 0; k0++ { }
			if k0 == n0 { break }
			for ; k1 < n1 && flag[k1] == 1; k1++ { }
			if k1 == n1 { break }
		
			flag[k0], flag[k1] = flag[k1], flag[k0]
			snap[k0], snap[k1] = snap[k1], snap[k0]
			id[k0], id[k1] = id[k1], id[k0]
		}
	}
}

// TrimFloat filters an array of floating point values to only contain values
// with infall snapshots at or before the given snapshot
func (tags *Tags) TrimFloat(iHalo int, x []float32, snap int) []float32 {
	j := 0
	snap16 := int16(snap)
	for i := range x {
		if tags.Snap[iHalo][i] <= snap16 {
			x[j] = x[i]
			j++
		}
	}
	return x[:j]
	
}

// TrimFloat filters an array of vectors to only contain values
// with infall snapshots at or before the given snapshot.
func (tags *Tags) TrimVector(iHalo int, x [][3]float32, snap int) [][3]float32 {
	j := 0
	snap16 := int16(snap)
	for i := range x {
		if tags.Snap[iHalo][i] <= snap16 {
			x[j] = x[i]
			j++
		}
	}
	return x[:j]
}

func (tags *Tags) ExpandFloat(
	iHalo int, x []float32, snap int, fill float32) []float32 {
	j := len(x) - 1
	snap16 := int16(snap)
	x = append(x, make([]float32, len(tags.Snap[iHalo]) - len(x))...)
	
	for i := len(x) - 1; i >= 0; i-- {
		if tags.Snap[iHalo][i] <= snap16 {
			x[i] = x[j]
			j--
		} else {
			x[i] = fill
		}
	}

	return x
}

func (tags *Tags) ExpandVector(
	iHalo int, x [][3]float32, snap int, fill [3]float32) [][3]float32 {

	j := len(x) - 1
	snap16 := int16(snap)
	x = append(x, make([][3]float32, len(tags.Snap[iHalo]) - len(x))...)
	
	for i := len(x) - 1; i >= 0; i-- {
		if tags.Snap[iHalo][i] <= snap16 {
			x[i] = x[j]
			j--
		} else {
			x[i] = fill
		}
	}

	return x
}

type TagLookup struct {
	Halo []int16
	Index []int32
}

func NewTagLookup(np int) *TagLookup {
	return &TagLookup{
		Halo: make([]int16, np),
		Index: make([]int32, np),
	}
}


func (l *TagLookup) AddTags(tags *Tags) {
	for iHalo := range tags.ID {
		for i := range tags.ID[iHalo] {
			idIdx := tags.ID[iHalo][i] - 1
			l.Halo[idIdx] = int16(iHalo)
			l.Index[idIdx] = int32(i)
		}
	}
}
