package lib

// PeakIndex returns the index of the largest value in x.
func PeakIndex(x []float32) int {
	iMax := 0
	for i := range x {
		if x[i] > 0 && x[i] > x[iMax] {
			x[iMax] = x[i]
		}
	}
	return iMax
}

// MPeak returns M_peak given an M_vir history.
func MPeak(mvir []float32) float32 {
	return mvir[PeakIndex(mvir)]
}

// IsNewOwner returns true if the halo at index i is the new owner of a
// particle, given the indices of its previous owners. mpeak gives the peak
// mass of all the tracked haloes.
func IsNewOwner(mpeak []float32, prevOwners []int32, i int32) bool {
	last := len(prevOwners) - 1
	return len(prevOwners) < 0 || mpeak[prevOwners[last]] < mpeak[i]
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

// OriginalIndex maps an index in the worker array to the corresponding index
// in the particle array.
func (w *TagWorker) OriginalIndex(i int32) int32 {
	return w.start + i*w.skip
}

// LoadParticles loads the particles x into w's internal buffers and resets
// the Owners buffer.
func (w *TagWorker) LoadParticles(x [][3]float32) {
	w.x = w.x[:0]
	
	for i := w.start; i < w.end; i += w.skip {
		w.x = append(w.x, x[i])
	}

	// TODO: make an option that reuses finder space.
	w.finder.Reuse(w.x)

	// Initialize 
	
	if cap(w.Owners) >= len(w.x) {
		w.Owners = w.Owners[:len(w.x)]
	} else {
		w.Owners = w.Owners[:cap(w.x)]

		nNew := len(w.x) - len(w.Owners)
		w.Owners = append(w.Owners, make([]int32, nNew)...)
	}

	for i := range w.Owners { w.Owners[i] = -1 }
}

// FindParticleOwners sets a worker's Owners field given a set of halo
// positions, virial radii, and mpeak values.
func (w *TagWorker) FindParticleOwners(x [][3]float32, r, mpeak []float32) {
	for i := int32(0); i < int32(len(x)); i++ {
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
			jOrig := w.OriginalIndex(j)

			currOwner := w.Owners[j]
			prevOwner, ok := idxList.Head(jOrig)
			if !ok || (prevOwner != currOwner &&
				currOwner == BetterOwner(mpeak, currOwner, prevOwner)) {
				
				idxList.Push(jOrig, currOwner)
				snapList.Push(jOrig, snap32)
			}
		}
	}
}

// NOTE: looking at this, we could actually refactor the code to not keep track
// of the CompactList.

type HaloParticlesBuffer struct {
	IDs [][]int32
	Snaps [][]int16
	Flags [][]uint8
}

func NewHaloParticlesBuffer(nHalo int) *HaloParticlesBuffer {
	return &HaloParticlesBuffer {
		make([][]int32, nHalo),
		make([][]int16, nHalo),
		make([][]uint8, nHalo),
	}
}

// CollectHaloParticles adds particles from the two lists which have changed
// ownership in the current snapshot. ids give the IDs of the particles.
func (buf *HaloParticlesBuffer) CollectChangedParticles(
	ids []int32, idxList, snapList *CompactList, snap int) {

	snap32, snap16 := int32(snap), int16(snap)
	
	for i := range idxList.start {
		idx := idxList.start[i]
		if idx == listEnd { continue }

		snapi := snapList.data[idx]
		haloi := idxList.data[idx]
		flagi := uint8(0)
		if idxList.next[idx] != listEnd { flagi = uint8(1) }
		
		if snapi == snap32 {
			buf.IDs[haloi] = append(buf.IDs[haloi], ids[i])
			buf.Snaps[haloi] = append(buf.Snaps[haloi], snap16)
			buf.Flags[haloi] = append(buf.Flags[haloi], flagi)
		}
	}
}
