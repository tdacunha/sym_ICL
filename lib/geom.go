package lib

func Bound(dx, L float32) float32 {
	if dx > L { return dx - L }
	if dx < 0 { return dx + L }
	return dx
}

func SymBound(dx, L float32) float32 {
	if dx > +L/2 { return dx - L }
	if dx < -L/2 { return dx + L }
	return dx
}

type FloatBounds struct {
	Origin, Span [3]float32
}

// PointBoundsNonPeriodic calculates a bounding box around a set of points
// which does not span box boundaries.
func PointBoundsNonPeriodic(x [][3]float32) *FloatBounds {
	min, max := x[0], x[0]
	for i := 1; i < len(x); i++ {
		for k := 0; k < 3; k++ {
			xx := x[i][k]
			if min[k] > xx {
				min[k] = xx
			} else if max[k] < xx {
				max[k] = xx
			}
		}
	}
	
	b := &FloatBounds{ }
	b.Origin = min
	for k := 0; k < 3; k++ {
		b.Span[k] = max[k] - min[k]
	}
	
	return b
}

// Bounds is a cell-aligned bounding box.
type Bounds struct {
	Origin, Span [3]int
}

// FromFloatBounds converts a FloatBounds object to a Bounds object within
// a grid with a given cell width, cw.
func FloatBoundsToIntBounds(fb *FloatBounds, cw float32) *Bounds {
	min, max := fb.Origin, [3]float32{ }
	for k := 0; k < 3; k++ { max[k] = min[k] + fb.Span[k] }

	b := &Bounds{ }
	for k := 0; k < 3; k++ {
		intMin, intMax := int(min[k]/cw), int(max[k]/cw)
		b.Origin[k] = intMin
		b.Span[k] = intMax - intMin + 1
	}

	return b
}

// SphereBounds creates a cell-aligned bounding box around a non-aligned
// sphere within a box with periodic boundary conditions.
func (b *Bounds) SphereBounds(pos [3]float32, r, cw, width float32) {
	for i := 0; i < 3; i++ {
		min, max := pos[i]-r, pos[i]+r
		if min < 0 {
			min += width
			max += width
		}
		
		minCell, maxCell := int(min/cw), int(max/cw)
		b.Origin[i] = minCell
		b.Span[i] = maxCell - minCell + 1
	}
}

// ConvertIndices converts non-periodic indices to periodic indices.
func (b *Bounds) ConvertIndices(x, y, z, width int) (bx, by, bz int) {
	bx = x - b.Origin[0]
	if bx < 0 {
		bx += width
	}
	by = y - b.Origin[1]
	if by < 0 {
		by += width
	}
	bz = z - b.Origin[2]
	if bz < 0 {
		bz += width
	}
	return bx, by, bz
}

// Inside returns true if the given value is within the bounding box along the
// given dimension. The periodic box width is given by width.
func (b *Bounds) Inside(val int, width int, dim int) bool {
	lo, hi := b.Origin[dim], b.Origin[dim]+b.Span[dim]
	if val >= hi {
		val -= width
	} else if val < lo {
		val += width
	}
	return val < hi && val >= lo
}
