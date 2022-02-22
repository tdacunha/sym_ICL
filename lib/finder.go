package lib

const (
	defaultFinderCells = 100
)

// Finder finds the halos in group A that are within the radii of
// group B. It is optimized for the case where A is extremely large and search
// radii are also very large.
//
// To do this, it does not allow for host list identifications and does not
// memoize results.
type Finder struct {
	g      *Grid
	gBuf   []int32
	idxBuf []int32
	dr2Buf []float32
	x      [][3]float32
	bufi   int
	cells  int
}

func getBounds(x [][3]float32, L float32) (*Bounds, int) {
	fb := PointBoundsNonPeriodic(x)
	maxSpan := fb.Span[0]
	if maxSpan < fb.Span[1] { maxSpan = fb.Span[1] }
	if maxSpan < fb.Span[2] { maxSpan = fb.Span[2] }

	cw := maxSpan / defaultFinderCells
	cells := int(L / cw)
	
	return FloatBoundsToIntBounds(fb, L/float32(cells)), cells
}

// NewFinder creates a new Finder corresponding to the given
// Grid. The Grid contains halos from group A.
func NewFinder(L float32, x [][3]float32) *Finder {
	b, cells  := getBounds(x, L)
	g := NewGrid(b, cells, L, len(x))
	g.Insert(x)
	
	f := &Finder{
		g: g,
		gBuf: make([]int32, g.MaxLength()),
		idxBuf: []int32{ },
		x: x, cells: cells,
	}

	return f
}

// FindSubhalos links grid halos (from group A) to a target halo (from group B).
// Returned array is an internal buffer, so please treat it kindly.
func (sf *Finder) Find(pos [3]float32, r0 float32) []int32 {
	sf.idxBuf = sf.idxBuf[:0]

	b := &Bounds{}
	c := sf.cells
	g := sf.g
	
	b.SphereBounds(pos, r0, g.cw, sf.g.Width)

	for dz := 0; dz < b.Span[2]; dz++ {
		z := b.Origin[2] + dz
		if z >= c { z -= c }
		
		if z < g.Origin[2] || z >= g.Span[2] + g.Origin[2] { continue }
		zOff := (z - g.Origin[2])*g.Span[0]*g.Span[1]
		
		for dy := 0; dy < b.Span[1]; dy++ {
			y := b.Origin[1] + dy
			if y >= c { y -= c }
			
			if y < g.Origin[1] || y >= g.Span[1] + g.Origin[1] { continue }
			yOff := (y - g.Origin[1])*g.Span[0]
			
			for dx := 0; dx < b.Span[0]; dx++ {
				x := b.Origin[0] + dx
				if x >= c { x -= c }

				if x < g.Origin[0] || x >= g.Span[0] + g.Origin[0] { continue }
				xOff := x - g.Origin[0]
				
				idx := zOff + yOff + xOff

				sf.gBuf = sf.g.ReadIndices(idx, sf.gBuf)
				sf.addSubhalos(sf.gBuf, pos[0], pos[1], pos[2], r0, sf.g.Width)
			}
		}
	}

	return sf.idxBuf
}

func (sf *Finder) addSubhalos(
	idxs []int32, xh, yh, zh, rh float32, L float32,
) {
	pL2, nL2 := L/2, -L/2
	for _, j := range idxs {
		sx, sy, sz := sf.x[j][0], sf.x[j][1], sf.x[j][2]
		dx, dy, dz, dr := xh-sx, yh-sy, zh-sz, rh

		if dx > pL2 {
			dx -= L
		} else if dx < nL2 {
			dx += L
		}

		if dy > pL2 {
			dy -= L
		} else if dy < nL2 {
			dy += L
		}

		if dz > pL2 {
			dz -= L
		} else if dz < nL2 {
			dz += L
		}
		
		dr2 := dx*dx + dy*dy + dz*dz

		if dr*dr >= dr2 {
			sf.idxBuf = append(sf.idxBuf, j)
		}
	}
}
