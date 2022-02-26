package lib

type Grid struct {
	Bounds
	Cells     int
	cw, Width float32

	// Grid-sized
	Heads []int32
	// Data-sized
	Next []int32
}

func NewGrid(b *Bounds, cells int, width float32, dataLen int) *Grid {
	g := &Grid{
		Bounds: *b,
		Cells: cells,
		cw:    width / float32(cells),
		Width: width,
		Heads: make([]int32, b.Span[0]*b.Span[1]*b.Span[2]),
		Next:  make([]int32, dataLen),
	}

	for i := range g.Heads {
		g.Heads[i] = listEnd
	}
	
	return g
}

func (g *Grid) Reuse(b *Bounds, cells int, dataLen int) {
	newHeads := g.Heads[:0]
	nHeads := b.Span[0]*b.Span[1]*b.Span[2]
	newNext := g.Next[:0]
	
	if cap(newHeads) > nHeads {
		newHeads = newHeads[:nHeads]
	} else {
		newHeads = newHeads[:cap(newHeads)]
		newHeads = append(newHeads, make([]int32, nHeads - cap(newHeads))...)
	}

	if cap(newNext) > dataLen {
		newNext = newNext[:dataLen]
	} else {
		newNext = newNext[:cap(newNext)]
		newNext = append(newNext, make([]int32, dataLen - cap(newNext))...)
	}
	
	*g = Grid{
		Bounds: *b,
		Cells: cells,
		cw: g.Width / float32(cells),
		Width: g.Width,
		Heads: newHeads,
		Next: newNext,
	}

	for i := range g.Heads {
		g.Heads[i] = listEnd
	}
	
}

func (g *Grid) Length(idx int) int {
	next := g.Heads[idx]
	n := 0
	for next != listEnd {
		n++
		next = g.Next[next]
	}
	return n
}

func (g *Grid) Insert(xs [][3]float32) {
	for i := range xs {
		x, y, z := xs[i][0], xs[i][1], xs[i][2]

		ix, iy, iz := int(x/g.cw), int(y/g.cw), int(z/g.cw)
		idx := (ix - g.Origin[0]) +
			(iy - g.Origin[1])*g.Span[0] +
			(iz - g.Origin[2])*g.Span[0]*g.Span[1]
		
		g.Next[i] = g.Heads[idx]
		g.Heads[idx] = int32(i)
	}
}

func (g *Grid) TotalCells() int {
	return len(g.Heads)
}

func (g *Grid) MaxLength() int {
	max := 0
	for i := 0; i < g.TotalCells(); i++ {
		l := g.Length(i)
		if l > max {
			max = l
		}
	}
	return max
}

func (g *Grid) ReadIndices(idx int, buf []int32) []int32 {
	buf = buf[:cap(buf)]

	next := g.Heads[idx]
	n := 0
	for next != listEnd {
		buf[n] = next
		n++
		next = g.Next[next]
	}

	return buf[:n]
}
