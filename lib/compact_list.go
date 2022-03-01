package lib

const (
	listEnd = -1
)

// Compact list is a compact data structure for storing many small lists. Avoids
// both excessive malloc overhead and heap fragmentation.
type CompactList struct {
	start, next, data []int32
}

// NewCompactList creates a new CompactList associated with n objects that have
// IDs ranging across [0, n). tailSpeed should be set to FastTailCalculation
// if calcualing Tail() quickly is more important than memory savings.
func NewCompactList(n int32) *CompactList {
	l := &CompactList{ make([]int32, n), nil, nil, }
	for i := range l.start { l.start[i] = listEnd }
	return l
}


// Push adds a piece of data to the object with the given ID.
func (l *CompactList) Push(id, data int32) {
	if l.start[id] == listEnd {
		n := int32(len(l.next))
		l.start[id] = n
		l.next = append(l.next, listEnd)
		l.data = append(l.data, data)
	} else {	
		i := l.start[id]
		l.start[id] = int32(len(l.next))
		l.next = append(l.next, i)
		l.data = append(l.data, data)
	}
}

// Head returns the first data item in a list.
func (l *CompactList) Head(id int32) (int32, bool) {
	if l.start[id] == listEnd { return listEnd, false }
	return l.data[l.start[id]], true
}

// GetArray returns all the elements associated with the given ID in order
// as an array.
func (l *CompactList) GetArray(id int32, buf []int32) []int32 {
	buf = buf[:0]
	if l.start[id] == listEnd { return buf }
	
	for i := l.start[id]; ; i = l.next[i] {
		buf = append(buf, l.data[i])
		if l.next[i] == listEnd { break }
	}
	return buf
}
