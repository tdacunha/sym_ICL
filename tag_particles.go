package main

import (
	"encoding/binary"
	"fmt"
	"os"
	"runtime"

	"github.com/phil-mansfield/guppy/lib/catio"
	"github.com/phil-mansfield/read_gadget"

	//"research/lmc_ges_tracking/lib"
)

var (
	MaxSnap = 235
	Blocks = 8
	NHaloes = 3
	HRLevel = 1
	
	SnapshotFormat = "path/to/snapshot%d.%d"
	TrackFileName = "track/file/name"

	OutputFormat = "path/to/output/format.%d"
)

// HaloTrack contains the evolution history of a single halo.
type HaloTrack struct {
	X, Y, Z, Rvir []float32
}

func main() {
	ids := InitIDs(NHaloes)

	tracks := ReadTracks(TrackFileName, NHaloes)

	for snap := 0; snap <= MaxSnap; snap++ {
		runtime.GC()
		
		for b := 0; b < Blocks; b++ {
			fileName := fmt.Sprintf(SnapshotFormat, snap, b)
			xp, _, idp := ReadLevel(fileName, HRLevel)

			AddIDs(ids, tracks, snap, xp, idp)
		}
	}
	
	WriteIDs(ids)
}

func InitIDs(n int) []map[int32]bool {
	out := make([]map[int32]bool, n)
	for i := range out {
		out[i] = map[int32]bool{ }
	}
	return out
}

func ReadTracks(fileName string, n int) []HaloTrack {
	f := catio.TextFile(fileName)
	
	idxs := make([]int, n*4)
	for i := range idxs { idxs[i] = i+1 }
	cols := f.ReadFloat32s(idxs)

	out := make([]HaloTrack, n)

	for i := range out {
		out[i].X = cols[4*i]
		out[i].Y = cols[4*i+1]
		out[i].Z = cols[4*i+2]
		out[i].Rvir = cols[4*i+3]
	}

	return out
}

func ReadLevel(fileName string, level int) (xp, vp [][3]float32, idp []int32) {
	f := read_gadget.OpenGadget2Zoom(fileName, []string{"x", "v", "id32"})
	xp = make([][3]float32, f.N[level])
	vp = make([][3]float32, f.N[level])
	idp = make([]int32, f.N[level])

	f.Read("x", level, xp)
	f.Read("v", level, vp)
	f.Read("id32", level, idp)

	return xp, vp, idp
}

func AddIDs(
	ids []map[int32]bool, tracks []HaloTrack, snap int,
	xp [][3]float32, idp []int32,
) {
	r2, x := make([]float32, len(tracks)), make([][3]float32, len(tracks))
	for i := range r2 {
		r2[i] = tracks[i].Rvir[snap] * tracks[i].Rvir[snap]
		x[i] = [3]float32{
			tracks[i].X[snap], tracks[i].Y[snap], tracks[i].Z[snap],
		}
	}
	
	for j := range xp {
		for i := range tracks {
			if _, ok := ids[i][idp[j]]; ok { break }
			
			dr2 := DR2(xp[j], x[i])
			if dr2 <= r2[i] {
				AddID(ids, i, idp[j])
				RemoveIDs(ids, i, idp[j])
				break
			}
		}
	}
}

func AddID(ids []map[int32]bool, i int, id int32) {
	ids[i][id] = true
}

func RemoveIDs(ids []map[int32]bool, i int, id int32) {
	for j := i+1; j < len(ids); j++ {
		delete(ids[j], id)
	}
}

func DR2(x1, x2 [3]float32) float32 {
	dr2 := float32(0)
	for i := range x1 {
		dx := x1[i] - x2[i]
		dr2 += dx*dx
	}
	return dr2
}

func WriteIDs(ids []map[int32]bool) {
	for i := range ids {
		fileName := fmt.Sprintf(OutputFormat, i)

		f, err := os.Create(fileName)
		if err != nil { panic(err.Error()) }

		idSlice := []int32{ }
		for id, _ := range ids[i] { idSlice = append(idSlice, id) }
		
		err = binary.Write(f, binary.LittleEndian, int32(len(idSlice)))
		if err != nil { panic(err.Error()) }
		err = binary.Write(f, binary.LittleEndian, idSlice)
		if err != nil { panic(err.Error()) }

		f.Close()
	}
}
