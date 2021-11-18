package main

import (
	"encoding/binary"
	"fmt"
	"log"
	"os"
	"runtime"
	"time"

	"github.com/phil-mansfield/lmc_ges_tracking/lib"
)

var (
	MaxSnap = 235
	Blocks = 8
	NHaloes = 3
	HRLevel = 1
	
	SnapshotFormat = "/scratch/users/enadler/Halo416/output/snapshot_%03d.%d"
	TrackFileName = "/scratch/users/phil1/lmc_ges_tracking/Halo416/halo_tracks.txt"
	OutputFormat = "/scratch/users/phil1/lmc_ges_tracking/Halo416/ids.%d"
)

// HaloTrack contains the evolution history of a single halo.
type HaloTrack struct {
	X, Y, Z, Rvir []float32
}

func main() {
	ids := InitIDs(NHaloes)

	tracks := ReadTracks(TrackFileName, NHaloes)

	for snap := 0; snap <= MaxSnap; snap++ {
		fmt.Println("tag", snap)
		runtime.GC()
		
		dtRead := 0.0
		dtAdd := 0.0

		for b := 0; b < Blocks; b++ {
			fileName := fmt.Sprintf(SnapshotFormat, snap, b)
			t1 := time.Now()
			xp, idp := ReadLevel(fileName, HRLevel)
			t2 := time.Now()
			AddIDs(ids, tracks, snap, xp, idp)
			t3 := time.Now()

			dtAdd += t3.Sub(t2).Seconds()
			dtRead += t2.Sub(t1).Seconds()
		}
		log.Printf("%3d %7d %7d %7d\n", 
			snap, len(ids[0]), len(ids[1]), len(ids[2]))
		log.Printf("%5.3f %5.3f\n", dtRead, dtAdd)
	}

	WriteIDs(ids)
}

func InitIDs(n int) []map[int32]int {
	out := make([]map[int32]int, n)
	for i := range out {
		out[i] = map[int32]int{ }
	}
	return out
}

func ReadTracks(fileName string, n int) []HaloTrack {
	fmt.Println(fileName, n)
	f := lib.TextFile(fileName)
	
	idxs := make([]int, n*7)
	for i := range idxs { idxs[i] = i }
	cols := f.ReadFloat32s(idxs)

	out := make([]HaloTrack, n)

	for i := range out {
		out[i].X = cols[7*i]
		out[i].Y = cols[7*i+1]
		out[i].Z = cols[7*i+2]
		out[i].Rvir = cols[7*i+6]
	}

	return out
}

func ReadLevel(fileName string, level int) (xp [][3]float32, idp []int32) {
	f := lib.OpenGadget2Zoom(fileName, []string{"x", "v", "id32"})
	xp = make([][3]float32, f.N[level])
	idp = make([]int32, f.N[level])

	f.Read("x", level, xp)
	f.Read("id32", level, idp)

	return xp, idp
}

func AddIDs(
	ids []map[int32]int, tracks []HaloTrack, snap int,
	xp [][3]float32, idp []int32,
) {
	r2, x := make([]float32, len(tracks)), make([][3]float32, len(tracks))
	for i := range r2 {
		r2[i] = tracks[i].Rvir[snap] * tracks[i].Rvir[snap]
		x[i] = [3]float32{
			tracks[i].X[snap], tracks[i].Y[snap], tracks[i].Z[snap],
		}
	}

ParticleLoop:	
	for j := range xp {
		for i := range tracks {
			if _, ok := ids[i][idp[j]]; ok {
				continue ParticleLoop
			}
		}
		for i := range tracks {
			if tracks[i].Rvir[snap] < 0 { continue }

			id := idp[j]
			dr2 := DR2(xp[j], x[i])

			if dr2 <= r2[i] {
				ids[i][id] = snap
				break
			}
		}
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

func WriteIDs(ids []map[int32]int) {
	for i := 0; i < len(ids); i++ {
		fileName := fmt.Sprintf(OutputFormat, i)
		f, err := os.Create(fileName)
		if err != nil { panic(err.Error()) }

		idSlice := []int32{ }
		snapSlice := []int16{ }
		for id, snap := range ids[i] {
			idSlice = append(idSlice, id)
			snapSlice = append(snapSlice, int16(snap))
		}

		err = binary.Write(f, binary.LittleEndian, int32(len(idSlice)))
		if err != nil { panic(err.Error()) }
		err = binary.Write(f, binary.LittleEndian, idSlice)
		if err != nil { panic(err.Error()) }
		err = binary.Write(f, binary.LittleEndian, snapSlice)
		if err != nil { panic(err.Error()) }

		f.Close()
	}
}
