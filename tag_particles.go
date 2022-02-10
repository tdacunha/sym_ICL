package main

import (
	"encoding/binary"
	"fmt"
	"log"
	"os"
	"runtime"
	"time"
	"strconv"
	"io/ioutil"
	"strings"

	"github.com/phil-mansfield/lmc_ges_tracking/lib"
)

var (
	HRLevel = 1

	/*
	MaxSnap = 235
	Blocks = 8
	
	SnapshotFormat = "/scratch/users/enadler/Halo416/output/snapshot_%03d.%d"
	MergerFileName = "/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/Halo416/mergers.dat"
	OutputFormat = "/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/Halo416/particles/ids.%d"
*/
)

// HaloTrack contains the evolution history of a single halo.
type HaloTrack struct {
	X, Y, Z, Rvir []float32
}

func main() {
	haloIndex := -1
    if len(os.Args) > 3 {
        panic("You must supply a file with (1) comovoing force softneign scales, (2) particle masses, (3) block numbers, (4) snapshot formats, (5) merger names, (6) ID-file formats, (7) output formats.")
    } else if len(os.Args) == 3 {
		var err error
		haloIndex, err = strconv.Atoi(os.Args[2])
		if err != nil { panic(err.Error()) }
	}

	inputName := os.Args[1]
	
	_, _, blocks, snapFmt, mergerName, idFmt, _ :=
		ParseInputFile(inputName)

	for i := range blocks {
		if haloIndex != -1 && haloIndex != i { continue }
		log.Println("Analyzing halo", i)
		AnalyzeHalo(blocks[i], snapFmt[i], mergerName[i], idFmt[i])
	}
}


func ParseInputFile(fname string) (eps, mp []float64, blocks []int, snapFmt, mergerName, idFmt, outFmt []string) {

	b, err := ioutil.ReadFile(fname)
	if err != nil {
		panic(fmt.Sprintf("Could not open input file: %s", err.Error()))
	}
	s := string(b)
	
	lines := strings.Split(s, "\n")
	for i := range lines {
		line := strings.Trim(lines[i], " ")
		if len(line) == 0 { continue }

		tok := strings.Split(lines[i], " ")
		cols := []string{ }
		for i := range tok {
			if len(tok[i]) > 0 {
				cols = append(cols, tok[i])
			}
		}
		
		if len(cols) != 7 {
			panic(fmt.Sprintf("Line %d of %s is '%s', but you need there " +
				"to be four columns.", i+1, fname, line))
		}

		epsi, err := strconv.ParseFloat(cols[0], 64)
		if err != nil {
			panic(fmt.Sprintf("Could not parse eps on line %d of " +
				"%s: %s", i+1, fname, cols[0]))
		}
		eps = append(eps, epsi)

		mpi, err := strconv.ParseFloat(cols[1], 64)
		if err != nil {
			panic(fmt.Sprintf("Could not parse mp on line %d of " +
				"%s: %s", i+1, fname, cols[1]))
		}
		mp = append(mp, mpi)

		blocki, err := strconv.Atoi(cols[2])
		if err != nil {
			panic(fmt.Sprintf("Could not parse block number on line %d of " +
				"%s: %s", i+1, fname, cols[2]))
		}
		blocks = append(blocks, blocki)	
		
		snapFmt = append(snapFmt, cols[3])
		mergerName = append(mergerName, cols[4])
		idFmt = append(idFmt, cols[5])
		outFmt = append(outFmt, cols[6])
	}
	
	return eps, mp, blocks, snapFmt, mergerName, idFmt, outFmt
}

func AnalyzeHalo(blocks int, snapFmt, mergerName, outFmt string) {
	log.Println("Blocks:", blocks)
	log.Println("SnapshotFormat:", snapFmt)
	log.Println("MergerFileName:", mergerName)
	log.Println("OutputFormat:", outFmt)
	mergers := lib.ReadMergers(mergerName)

	ids := InitIDs(mergers.Haloes)

	// .____.
	tracks := make([]HaloTrack, mergers.Haloes)
	for i := range tracks {
		x := [3][]float32{  }
		for k := 0; k < 3; k++ {
			x[k] = make([]float32, mergers.Snaps)
			for j := range x[k] {
				x[k][j] =  mergers.X[i][j][k]
			}
		}

		tracks[i] = HaloTrack{ 
			x[0], x[1], x[2], mergers.Rvir[i],
		}
	}

	for snap := 0; snap < len(mergers.Mvir[0]); snap++ {
		log.Printf("    Snap %d", snap)
		runtime.GC()
		
		dtRead := 0.0
		dtAdd := 0.0
		nRead := 0

		for b := 0; b < blocks; b++ {
			fileName := fmt.Sprintf(snapFmt, snap, b)
			t1 := time.Now()
			xp, idp := ReadLevel(fileName, HRLevel)
			nRead += len(xp)
			t2 := time.Now()
			AddIDs(ids, tracks, snap, xp, idp)
			t3 := time.Now()
			dtAdd += t3.Sub(t2).Seconds()
			dtRead += t2.Sub(t1).Seconds()
		}

		log.Printf("    %5.3f %5.3f\n", dtRead, dtAdd)
	}
	sizes := make([]int, len(tracks)) 
	for i := range sizes { sizes[i] = len(ids[i]) }
	log.Printf("halo sizes: %d\n", sizes)

	WriteIDs(outFmt, ids)
}

func BoundingBox(xp [][3]float32) (low, high [3]float32) {
	low, high = xp[0], xp[0]
	for i := range xp {
		for k := 0; k < 3; k++ {
			if xp[i][k] < low[k] {
				low[k] = xp[i][k]
			} else if xp[i][k] > high[k] {
				high[k] = xp[i][k]
			}
		}
	}
	return low, high
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

func WriteIDs(outFmt string, ids []map[int32]int) {
	for i := 0; i < len(ids); i++ {
		fileName := fmt.Sprintf(outFmt, i)
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
