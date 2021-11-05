package main

import (
	"fmt"
	"log"
	"os"
	"runtime"

	"research/lmc_ges_tracking/lib"
	"github.com/phil-mansfield/guppy/lib/catio"
)

var (
	MaxSnap = 235
	TreeFileNames = []string{ "/scratch/users/enadler/Halo416/rockstar/trees/tree_0_0_0.dat" }
	CentralID = 9285216
	SubhaloIDs = []int{ 9278166, 9285216 }
	
	HaloIDs = append([]int{CentralID}, SubhaloIDs...)

	OutputFileName = "/scratch/users/phil1/lmc_ges_tracking/Halo416/halo_tracks.txt"
)

// Haloes is a collection of raw columns from the tree.dat files.
type Haloes struct {
	ID, DFID, Snap []int
	Rvir, X, Y, Z []float32
}

// HaloTrack contains the evolution history of a single halo.
type HaloTrack struct {
	RootID int
	X, Y, Z, Rvir []float32
}

func main() {
	tracks := []HaloTrack{ }

	cfg := catio.DefaultConfig
	cfg.SkipLines = 46

	// Loop over different tree files
	for i := range TreeFileNames {
		log.Println("Reading", TreeFileNames[i])
		rd := catio.TextFile(TreeFileNames[i], cfg)

		// If the tree files are too large to load into memory at once,
		// break it up and read each block separately.
		for b := 0; b < rd.Blocks(); b++ {
			runtime.GC()

			// Read the columns
			log.Printf("Read ints for block %d", b)
			iCols := rd.ReadIntBlock([]int{1, 28, 31}, b)
			id, dfid, snap := iCols[0], iCols[1], iCols[2]
			
			log.Printf("Read ints for block %d", b)
			fCols := rd.ReadFloat32Block([]int{11, 17, 18, 19}, b)
			rvir, x, y, z := fCols[0], fCols[1], fCols[2], fCols[3]

			h := &Haloes{ id, dfid, snap, x, y, z, rvir }

			log.Println("Getting tracks")
			// Find tracks from the currently-read block
			tracks = append(tracks, GetTracks(h, HaloIDs)...)
		}
	}

	log.Println("Writing")
	WriteTracks(tracks, CentralID, OutputFileName)
}

// GetTracks returns all the tracks for haloes whose roots are within rootIDs.
func GetTracks(h *Haloes, haloIDs []int) []HaloTrack {
	tracks := []HaloTrack{ }

	order := IDOrder(h.DFID)
	log.Println("Finished sorting")

	for j := range order {
		i := order[j]

		for k := range haloIDs {
			if h.DFID[i] == haloIDs[k] {
				log.Println("Finished halo track", k)
				tracks = append(tracks, GetSingleTrack(h, order, j))
			}
		}
	}

	return tracks
}

func IDOrder(id []int) []int {
	fid := make([]float64, len(id))
	for i := range id {
		// This works as long as the maximum ID is smaller than the manstissa
		// of a double-precision float.
		if id[i] > (1<<52) {
			panic("You need to write array.QuickSort so it works on ints: " +
				"there's an ID larger than the double-precision manstissa.")
		}
		fid[i] = float64(id[i])
	}

	return lib.QuickSortIndex(fid)
}

func GetSingleTrack(h *Haloes, order []int, j0 int) HaloTrack {
	n := MaxSnap + 1
	rootID := h.ID[order[j0]]
	track := HaloTrack{
		rootID, make([]float32, n), make([]float32, n),
		make([]float32, n), make([]float32, n),
	}
	
	for i := range track.X {
		track.X[i], track.Y[i], track.Z[i], track.Rvir[i] = -1, -1, -1, -1
	}
	
	jStart := FirstBranchIndex(h.DFID, h.Snap, order, j0)
	jEnd := LastBranchIndex(h.DFID, h.Snap, order, j0)

	for j := jStart; j <= jEnd; j++ {
		i := order[j]
		snap := h.Snap[i]

		track.Z[snap] = h.Z[i]
		track.Y[snap] = h.Y[i]
		track.Z[snap] = h.Z[i]
		track.Rvir[snap] = h.Rvir[i]
	}

	return track
}

// j0 indexes into order. Returns an index into order
func LastBranchIndex(dfid, snap, order []int, j0 int) int {
	for j := j0; ; j++ {
		if j == len(order) - 1 { return j }

		currDFID, currSnap := dfid[order[j]], snap[order[j]]
		nextDFID, nextSnap := dfid[order[j + 1]], snap[order[j + 1]]

		if nextDFID - 1 != currDFID || nextSnap >= currSnap { return j }
	}
	panic("Impossible")
}

// j0 indexes into order. Returns an index into order
func FirstBranchIndex(dfid, snap, order []int, j0 int) int {
	for j := j0;  ; j-- {
		if j == 0 { return j }

		currDFID, currSnap := dfid[order[j]], snap[order[j]]
		nextDFID, nextSnap := dfid[order[j - 1]], snap[order[j - 1]]

		if nextDFID + 1 != currDFID || nextSnap <= currSnap { return j }
	}
	panic("Impossible")
}

// WriteTracks writes a given set of tracks to the given output file. The
// central is written first.
func WriteTracks(tracks []HaloTrack, centralRootID int, fileName string) {
	f, err := os.Create(fileName)
	if err != nil { panic(err.Error()) }
	defer f.Close()
	
	// Find the index of the central. I think it will always be zero.
	var iCentral int
	for iCentral = 0; iCentral < len(tracks); iCentral++ {
		if tracks[iCentral].RootID == centralRootID { break }
	}


	fmt.Fprintf(f, 
`# This file contains the evolutionary tracks of a central halo  and several of
# its subhaloes. The first column gives the snapshot number, and every
# subsequent four columns give the x, y, and z coordinates and Rvir of a
# different halo during that snapshot. The first set of four columns belong to
# the central and each subsequent group of four is another subhalo. All values
# are given in comoving Mpc/h and are set to -1 in snapshots where the subhalo
# does not exist.\n`)
	
	for snap := 0; snap <= MaxSnap; snap++ {
		c := tracks[iCentral]
		fmt.Fprintf(f, "%3d ", snap)
		fmt.Fprintf(f, "%8.5f %8.5f %8.5f %7.5f ",
			c.X[snap], c.Y[snap], c.Z[snap], c.Rvir[snap])
		
		for i := range tracks {
			if i == iCentral { continue }

			s := tracks[i]
			fmt.Fprintf(f, "%8.5f %8.5f %8.5f %7.5f ",
				s.X[snap], s.Y[snap], s.Z[snap], s.Rvir[snap])
		}
		
		fmt.Fprintf(f, "\n")
	}
}
