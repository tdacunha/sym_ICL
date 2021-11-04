package main

import (
	"fmt"
	"os"
	"runtime"

	"research/lmc_ges_tracking/lib"
	"github.com/phil-mansfield/guppy/lib/catio"
)

var (
	MaxSnap = 235
	TreeFileNames = []string{ "path/to/file" }
	CentralRootID = 9285216
	SubhaloRootIDs = []int{ 9278166, 9285216 }
	
	RootIDs = append([]int{CentralRootID}, SubhaloRootIDs...)

	OutputFileName = "halo_tracks"
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

	// Loop over different tree files
	for i := range TreeFileNames {
		rd := catio.TextFile(TreeFileNames[i])

		// If the tree files are too large to load into memory at once,
		// break it up and read each block separately.
		for b := 0; b < rd.Blocks(); b++ {
			runtime.GC()

			// Read the columns
			iCols := rd.ReadIntBlock([]int{1, 28, 31}, b)
			id, dfid, snap := iCols[0], iCols[1], iCols[2]
			
			fCols := rd.ReadFloat32Block([]int{11, 17, 18, 19}, b)
			rvir, x, y, z := fCols[0], fCols[1], fCols[2], fCols[3]

			h := &Haloes{ id, dfid, snap, x, y, z, rvir }

			// Find tracks from the currently-read block
			tracks = append(tracks, GetTracks(h, RootIDs)...)
		}
	}

	WriteTracks(tracks, CentralRootID, OutputFileName)
}

// GetTracks returns all the tracks for haloes whose roots are within rootIDs.
func GetTracks(h *Haloes, rootIDs []int) []HaloTrack {
	tracks := []HaloTrack{ }
	
	order := IDOrder(h.DFID)

	for j := range order {
		i := order[j]

		for k := range rootIDs {
			if h.DFID[i] == rootIDs[k] {
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
	
	for j := j0; ; j++ {
		i := order[j]
		snap := h.Snap[i]

		track.Z[snap] = h.Z[i]
		track.Y[snap] = h.Y[i]
		track.Z[snap] = h.Z[i]
		track.Rvir[snap] = h.Rvir[i]

		if j == len(order) - 1 { break }
		iNext := order[j+1]
		if h.DFID[i]+1 != h.DFID[iNext] || h.Snap[i] <= h.Snap[iNext] { break }
	}

	return track
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
