package main

import (
	"fmt"
	"log"
	"math"
	"os"
	"runtime"
	"strconv"
	"github.com/phil-mansfield/symphony_pipeline/lib"
)

var (
	HRLevel = 1
)

func main() {
	haloIndex := -1
    if len(os.Args) > 3 || len(os.Args) <= 1 {
        panic("You must supply a config file and an (optional) index to " +
			"restrict analysis to.")
    } else if len(os.Args) == 3 {
		var err error
		haloIndex, err = strconv.Atoi(os.Args[2])
		if err != nil { panic(err.Error()) }
	}

	inputName := os.Args[1]
	
	cfg := lib.ParseConfig(inputName)
	
	for i := range cfg.BaseDir {
		if haloIndex != -1 && i != haloIndex { continue }
		XV(cfg, i)
		runtime.GC()
	}
}

func XV(cfg *lib.Config, cfgi int) {
	log.Printf("Assigning positions and velocities for halo %d (%s)",
		cfgi, cfg.BaseDir[cfgi])

	baseDir, snapFmt := cfg.BaseDir[cfgi], cfg.SnapFormat[cfgi]
	mergers := lib.ReadMergers(lib.MergerFileName(baseDir))

	maxSnap := len(mergers.Rvir[0]) - 1

	fileName := fmt.Sprintf(snapFmt, maxSnap, 0)
	header := lib.OpenGadget2Zoom(fileName, []string{"x", "v", "id32"})

	x := make([][3]float32, header.NTot[HRLevel])
	v := make([][3]float32, header.NTot[HRLevel])
	
	pHeader := lib.ReadParticleHeader(cfg.BaseDir[cfgi])

	tags := lib.ReadTags(cfg.BaseDir[cfgi], pHeader)
	
	xh := NewHaloVector(tags)
	vh := NewHaloVector(tags)

	for snap := 0; snap <= maxSnap; snap++ {
		log.Printf("Snapshot %3d", snap)
		dirName := lib.SnapDirName(cfg.BaseDir[cfgi], snap)
		lib.MaybeMkdir(dirName)

		ClearVectors(x)
		ClearVectors(v)

		for b := 0; b < cfg.Blocks[cfgi]; b++ {
			fileName := fmt.Sprintf(snapFmt, snap, b)
			ReadToIDGrid(fileName, x, v)
			runtime.GC()
		}

		CheckVectors(x)
		CheckVectors(v)

		FillHaloVector(tags, x, xh, snap)
		FillHaloVector(tags, v, vh, snap)

		lib.WriteVector(pHeader, cfg.BaseDir[cfgi], "x", snap, xh)
		lib.WriteVector(pHeader, cfg.BaseDir[cfgi], "v", snap, vh)

		runtime.GC()
	}
}

func CountFiles(pHeader *lib.ParticleHeader) int {
	max := pHeader.FileIdxs[0]
	for i := range pHeader.FileIdxs {
		if pHeader.FileIdxs[i] > max {
			max = pHeader.FileIdxs[i]
		}
	}
	return int(max)
}

func ReadToIDGrid(fileName string, x, v [][3]float32) {
	f := lib.OpenGadget2Zoom(fileName, []string{"x", "v", "id32"})
	
	xp := make([][3]float32, f.N[HRLevel])
	idp := make([]int32, f.N[HRLevel])

	f.Read("x", HRLevel, xp)
	f.Read("id32", HRLevel, idp)

	for i := range xp {
		x[idp[i]] = xp[i]
	}

	vp := xp
	f.Read("v", HRLevel, vp)
	
	for i := range vp {
		v[idp[i]] = vp[i]
	}
}

func ClearVectors(x [][3]float32) {
	for i := range x {
		for dim := 0; dim < 3; dim++ {
			x[i][dim] = float32(math.NaN())
		}
	}
}

func CheckVectors(x [][3]float32) {
	for i := range x {
		for dim := 0; dim < 3; dim++ {
			if math.IsNaN(float64(x[i][dim])) {
				panic(fmt.Sprintf("Vector with ID %d not set.", i + 1))
			}
		}
	}
}

func NewHaloVector(tags *lib.Tags) [][][3]float32 {
	vec := make([][][3]float32, len(tags.N0))
	for i := range vec {
		vec[i] = make([][3]float32, tags.N0[i])
	}
	return vec
}

func FillHaloVector(tags *lib.Tags, x [][3]float32,
	xh [][][3]float32, snap int) {

	for i := range xh {
		xh[i] = xh[i][:cap(xh[i])]
		for j := range xh[i] {
			xh[i][j] = x[tags.ID[i][j] - 1]
		}

		xh[i] = tags.TrimVector(i, xh[i], snap)
	}
}
