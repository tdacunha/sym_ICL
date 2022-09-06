package main

import (
	//"fmt"
	"log"
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
		Phi(cfg, i)
		runtime.GC()
	}
}

func Phi(cfg *lib.Config, cfgi int) {
	log.Printf("Computing potential energies for halo %d (%s)",
		cfgi, cfg.BaseDir[cfgi])

	baseDir, snapFmt := cfg.BaseDir[cfgi], cfg.SnapFormat[cfgi]

	mergers := lib.ReadMergers(lib.MergerFileName(baseDir))
	pHeader := lib.ReadParticleHeader(cfg.BaseDir[cfgi])
	tags := lib.ReadTags(cfg.BaseDir[cfgi], pHeader)

	_, _ = snapFmt, tags
	maxSnap := len(mergers.Rvir[0])
	for snap := 0; snap <= maxSnap; snap++ {
		log.Printf("Snapshot %03d", snap)
		xh := lib.ReadVector(baseDir, "x", snap, pHeader)
		vh := lib.ReadVector(baseDir, "v", snap, pHeader)
		_ = vh
		for i := range xh {
			if len(xh[i]) > 0 {
				log.Printf("    %d: %d", i, len(xh[i]))
			}
		}
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
