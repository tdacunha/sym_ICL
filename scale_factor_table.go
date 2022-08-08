package main

import (
	"encoding/binary"
	"fmt"
	"os"
	"path"
	"strconv"
	"github.com/phil-mansfield/symphony_pipeline/lib"
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
		ScaleFactorTable(cfg, i)
	}
}

func ScaleFactorTable(cfg *lib.Config, cfgi int) {
	fmt.Println(cfgi)
	baseDir, snapFmt := cfg.BaseDir[cfgi], cfg.SnapFormat[cfgi]
	mergers := lib.ReadMergers(lib.MergerFileName(baseDir))

	maxSnap := len(mergers.Rvir[0]) - 1

	outFile := path.Join(baseDir, "halos", "snap_scale.dat")

	scales := make([]float64, maxSnap + 1)

	for snap := 0; snap <= maxSnap; snap++ {
		fileName := fmt.Sprintf(snapFmt, snap, 0)
		header := lib.OpenGadget2Zoom(fileName, []string{"x", "v", "id32"})
		scales[snap] = header.Scale
	}

	f, err := os.Create(outFile)
	if err != nil { panic(err.Error()) }
	defer f.Close()

	err = binary.Write(f, binary.LittleEndian, int64(len(scales)))
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, binary.LittleEndian, scales)
	if err != nil { panic(err.Error()) }
}
