package main

import (
	"encoding/binary"
	"fmt"
	"os"
	"log"
	
	"github.com/phil-mansfield/guppy/lib/catio"
)

var (
	TreeFileName ="/scratch/users/enadler/Halo416/rockstar/trees/tree_0_0_0.dat"
	OutFileName ="/scratch/users/phil1/lmc_ges_tracking/Halo416/tree.bin"
)

func main() {
	log.Println("Starting")
	cfg := catio.DefaultConfig
	cfg.SkipLines = 45
	rd := catio.TextFile(TreeFileName, cfg)
	
	log.Println("Reading")
	cols := rd.ReadInts([]int{ 1, 3, 6, 28, 31 })
	id, descid, upid, dfid, snap := cols[0], cols[1], cols[2], cols[3], cols[4]
	
	mvir := rd.ReadFloat64s([]int{ 10 })[0]

	log.Println("Writing")
	Write(id, descid, upid, dfid, snap, mvir)
	log.Println("Done")
}

func Write(id, descid, upid, dfid, snap []int, mvir []float64) {
	n := len(id)
	id32, descid32, upid32 := make([]int32,n), make([]int32,n), make([]int32,n)
	dfid32, snap32 := make([]int32, n), make([]int32, n)
	maxID, maxDFID := 0, 0
	for i := range id32 {
		id32[i], descid32[i] = int32(id[i]), int32(descid[i])
		upid32[i] = int32(upid[i])
		dfid32[i], snap32[i] = int32(dfid[i]), int32(snap32[i])

		if id[i] > maxID {
			maxID = id[i]
		}
		if dfid[i] > maxDFID {
			maxDFID = dfid[i]
		}
	}

	fmt.Println(len(id), maxID, maxDFID)

	f, err := os.Create(OutFileName)
	defer f.Close()
	if err != nil { panic(err.Error()) }

	order := binary.LittleEndian
	err = binary.Write(f, order, id32)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, descid32)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, upid32)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, dfid32)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, snap32)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, mvir)
	if err != nil { panic(err.Error()) }
}
