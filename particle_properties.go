package main

import (
	"encoding/binary"
	"fmt"
	"log"
	"math"
	"os"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"io/ioutil"

	"github.com/phil-mansfield/gravitree"
	"github.com/phil-mansfield/lmc_ges_tracking/lib"
)


var (
	HRLevel = 1
)

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

	eps, mp, blocks, snapFmt, mergerName, idFmt, outFmt :=
		ParseInputFile(inputName)
	for i := range eps {
		if haloIndex != -1 && haloIndex != i { continue }
		AnalyzeHalo(eps[i], mp[i], blocks[i], snapFmt[i],
			mergerName[i], idFmt[i], outFmt[i])
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


func AnalyzeHalo(eps, mp float64, blocks int,
	snapFmt, mergerName, idFmt, outFmt string) {
	
	mergers := lib.ReadMergers(mergerName)
	nh, ns := len(mergers.Mvir), len(mergers.Mvir[0])

	ids, snaps := make([][]int32, nh), make([][]int16, nh)
	for i := 0; i < nh; i++ {
		ids[i], snaps[i] = ReadIDs(fmt.Sprintf(idFmt, i))
	}

	x, v := make([][][3]float32, nh), make([][][3]float32, nh)
	phi := make([][]float32, nh)
	for i := 0; i < nh; i++ {
		x[i] = make([][3]float32, len(ids[i]))
		v[i] = make([][3]float32, len(ids[i]))
		phi[i] = make([]float32, len(ids[i]))
	}

	RvirMW0 := float64(mergers.Mvir[0][ns-1])
	MvirMW0 := float64(mergers.Rvir[0][ns-1])

	for snap := 0; snap < ns; snap++ {
		log.Println("part", snap)
		runtime.GC()

		for i := 0; i < nh; i++ {	
			for j := range x[i] {
				x[i][j] = [3]float32{ -1, -1, -1 }
				v[i][j] = [3]float32{ -1, -1, -1 }
				phi[i][j] = -1
			}
		}
		
		var (
			scale float64
			xp, vp [][3]float32
			idp []int32
		)
		for b := 0; b < blocks; b++ {
			fileName := fmt.Sprintf(snapFmt, snap, b)
			xp, vp, idp, scale = ReadLevel(fileName, HRLevel)

			order := IDOrder(idp)

			for i := 0; i < nh; i++ {
				SetXV(snap, xp, vp, idp, order,
					ids[i], snaps[i], x[i], v[i])
			}
		}
		
		for i := 0; i < nh; i++ {
			CalcPotential(eps, snap, scale, snaps[i],
				x[i], phi[i], RvirMW0, MvirMW0/mp)
			WriteParticles(fmt.Sprintf(outFmt, snap, i),
				x[i], v[i], phi[i], len(ids[i]))
		}
	}
}

func SetXV(snap int, xp, vp [][3]float32, idp []int32,
	order []int, id []int32, snaps []int16, x, v [][3]float32) {

	snap16 := int16(snap)
	for i := range id {
		if snaps[i] > snap16 { continue }
		j := sort.Search(len(idp), func (j int) bool {
			return idp[order[j]] >= id[i]
		})

		if j == len(order) || idp[order[j]] != id[i] { continue }
		x[i] = xp[order[j]]
		v[i] = vp[order[j]]
	}
}

func CalcPotential(eps float64, snap int, scale float64,
	snaps []int16, x [][3]float32, phi []float32, rvir, nvir float64) {

	x0, snap16 := [][3]float64{ }, int16(snap)
	for i := range snaps {
		if snaps[i] <= snap16 {
			x0 = append(x0, [3]float64{
				float64(x[i][0])*scale,
				float64(x[i][1])*scale,
				float64(x[i][2])*scale})
		}
	}
	if len(x0) == 0 { return }

	phi0 := make([]float64, len(x0))
	ComputePhi(eps, x0, phi0, rvir, nvir)
	
	j := 0
	for i := range phi0 {
		if snap16 >= snaps[i] {
			phi[i] = float32(phi0[j])
			j++
		}
	}
}


type Vec64Sort [][3]float64
func (v Vec64Sort) Len() int { return len(v) }
func (v Vec64Sort) Less(i, j int) bool { return v[i][2] < v[j][2] }
func (v Vec64Sort) Swap(i, j int) { v[i], v[j] = v[j], v[i] }


func ComputePhi(eps float64, x [][3]float64, phi []float64,
	rvir, nvir float64) {

    cm := CenterOfMass(x)

	scale := rvir
    for i := range x {
        for k := 0; k < 3; k++ {
            x[i][k] -= cm[k]
            x[i][k] *= scale
        }
    }

    tree := gravitree.NewTree(x)
    tree.Potential(eps * scale, phi)

    for i := range phi {
		phi[i] /= nvir
    }

	r := make([]float64, len(x))
	for i := range r {
		dx, dy, dz := x[i][0], x[i][1], x[i][2]
		r[i] = math.Sqrt(dx*dx + dy*dy + dz*dz)
	}
}

func CenterOfMass(x [][3]float64) [3]float64 {
    sum := [3]float64{ }

    for i := range x {
		for k := 0; k < 3; k++ {
            sum[k] += x[i][k]
        }
    }
    for k := 0; k < 3; k++ {
        sum[k] /= float64(len(x))
    }
    return sum
}


func IDOrder(id []int32) []int {
	fid := make([]float64, len(id))
	for i := range id {
		fid[i] = float64(id[i])
	}

	return lib.QuickSortIndex(fid)
}

func ReadLevel(fileName string, level int) (xp, vp [][3]float32, idp []int32, scale float64) {
	f := lib.OpenGadget2Zoom(fileName, []string{"x", "v", "id32"})
	xp = make([][3]float32, f.N[level])
	vp = make([][3]float32, f.N[level])
	idp = make([]int32, f.N[level])

	f.Read("x", level, xp)
	f.Read("v", level, vp)
	f.Read("id32", level, idp)

	return xp, vp, idp, 1/(1 + f.Z)
}

func ReadIDs(fileName string) (ids []int32, snaps []int16) {
	f, err := os.Open(fileName)
	if err != nil { panic(err.Error()) }

	n := int32(0)
	err = binary.Read(f, binary.LittleEndian, &n)
	if err != nil { panic(err.Error()) }

	ids, snaps = make([]int32, n), make([]int16, n)
	err = binary.Read(f, binary.LittleEndian, ids)
	if err != nil { panic(err.Error()) }
	err = binary.Read(f, binary.LittleEndian, snaps)
	if err != nil { panic(err.Error()) }

	return ids, snaps
}

func WriteParticles(fileName string, x, v [][3]float32,
	phi []float32, nPart int) {

	buf := make([]float32, len(x)*3)

	f, err := os.Create(fileName)
	if err != nil { panic(err.Error()) }
	defer f.Close()

	err = binary.Write(f, binary.LittleEndian, int32(nPart))
	if err != nil { panic(err.Error()) }

	for i := range x {
		buf[3*i]= x[i][0]
		buf[3*i+1]= x[i][1]
		buf[3*i+2]= x[i][2]
	}
	err = binary.Write(f, binary.LittleEndian, buf)
	if err != nil { panic(err.Error()) }

	for i := range v {
		buf[3*i]= v[i][0]
		buf[3*i+1]= v[i][1]
		buf[3*i+2]= v[i][2]
	}
	err = binary.Write(f, binary.LittleEndian, buf)
	if err != nil { panic(err.Error()) }

	err = binary.Write(f, binary.LittleEndian, phi)
	if err != nil { panic(err.Error()) }
}
