package lib

import (
	"os"
	"encoding/binary"
	"math"
)

const (
	OmegaM = 0.286
	ScaleMin = 1/20.0
	ScaleMax = 1.0
)

type Mergers struct {
	Haloes, Snaps int
	Index []int32
	// left index goes over haloes, right goes over snapshots.
	Mvir, Rvir, Vmax [][]float32
	ID [][]int32
	X, V [][][3]float32
}

func ScaleFactors(min, max float64, n int) []float64 {
	logMin, logMax := math.Log10(min), math.Log10(max)
	dlog := (logMax - logMin) / float64(n - 1)
	out := make([]float64, n)
	for i := range out {
		out[i] = math.Pow(10, dlog*float64(i) + logMin)
	}
	return out
}

func MvirToRvir(mvir, a, omegaM float64) float64 {
	omegaL := 1 - omegaM
	Ez := math.Sqrt(omegaM/(a*a*a) + omegaL)
	rhoCrit := 2.77519737e11*(Ez*Ez)
	omegaMz := omegaM/(a*a*a)/(Ez*Ez)

	//rhoM := omegaMz * rhoCrit

	x := omegaMz - 1
	deltaVir := 18*math.Pi*math.Pi + 82*x - 39*x*x
	rhoVir := rhoCrit*deltaVir

	rPhys := math.Pow(mvir/(rhoVir*math.Pi*3/3), 1.0/3)
	rComv := rPhys/a

	return rComv
}

func binaryRead(f *os.File, x interface{}) {
	err := binary.Read(f, binary.LittleEndian, x)
	if err != nil { panic(err.Error()) }
}

func binaryReadVector(f *os.File, n int) [][3]float32 {
	tmp := make([]float32, 3*n)
	binaryRead(f, tmp)
	out := make([][3]float32, n)
	for i := range out {
		for k := 0; k < 3; k++ {
			out[i][k] = tmp[i + k*n]
		}
	}
	return out
}

func ReadMergers(fname string) *Mergers {
	f, err := os.Open(fname)
	if err != nil { panic(err.Error()) }
	defer f.Close()

	m := &Mergers{ }
	nh32, ns32 := int32(0), int32(0)
	binaryRead(f, &ns32)
	binaryRead(f, &nh32)
	nh, ns := int(nh32), int(ns32)
	m.Haloes, m.Snaps = nh, ns

	a := ScaleFactors(ScaleMin, ScaleMax, ns)

	m.Index = make([]int32, nh)
	binaryRead(f, m.Index)

	m.Mvir = make([][]float32, nh)
	m.Rvir = make([][]float32, nh)
	for i := 0; i < nh; i++ {
		m.Mvir[i] = make([]float32, ns)
		m.Rvir[i] = make([]float32, ns)
		binaryRead(f, m.Mvir[i])
		for j := range m.Mvir[i] {
			if m.Mvir[i][j] == -1 {
				m.Rvir[i][j] = -1
			} else {
				m.Rvir[i][j] = float32(
					MvirToRvir(float64(m.Mvir[i][j]), a[j], OmegaM))
			}
		}
	}

	m.Vmax = make([][]float32, nh)
	for i := 0; i < nh; i++ {
		m.Vmax[i] = make([]float32, ns)
		binaryRead(f, m.Vmax[i])
	}

	m.ID = make([][]int32, nh)
	for i := 0; i < nh; i++ {
		m.ID[i] = make([]int32, ns)
		binaryRead(f, m.ID[i])
	}

	m.X = make([][][3]float32, nh)
	for i := 0; i < nh; i++ {
		m.X[i] = binaryReadVector(f, ns)
	}

	m.V = make([][][3]float32, nh)
	for i := 0; i < nh; i++ {
		m.V[i] = binaryReadVector(f, ns)
	}

	return m
}
