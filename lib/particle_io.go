package lib

import (
	"math"
	"math/rand"
	"os"
	"fmt"
	"encoding/binary"
)

const (
	VectorVarCode = iota
	FloatVarCode
)

// ParticleHeader is the contents of the header.dat file, and contains
// information on where a halo's particles are within the different files.
type ParticleHeader struct {
	NFile int32 // The numbers of files
	NHalo, NParticle int32 // Number of haloes and particles across all files
	FileLengths []int32 // The number of particles per file
	Offsets []int32 // Starting index of a given halo in its file
	Sizes []int32 // Number of particles associated with the halo
	FileIdxs []int32 // Which file the halo is in
	N0 []int32 // The number of flag=0 particles in the halo.
}

// ReadParticleHeader reads the header file containng the number of files,
// haloes, and particles, as well as where those haloes are within files and
// how many particles they have.
func ReadParticleHeader(baseDir string) *ParticleHeader {
	order := binary.LittleEndian
	
	fHeader, err := os.Open(ParticleHeaderName(baseDir))
	if err != nil { panic(err.Error()) }
	defer fHeader.Close()

	hd := &ParticleHeader{ }
	
	err = binary.Read(fHeader, order, &hd.NFile)
	if err != nil { panic(err.Error()) }
	err = binary.Read(fHeader, order, &hd.NHalo)
	if err != nil { panic(err.Error()) }
	err = binary.Read(fHeader, order, &hd.NParticle)
	if err != nil { panic(err.Error()) }

	hd.FileLengths = make([]int32, hd.NFile)
	hd.Offsets = make([]int32, hd.NHalo)
	hd.Sizes = make([]int32, hd.NHalo)
	hd.FileIdxs = make([]int32, hd.NHalo)
	hd.N0 = make([]int32, hd.NHalo)

	if err != nil { panic(err.Error()) }
	err = binary.Read(fHeader, order, hd.FileLengths)
	if err != nil { panic(err.Error()) }
	err = binary.Read(fHeader, order, hd.Offsets)
	if err != nil { panic(err.Error()) }
	err = binary.Read(fHeader, order, hd.Sizes)
	if err != nil { panic(err.Error()) }
	err = binary.Read(fHeader, order, hd.FileIdxs)
	if err != nil { panic(err.Error()) }
	err = binary.Read(fHeader, order, hd.N0)

	return hd
}

// Tags is a collection of informaiton
type Tags struct {
	N0 []int32 // Number of Flag = 0 particles for each halo
	ID [][]int32 // ID of each particles
	Snap [][]int16 // Snapshot of each particle
	Flag [][]uint8 // 0 -> belonged to this halo first, 1 ->
	               // belonged to a less massive halo first
}

// ReadAllTags reads the tags associated with all tracked haloes.
func ReadTags(baseDir string, hd *ParticleHeader) *Tags {
	order := binary.LittleEndian
	i := 0 // indexes over all haloes

	tags := &Tags{
		ID: make([][]int32, hd.NHalo),
		Snap: make([][]int16, hd.NHalo),
		Flag: make([][]uint8, hd.NHalo),
	}
	
	for iFile := 0; iFile < int(hd.NFile); iFile++ {
		f, err := os.Open(TagFileName(baseDir, iFile))
		if err != nil { panic(err.Error()) }

		for ; i < len(hd.FileIdxs) &&  hd.FileIdxs[i] == int32(iFile); i++ {
			n := hd.Sizes[i]

			id := make([]int32, n)
			snap := make([]int16, n)
			flag := make([]uint8, n)

			err = binary.Read(f, order, id)
			if err != nil { panic(err.Error()) }
			err = binary.Read(f, order, snap)
			if err != nil { panic(err.Error()) }
			err = binary.Read(f, order, flag)
			if err != nil { panic(err.Error()) }

			tags.ID[i] = id
			tags.Snap[i] = snap
			tags.Flag[i] = flag
		}

		f.Close()
	}

	tags.N0 = hd.N0

	return tags
}

func ReadTagLookup(baseDir string) *TagLookup {
	f, err := os.Open(TagLookupName(baseDir))
	if err != nil { panic(err.Error()) }

	order := binary.LittleEndian

	np := int32(0)
	err = binary.Read(f, order, &np)
	if err != nil { panic(err.Error()) }

	l := &TagLookup{
		Halo: make([]int16, np),
		Index: make([]int32, np),
	}

	err = binary.Read(f, order, l.Halo)
	if err != nil { panic(err.Error()) }
	err = binary.Read(f, order, l.Index)
	if err != nil { panic(err.Error()) }

	return l
}

// WriteTags writes tags and other particle information to nFiles files to the
// correct location in baseDir. Each element in the three arrays correspond to
// a different halo in the subhalos.dat file. idxs is a an array of indices to
// the particles in the halo. snaps is an array of the first snapshot that those
// particles entered that halo, and flags is flags containing information about
// the particle.
func WriteTags(baseDir string, nFiles int, tags *Tags, l *TagLookup) {
	MaybeMkdir(ParticleDirName(baseDir))
	
	ids, snaps, flags := tags.ID, tags.Snap, tags.Flag
	
	order := binary.LittleEndian
	
	// Find target number of particles per file.
	nTot := 0
	for i := range ids { nTot += len(ids[i]) }
	nPerFile := int(math.Ceil(float64(nTot)/float64(nFiles)))
	
	fileLengths := []int32{ }
	offsets, sizes := make([]int32, len(ids)), make([]int32, len(ids))
	fileIdxs := make([]int32, len(ids))

	j := 0 // Indexes over haloes
	totalLen := 0 // Total number of particles
	for iFile := 0; iFile < nFiles; iFile++ {
		fileLen := 0 // number of tags per file

		f, err := os.Create(TagFileName(baseDir, iFile))
		
		for ; j < len(ids) && totalLen < (iFile+1)*nPerFile; j++{
			err = binary.Write(f, order, ids[j])
			if err != nil { panic(err.Error()) }
			err = binary.Write(f, order, snaps[j])
			if err != nil { panic(err.Error()) }
			err = binary.Write(f, order, flags[j])
			if err != nil { panic(err.Error()) }
			
			// Log information into particle header arrays
			offsets[j] = int32(fileLen)
			sizes[j] = int32(len(ids[j]))
			fileIdxs[j] = int32(iFile)

			fileLen += len(ids[j])
			totalLen += len(ids[j])
		}

		fileLengths = append(fileLengths, int32(fileLen))

		f.Close()
	}

	fHeader, err := os.Create(ParticleHeaderName(baseDir))
	if err != nil { panic(err.Error()) }

	err = binary.Write(fHeader, order, int32(nFiles))
	if err != nil { panic(err.Error()) }
	err = binary.Write(fHeader, order, int32(len(ids)))
	if err != nil { panic(err.Error()) }
	err = binary.Write(fHeader, order, int32(totalLen))
	if err != nil { panic(err.Error()) }
	err = binary.Write(fHeader, order, fileLengths)
	if err != nil { panic(err.Error()) }
	err = binary.Write(fHeader, order, offsets)
	if err != nil { panic(err.Error()) }
	err = binary.Write(fHeader, order, sizes)
	if err != nil { panic(err.Error()) }
	err = binary.Write(fHeader, order, fileIdxs)
	if err != nil { panic(err.Error()) }
	err = binary.Write(fHeader, order, tags.N0)
	if err != nil { panic(err.Error()) }

	fHeader.Close()

	f, err := os.Create(TagLookupName(baseDir))
	if err != nil { panic(err.Error()) }

	err = binary.Write(f, order, int32(len(l.Halo)))
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, l.Halo)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, l.Index)
	if err != nil { panic(err.Error()) }

	f.Close()
}

func fileExists(name string) bool {
	// Based on this stack overflow answe
	// https://stackoverflow.com/questions/10510691/how-to-check-whether-a-file-or-directory-exists
	_, err := os.Stat(name)
	if err == nil { return true }
	if os.IsNotExist(err) { return false }
	panic(err.Error())
}

// Makes a directory if it doesn't already exist
func MaybeMkdir(name string) {
	if fileExists(name) { return }
	err := os.MkdirAll(name, 0744)
	if err != nil { panic(err.Error()) }
}

// WriteVector splits a vector array across files in the same way that WriteTags
// splits tags across files.
func WriteVector(
	hd *ParticleHeader, baseDir, varName string, snap int, x [][][3]float32,
) {
	MaybeMkdir(SnapDirName(baseDir, snap))

	order := binary.LittleEndian
	
	totalLen := 0 // Total number of particles
	for iFile := 0; iFile < int(hd.NFile); iFile++ {
		// Create files
		f, err := os.Create(VarFileName(baseDir, varName, snap, iFile))
		if err != nil { panic(err.Error()) }
		
		err = binary.Write(f, order, int32(VectorVarCode))
		if err != nil { panic(err.Error()) }
		
		for j := 0; j < int(hd.NHalo); j++ {
			if int32(iFile) != hd.FileIdxs[j] || len(x[j]) == 0 { continue }

			x16, min, max := Vector32ToUint16(x[j])

			err = binary.Write(f, order, int64(len(x[j])))
			if err != nil { panic(err.Error()) }
			err = binary.Write(f, order, min)
			if err != nil { panic(err.Error()) }
			err = binary.Write(f, order, max)
			if err != nil { panic(err.Error()) }
			err = binary.Write(f, order, x16)
			if err != nil { panic(err.Error()) }

			totalLen += len(x[j])
		}

		f.Close()
	}
}

// ReadVector reads the floats associated with a given halo.
func ReadVector(
	hd *ParticleHeader, baseDir, varName string, snap int,
) [][][3]float32 {
	
	order := binary.LittleEndian

	x := make([][][3]float32, hd.NHalo)
	
	for iFile := 0; iFile < int(hd.NFile); iFile++ {
		f, err := os.Open(VarFileName(baseDir, varName, snap, iFile))
		if err != nil { panic(err.Error()) }

		code := int32(-1)
		err = binary.Read(f, order, &code)
		if code != VectorVarCode {
			panic(fmt.Sprintf("%s is not vector variable file.",
				VarFileName(baseDir, varName, snap, iFile)))
		}

		for j := 0; j < int(hd.NHalo); j++ {
			if int32(iFile) != hd.FileIdxs[j] {
				continue
			} else if hd.N0[j] == 0 {
				x[j] = [][3]float32{ }
				continue
			}
			//panic("Conditioning in loop not updated")

			var size int64
			err = binary.Read(f, order, &size)
			if err != nil { panic(err.Error()) }

			xi16 := make([]uint16, 3*size)
			xf32 := make([][3]float32, size)
			var min, max [3]float32

			err = binary.Read(f, order, &min)
			if err != nil { panic(err.Error()) }
			err = binary.Read(f, order, &max)
			if err != nil { panic(err.Error()) }
			err = binary.Read(f, order, xi16)
			if err != nil { panic(err.Error()) }

			Uint16ToVector32(xi16, min, max, xf32)
			x[j] = xf32
		}
				
		f.Close()
	}

	return x
}

// WriteFloat splits a float32 array across files in the same way that WriteTags
// splits tags across files.
func WriteFloat(
	hd *ParticleHeader, baseDir, varName string, snap int, x [][]float32,
) {
	MaybeMkdir(SnapDirName(baseDir, snap))

	order := binary.LittleEndian

	totalLen := 0 // Total number of particles
	for iFile := 0; iFile < int(hd.NFile); iFile++ {
		// Create files
		f, err := os.Create(VarFileName(baseDir, varName, snap, iFile))
		if err != nil { panic(err.Error()) }
		
		err = binary.Write(f, order, int32(FloatVarCode))
		if err != nil { panic(err.Error()) }

		for j := 0; j < int(hd.NHalo); j++ {
			if int32(iFile) != hd.FileIdxs[j] || len(x[j]) == 0 {
				continue
			} 
			
			x16, min, max := Float32ToUint16(x[j])

			err = binary.Write(f, order, int64(len(x[j])))
			if err != nil { panic(err.Error()) }
			err = binary.Write(f, order, min)
			if err != nil { panic(err.Error()) }
			err = binary.Write(f, order, max)
			if err != nil { panic(err.Error()) }
			err = binary.Write(f, order, x16)
			if err != nil { panic(err.Error()) }

			totalLen += len(x[j])
		}

		f.Close()
	}
}

func Float32ToUint16(x []float32) (out []uint16, min, max float32) {
	min, max = x[0], x[0]
	out = make([]uint16, len(x))

	for i := 1; i < len(x); i++ {
		if x[i] < min {
			min = x[i]
		} else if x[i] > max {
			max = x[i]
		}
	}

	dx := max - min
	if dx == 0 { dx = 1 }
	for i := range x {
		ix := int(math.MaxUint16 * (x[i] - min) / dx)
		if ix < 0 { ix = 0 }
		if ix > math.MaxUint16 { ix = math.MaxUint16 }
		out[i] = uint16(ix)
	}

	return out, min, max
}

func Vector32ToUint16(x [][3]float32) (out []uint16, min, max [3]float32) {
	min, max = x[0], x[0]
	out = make([]uint16, len(x)*3)

	for i := 1; i < len(x); i++ {
		for dim := 0; dim < 3; dim++ {
			if x[i][dim] < min[dim] {
				min[dim] = x[i][dim]
			} else if x[i][dim] > max[dim] {
				max[dim] = x[i][dim]
			}
		}
	}

	dx := [3]float32{ }
	for dim := 0; dim < 3; dim++ {
		dx[dim] = max[dim] - min[dim]
		if dx[dim] == 0 { dx[dim] = 1 }
	}
	for i := range x {
		for dim := 0; dim < 3; dim++ {
			ix := int(math.MaxUint16 * (x[i][dim] - min[dim]) / dx[dim])
			if ix < 0 { ix = 0 }
			if ix > math.MaxUint16 { ix = math.MaxUint16 }
			out[3*i+dim] = uint16(ix)
		}
	}
	
	return out, min, max
}

func Uint16ToFloat32(x []uint16, min, max float32, out []float32) {	
	dx := max - min
	for i := range out {
		out[i] = dx/math.MaxUint16*(float32(x[i]) + rand.Float32()) + min
	}
}

func Uint16ToVector32(x []uint16, min, max [3]float32, out [][3]float32) {	
	dx := [3]float32{ }
	for dim := 0; dim < 3; dim++ {
		dx[dim] = max[dim] - min[dim]
	}
	for i := range out {
		out[i] = [3]float32{
			dx[0]/math.MaxUint16*(float32(x[3*i+0]) + rand.Float32()) + min[0],
			dx[1]/math.MaxUint16*(float32(x[3*i+1]) + rand.Float32()) + min[1],
			dx[2]/math.MaxUint16*(float32(x[3*i+2]) + rand.Float32()) + min[2],
		}
	}
}

// ReadFloat reads the floats associated with a given halo.
func ReadFloat(
	hd * ParticleHeader, baseDir, varName string, snap int,
) [][]float32 {
	
	order := binary.LittleEndian

	x := make([][]float32, hd.NHalo)
	
	for iFile := 0; iFile < int(hd.NFile); iFile++ {
		f, err := os.Open(VarFileName(baseDir, varName, snap, iFile))
		if err != nil { panic(err.Error()) }

		code := int32(-1)
		err = binary.Read(f, order, &code)
		if code != FloatVarCode {
			panic(fmt.Sprintf("%s is not a float variable file.",
				VarFileName(baseDir, varName, snap, iFile)))
		}

		for j := 0; j < int(hd.NHalo); j++ {
			if int32(iFile) != hd.FileIdxs[j] {
				continue
			} else if hd.N0[j] == 0 { 
				x[j] = []float32{ }
				continue
			}

			var size int64
			err = binary.Read(f, order, &size)
			if err != nil { panic(err.Error()) }

			xi16 := make([]uint16, size)
			xf32 := make([]float32, size)
			var min, max float32

			err = binary.Read(f, order, &min)
			if err != nil { panic(err.Error()) }
			err = binary.Read(f, order, &max)
			if err != nil { panic(err.Error()) }
			err = binary.Read(f, order, xi16)
			if err != nil { panic(err.Error()) }

			Uint16ToFloat32(xi16, min, max, xf32)
			x[j] = xf32
		}

		f.Close()
	}

	return x
}
