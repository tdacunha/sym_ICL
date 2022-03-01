package lib

import (
	"math"
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
	
	if err != nil { panic(err.Error()) }
	err = binary.Read(fHeader, order, hd.FileLengths)
	if err != nil { panic(err.Error()) }
	err = binary.Read(fHeader, order, hd.Offsets)
	if err != nil { panic(err.Error()) }
	err = binary.Read(fHeader, order, hd.Sizes)
	if err != nil { panic(err.Error()) }
	err = binary.Read(fHeader, order, hd.FileIdxs)
	if err != nil { panic(err.Error()) }

	return hd
}

// Tags is a collection of informaiton
type Tags struct {
	Idx [][]int32
	Snap [][]int16
	Flag [][]uint8
}

// ReadAllTags reads the tags associated with all tracked haloes.
func ReadTags(baseDir string, hd *ParticleHeader) *Tags {
	order := binary.LittleEndian
	i := 0 // indexes over all haloes

	tags := &Tags{
		Idx: make([][]int32, hd.NHalo),
		Snap: make([][]int16, hd.NHalo),
		Flag: make([][]uint8, hd.NHalo),
	}
	
	for iFile := 0; iFile < int(hd.NFile); iFile++ {
		fIdx, err := os.Open(TagFileName(baseDir, iFile))
		if err != nil { panic(err.Error()) }
		fSnap, err := os.Open(SnapFileName(baseDir, iFile))
		if err != nil { panic(err.Error()) }
		fFlag, err := os.Open(FlagFileName(baseDir, iFile))	
		if err != nil { panic(err.Error()) }

		idx := make([]int32, hd.FileLengths[iFile])
		snap := make([]int16, hd.FileLengths[iFile])
		flag := make([]uint8, hd.FileLengths[iFile])
		
		err = binary.Read(fIdx, order, idx)
		if err != nil { panic(err.Error()) }
		err = binary.Read(fSnap, order, snap)
		if err != nil { panic(err.Error()) }
		err = binary.Read(fFlag, order, flag)
		if err != nil { panic(err.Error()) }
		
		fIdx.Close()
		fSnap.Close()
		fFlag.Close()

		for ; i < len(hd.FileIdxs) &&  hd.FileIdxs[i] == int32(iFile); i++ {
			start, end := hd.Offsets[i], hd.Offsets[i] + hd.Sizes[i]
			tags.Idx[i] = idx[start: end]
			tags.Snap[i] = snap[start: end]
			tags.Flag[i] = flag[start: end]
		}
	}

	return tags
}

// WriteTags writes tags and other particle information to nFiles files to the
// correct location in baseDir. Each element in teh three arrays correspond to
// a different halo in the mergers.dat file. idxs is a an array of indices to
// the particles in the halo. snaps is an array of the first snapshot that those
// particles entered that halo, and flags is flags containing information about
// the particle.
func WriteTags(baseDir string, nFiles int, tags *Tags) {
	MaybeMkdir(ParticleDirName(baseDir))
	MaybeMkdir(TagDirName(baseDir))
	
	idxs, snaps, flags := tags.Idx, tags.Snap, tags.Flag
	
	order := binary.LittleEndian
	
	// Find target number of particles per file.
	nTot := 0
	for i := range idxs { nTot += len(idxs[i]) }
	nPerFile := int(math.Ceil(float64(nTot)/float64(nFiles)))
	
	fileLengths := []int32{ }
	offsets, sizes := make([]int32, len(idxs)), make([]int32, len(idxs))
	fileIdxs := make([]int32, len(idxs))

	j := 0 // Indexes over haloes
	totalLen := 0 // Total number of particles
	for iFile := 0; iFile < nFiles; iFile++ {
		fileLen := 0 // number of tags per file

		// Create files
		fTag, err := os.Create(TagFileName(baseDir, iFile))
		if err != nil { panic(err.Error()) }
		fSnap, err := os.Create(SnapFileName(baseDir, iFile))
		if err != nil { panic(err.Error()) }
		fFlag, err := os.Create(FlagFileName(baseDir, iFile))	
		if err != nil { panic(err.Error()) }
		
		for ; j < len(idxs) && totalLen < (iFile+1)*nPerFile; j++{
			err = binary.Write(fTag, order, idxs[j])
			if err != nil { panic(err.Error()) }
			err = binary.Write(fSnap, order, snaps[j])
			if err != nil { panic(err.Error()) }
			err = binary.Write(fFlag, order, flags[j])
			if err != nil { panic(err.Error()) }
			
			// Log information into particle header arrays
			offsets[j] = int32(fileLen)
			sizes[j] = int32(len(idxs[j]))
			fileIdxs[j] = int32(iFile)

			fileLen += len(idxs[j])
			totalLen += len(idxs[j])
		}

		fileLengths = append(fileLengths, int32(fileLen))

		fTag.Close()
		fSnap.Close()
		fFlag.Close()
	}

	fHeader, err := os.Create(ParticleHeaderName(baseDir))
	if err != nil { panic(err.Error()) }
	defer fHeader.Close()

	err = binary.Write(fHeader, order, int32(nFiles))
	if err != nil { panic(err.Error()) }
	err = binary.Write(fHeader, order, int32(len(idxs)))
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
	baseDir string, nFiles int, varName string, snap int, x [][][3]float32,
) {
	MaybeMkdir(SnapDirName(baseDir, snap))

	order := binary.LittleEndian

	nTot := 0
	for i := range x { nTot += len(x[i]) }
	nPerFile := int(math.Ceil(float64(nTot)/float64(nFiles)))
	
	j := 0 // Indexes over haloes
	totalLen := 0 // Total number of particles
	for iFile := 0; iFile < nFiles; iFile++ {
		// Create files
		f, err := os.Create(VarFileName(baseDir, varName, snap, iFile))
		if err != nil { panic(err.Error()) }
		
		err = binary.Write(f, order, int32(VectorVarCode))
		if err != nil { panic(err.Error()) }
		
		for ; j < len(x) && totalLen < (iFile+1)*nPerFile; j++{
			err = binary.Write(f, order, x[j])
			if err != nil { panic(err.Error()) }

			totalLen += len(x[j])
		}

		f.Close()
	}
}

// ReadVector reads the floats associated with a given halo.
func ReadVector(baseDir, varName string,
	snap int, hd *ParticleHeader) [][][3]float32 {
	
	order := binary.LittleEndian
	i := 0 // indexes over all haloes

	x := make([][][3]float32, hd.NHalo)
	
	for iFile := 0; iFile < int(hd.NFile); iFile++ {
		f, err := os.Open(VarFileName(baseDir, varName, snap, iFile))
		if err != nil { panic(err.Error()) }

		code := int32(-1)
		err = binary.Read(f, order, &code)
		if code != VectorVarCode {
			panic(fmt.Sprintf("%s is not vector variable file.",
				VarFileName(baseDir, varName, snap, i)))
		}

		
		xi := make([][3]float32, hd.FileLengths[iFile])
		
		err = binary.Read(f, order, xi)
		if err != nil { panic(err.Error()) }
		
		f.Close()

		for ; i < len(hd.FileIdxs) &&  hd.FileIdxs[i] == int32(iFile); i++ {
			start, end := hd.Offsets[i], hd.Offsets[i] + hd.Sizes[i]
			x[i] = xi[start: end]
		}
	}

	return x
}

// WriteFloat splits a float32 array across files in the same way that WriteTags
// splits tags across files.
func WriteFloat(
	baseDir string, nFiles int, varName string, snap int, x [][]float32,
) {
	MaybeMkdir(SnapDirName(baseDir, snap))

	order := binary.LittleEndian

	nTot := 0
	for i := range x { nTot += len(x[i]) }
	nPerFile := int(math.Ceil(float64(nTot)/float64(nFiles)))
	
	j := 0 // Indexes over haloes
	totalLen := 0 // Total number of particles
	for iFile := 0; iFile < nFiles; iFile++ {
		// Create files
		f, err := os.Create(VarFileName(baseDir, varName, snap, iFile))
		if err != nil { panic(err.Error()) }
		
		err = binary.Write(f, order, int32(FloatVarCode))
		if err != nil { panic(err.Error()) }
		
		for ; j < len(x) && totalLen < (iFile+1)*nPerFile; j++{
			err = binary.Write(f, order, x[j])
			if err != nil { panic(err.Error()) }

			totalLen += len(x[j])
		}

		f.Close()
	}
}

// ReadFloat reads the floats associated with a given halo.
func ReadFloat(baseDir string, varName string,
	snap int, hd *ParticleHeader) [][]float32 {
	
	order := binary.LittleEndian
	i := 0 // indexes over all haloes

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

		xi := make([]float32, hd.FileLengths[iFile])
		
		err = binary.Read(f, order, xi)
		if err != nil { panic(err.Error()) }
		
		f.Close()

		for ; i < len(hd.FileIdxs) &&  hd.FileIdxs[i] == int32(iFile); i++ {
			start, end := hd.Offsets[i], hd.Offsets[i] + hd.Sizes[i]
			x[i] = xi[start: end]
		}
	}

	return x
}
