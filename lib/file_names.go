package lib

import (
	"path"
	"fmt"
	"io/ioutil"
)


func ParticleDirName(baseDir string) string {
	return path.Join(baseDir, "particles")
}

func HaloDirName(baseDir string) string {
	return path.Join(baseDir, "haloes")
}

func SnapDirName(baseDir string, snap int) string {
	return path.Join(ParticleDirName(baseDir), fmt.Sprintf("snap_%03d", snap))
}

func ParticleHeaderName(baseDir string) string {
	return path.Join(ParticleDirName(baseDir), "particle_header.dat")
}

func MergerFileName(baseDir string) string {
	return path.Join(HaloDirName(baseDir), "mergers.dat")
}

func BranchesFileName(baseDir string) string {
	return path.Join(HaloDirName(baseDir), "branches.dat")
}

func TreeFileNames(baseDir string) []string {
	dir := HaloDirName(baseDir)
	
	files, err := ioutil.ReadDir(dir)
	if err != nil { panic(err.Error()) }
	
	out := []string{ }
	for i := range files {
		name := files[i].Name()
		if len(name) >= 7 &&
			name[len(name)-7:] == ".df.bin" {
			out = append(out, path.Join(dir, name))
		}
	}

	return out
}

func TagDirName(baseDir string) string {
	return path.Join(ParticleDirName(baseDir), "tags")
}

func TagFileName(baseDir string, i int) string {
	return path.Join(TagDirName(baseDir), fmt.Sprintf("tags.%d.dat", i))
}

func SnapFileName(baseDir string, i int) string {
	return path.Join(TagDirName(baseDir), fmt.Sprintf("snaps.%d.dat", i))
}

func FlagFileName(baseDir string, i int) string {
	return path.Join(TagDirName(baseDir), fmt.Sprintf("flags.%d.dat", i))
}

func VarFileName(baseDir, varName string, snap, i int) string {
	return path.Join(SnapDirName(baseDir, snap),
		fmt.Sprintf("%s.%d.dat", varName, i))
}
