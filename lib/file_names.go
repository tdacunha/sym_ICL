package lib

import (
	"path"
	"fmt"
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
