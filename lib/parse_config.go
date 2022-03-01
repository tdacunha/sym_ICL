package lib

import (
	"io/ioutil"
	"fmt"
	"strings"
	"strconv"
)

// Config contains all the fields inside a particle-tracking config file.
type Config struct {
	Eps, Mp []float64
	Blocks []int
	// Input directories.
	SnapFormat, TreeDir []string
	// The directory where all the library's files are stored
	BaseDir []string
}

// ParseConfig parses a six-column particle-tracking config file. The contents
// of this file are documented in the README.
func ParseConfig(fname string) *Config {
	cfg := &Config{ }
	
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
		
		if len(cols) != 6 {
			panic(fmt.Sprintf("Line %d of %s is '%s', but you need there " +
				"to be six columns.", i+1, fname, line))
		}

		eps, err := strconv.ParseFloat(cols[0], 64)
		if err != nil {
			panic(fmt.Sprintf("Could not parse eps on line %d of " +
				"%s: %s", i+1, fname, cols[0]))
		}
		cfg.Eps = append(cfg.Eps, eps)

		mp, err := strconv.ParseFloat(cols[1], 64)
		if err != nil {
			panic(fmt.Sprintf("Could not parse mp on line %d of " +
				"%s: %s", i+1, fname, cols[1]))
		}
		cfg.Mp = append(cfg.Mp, mp)

		blocks, err := strconv.Atoi(cols[2])
		if err != nil {
			panic(fmt.Sprintf("Could not parse block number on line %d of " +
				"%s: %s", i+1, fname, cols[2]))
		}
		cfg.Blocks = append(cfg.Blocks, blocks)
		
		cfg.SnapFormat = append(cfg.SnapFormat, cols[3])
		cfg.TreeDir = append(cfg.TreeDir, cols[4])
		cfg.BaseDir = append(cfg.BaseDir, cols[5])
	}
	
	return cfg
}
