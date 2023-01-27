package lib

import (
	"io/ioutil"
	"fmt"
	"strings"
	"strconv"
)

// Config contains all the fields inside a particle-tracking config file.
type Config struct {
	MatchID, MaxSnap []int32
	Eps, Mp []float64
	Blocks []int
	// Input directories.
	SnapFormat, TreeDir []string
	// The directory where all the library's files are stored
	BaseDir []string
	TreeStyle []string
	UMFormat []string
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
		
		if len(cols) != 10 {
			panic(fmt.Sprintf("Line %d of %s is '%s', but you need there " +
				"to be ten columns.", i+1, fname, line))
		}

		matchID, err := strconv.Atoi(cols[0])
		if err != nil {
			panic(fmt.Sprintf("Could not parse ID on line %d of " +
				"%s: %s", i+1, fname, cols[0]))
		}
		cfg.MatchID = append(cfg.MatchID, int32(matchID))


		maxSnap, err := strconv.Atoi(cols[1])
		if err != nil {
			panic(fmt.Sprintf("Could not parse snapshot on line %d of " +
				"%s: %s", i+1, fname, cols[1]))
		}
		cfg.MaxSnap = append(cfg.MaxSnap, int32(maxSnap))
		
		eps, err := strconv.ParseFloat(cols[2], 64)
		if err != nil {
			panic(fmt.Sprintf("Could not parse eps on line %d of " +
				"%s: %s", i+1, fname, cols[2]))
		}
		cfg.Eps = append(cfg.Eps, eps)

		mp, err := strconv.ParseFloat(cols[3], 64)
		if err != nil {
			panic(fmt.Sprintf("Could not parse mp on line %d of " +
				"%s: %s", i+1, fname, cols[3]))
		}
		cfg.Mp = append(cfg.Mp, mp)

		blocks, err := strconv.Atoi(cols[4])
		if err != nil {
			panic(fmt.Sprintf("Could not parse block number on line %d of " +
				"%s: %s", i+1, fname, cols[4]))
		}
		cfg.Blocks = append(cfg.Blocks, blocks)
		
		cfg.SnapFormat = append(cfg.SnapFormat, cols[5])
		cfg.TreeDir = append(cfg.TreeDir, cols[6])
		cfg.BaseDir = append(cfg.BaseDir, cols[7])
		
		treeStyle := cols[8]
		switch treeStyle {
		case "ct_rvmax", "ct_rhapsody":
		default:
			panic(fmt.Sprintf("The tree style '%s' is not recognized.",
				treeStyle))
		}
		cfg.TreeStyle = append(cfg.TreeStyle, treeStyle)
		cfg.UMFormat = append(cfg.UMFormat, cols[9])
	}
	
	return cfg
}
