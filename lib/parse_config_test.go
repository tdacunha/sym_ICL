package lib

import (
	"testing"
)

func TestParseConfig(t *testing.T) {
	cfgName := "test_dir/test_config.1.txt"
	cfg := ParseConfig(cfgName)

	matchID:= []int32{ 9, 8 }
	matchSnap := []int32{ 99, 88 }
	eps := []float64{ 1.0, 3.5 }
	mp := []float64{ 2.0, 5.5 }
	blocks := []int{ 10, 100 }
	snapFmt := []string{"a", "aa"}
	treeDir := []string{"b", "bb"}
	baseDir := []string{"c", "cc"}
	musicConfig := []string{"d", "dd"}
	
	if !Int32Eq(matchID, cfg.MatchID) {
		t.Errorf("Expected MatchID = %d, got %d",
			matchID, cfg.MatchID)
	}

	if !Int32Eq(matchSnap, cfg.MatchSnap) {
		t.Errorf("Expected MatchSnap = %d, got %d",
			matchSnap, cfg.MatchSnap)
	}
	
	if !Float64Eq(mp, cfg.Mp) {
		t.Errorf("Expected Mp = %g, got %g",
			mp, cfg.Mp)
	}

	if !Float64Eq(eps, cfg.Eps) {
		t.Errorf("Expected Eps = %g, got %g",
			eps, cfg.Eps)
	}

	if !IntEq(blocks, cfg.Blocks) {
		t.Errorf("Expected Blocks = %d, got %d",
			blocks, cfg.Blocks)
	}
	
	if !StringEq(snapFmt, cfg.SnapFormat) {
		t.Errorf("Expected SnapFormat = %s, got %s",
			snapFmt, cfg.SnapFormat)
	}

	if !StringEq(treeDir, cfg.TreeDir) {
		t.Errorf("Expected TreeDir = %s, got %s",
			treeDir, cfg.TreeDir)
	}
	
	if !StringEq(baseDir, cfg.BaseDir) {
		t.Errorf("Expected BaseDir = %s, got %s",
			baseDir, cfg.BaseDir)
	}

	if !StringEq(musicConfig, cfg.MusicConfig) {
		t.Errorf("Expected MusicConfig = %s, got %s",
			musicConfig, cfg.MusicConfig)
	}
}

func Float64Eq(x, y []float64) bool {
	if len(x) != len(y) { return false }
	for i := range x {
		if x[i] != y[i] { return false }
	}
	return true
}

func IntEq(x, y []int) bool {
	if len(x) != len(y) { return false }
	for i := range x {
		if x[i] != y[i] { return false }
	}
	return true
}

func StringEq(x, y []string) bool {
	if len(x) != len(y) { return false }
	for i := range x {
		if x[i] != y[i] { return false }
	}
	return true
}
