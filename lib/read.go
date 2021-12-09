package lib

type fixedHeader struct {
	NHalo, NInt, NFlaot, NVector int32
}

type Header struct {
	fixedHeader
	Cols []int32
}

func WriteTreeHeader(f *os.File, hd *Header) {
	
}

func ReadTreeHeader(f *os.File, hd *Header) {
	
}

func ReadVar(f *os.File, varName string) {
}
