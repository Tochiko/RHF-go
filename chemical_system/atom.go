package chemical_system

type Atom struct {
	data  *AtomicData
	coord [3]float32
}

func NewAtom(symbol string, coord [3]float32) *Atom {
	return &Atom{
		data:  PSE_BY_SYMBOL[symbol],
		coord: coord,
	}
}
