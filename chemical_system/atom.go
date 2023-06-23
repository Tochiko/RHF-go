package chemical_system

const A0 = 0.529177210903 // 1 Bohr

type Atom struct {
	data  *AtomicData
	coord [3]float64
}

func NewAtom(symbol string, coord [3]float64) *Atom {
	result := &Atom{
		data: PSE_BY_SYMBOL[symbol],
	}
	for i, c := range coord {
		result.coord[i] = c / A0
	}
	return result
}
