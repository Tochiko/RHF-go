package chemical_system

import (
	"RHF-go/integrals"
	"math"
)

type Contracted3Gaussian struct {
	coefs     [3]float64
	normconst [3]float64
	exps      [3]float64
	location  [3]float64
	ijk       [3]int8
	atom      *AtomicData
}

func NewContracted3Gaussian(coefs, exps, location [3]float64, ijk [3]int8, atom *AtomicData) *Contracted3Gaussian {
	result := &Contracted3Gaussian{
		coefs:    coefs,
		exps:     exps,
		location: location,
		ijk:      ijk,
		atom:     atom,
	}
	result.normalize()
	return result
}

func (cg *Contracted3Gaussian) normalize() {
	for i, exp := range cg.exps {
		a := integrals.S_ij(cg.ijk[0], cg.ijk[0], exp, exp, cg.location[0], cg.location[0])
		b := integrals.S_ij(cg.ijk[1], cg.ijk[1], exp, exp, cg.location[1], cg.location[1])
		c := integrals.S_ij(cg.ijk[2], cg.ijk[2], exp, exp, cg.location[2], cg.location[2])
		cg.normconst[i] = 1.0 / math.Sqrt(a*b*c)
	}
}

func (cg *Contracted3Gaussian) Copy() *Contracted3Gaussian {
	return &Contracted3Gaussian{
		coefs:     cg.coefs,
		exps:      cg.exps,
		location:  cg.location,
		ijk:       cg.ijk,
		atom:      cg.atom,
		normconst: cg.normconst,
	}
}
func (cg *Contracted3Gaussian) SetLocation(location [3]float64) {
	cg.location = location
}
