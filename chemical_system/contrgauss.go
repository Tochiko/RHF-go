package chemical_system

import (
	"RHF-go/integrals"
	"math"
)

type Contracted3Gaussian struct {
	SeqNum    int
	coefs     [3]float64
	normconst [3]float64
	exps      [3]float64
	shift     [3]float64
	ijk       [3]int8
	angMom    int8
	atom      *Atom
}

func NewContracted3Gaussian(coefs, exps, shift [3]float64, ijk [3]int8, oType int8) *Contracted3Gaussian {
	result := &Contracted3Gaussian{
		coefs:  coefs,
		exps:   exps,
		shift:  shift,
		ijk:    ijk,
		angMom: oType,
	}
	result.normalize()
	return result
}

func (cg *Contracted3Gaussian) normalize() {
	for i, exp := range cg.exps {
		a := integrals.SIJ(cg.ijk[0], cg.ijk[0], exp, exp, cg.shift[0], cg.shift[0])
		b := integrals.SIJ(cg.ijk[1], cg.ijk[1], exp, exp, cg.shift[1], cg.shift[1])
		c := integrals.SIJ(cg.ijk[2], cg.ijk[2], exp, exp, cg.shift[2], cg.shift[2])
		cg.normconst[i] = 1.0 / math.Sqrt(a*b*c)
	}
}

func (cg *Contracted3Gaussian) Copy() *Contracted3Gaussian {
	return &Contracted3Gaussian{
		coefs:     cg.coefs,
		exps:      cg.exps,
		shift:     cg.shift,
		ijk:       cg.ijk,
		atom:      cg.atom,
		normconst: cg.normconst,
		angMom:    cg.angMom,
	}
}
func (cg *Contracted3Gaussian) SetSeqNum(seqNum int) {
	cg.SeqNum = seqNum
}

func (cg *Contracted3Gaussian) SetAtom(atom *Atom) {
	cg.atom = atom
	cg.shift = atom.coord
}

func (cg *Contracted3Gaussian) SetCoefs(coefs [3]float64) {
	cg.coefs = coefs
}

func (cg *Contracted3Gaussian) GetCoefs() [3]float64 {
	return cg.coefs
}

func (cg *Contracted3Gaussian) S(other *Contracted3Gaussian) float64 {
	result := 0.
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			result += integrals.SIJ(cg.ijk[0], other.ijk[0], cg.exps[i], other.exps[j], cg.shift[0], other.shift[0]) *
				integrals.SIJ(cg.ijk[1], other.ijk[1], cg.exps[i], other.exps[j], cg.shift[1], other.shift[1]) *
				integrals.SIJ(cg.ijk[2], other.ijk[2], cg.exps[i], other.exps[j], cg.shift[2], other.shift[2]) *
				cg.coefs[i] * other.coefs[j] * cg.normconst[i] * other.normconst[j]
		}
	}
	return result
}

func (cg *Contracted3Gaussian) T(other *Contracted3Gaussian) float64 {
	result := 0.
	indices := [5]int8{0, 1, 2, 0, 1}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			for k := 0; k < 3; k++ {
				ind := indices[i]
				ind_p1 := indices[i+1]
				ind_p2 := indices[i+2]
				result += integrals.TIJ(cg.ijk[ind], other.ijk[ind], cg.exps[j], other.exps[k], cg.shift[ind], other.shift[ind]) *
					integrals.SIJ(cg.ijk[ind_p1], other.ijk[ind_p1], cg.exps[j], other.exps[k], cg.shift[ind_p1], other.shift[ind_p1]) *
					integrals.SIJ(cg.ijk[ind_p2], other.ijk[ind_p2], cg.exps[j], other.exps[k], cg.shift[ind_p2], other.shift[ind_p2]) *
					cg.coefs[j] * other.coefs[k] * cg.normconst[j] * other.normconst[k]
			}
		}
	}
	return result
}

func (cg *Contracted3Gaussian) VNuc(other *Contracted3Gaussian, position [3]float64) float64 {
	result := 0.
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			result += integrals.VNucIJ(cg.ijk[0], cg.ijk[1], cg.ijk[2], other.ijk[0], other.ijk[1], other.ijk[2], cg.exps[i], other.exps[j], cg.shift, other.shift, position) *
				cg.coefs[i] * other.coefs[j] * cg.normconst[i] * other.normconst[j]
		}
	}
	return result
}

func (cg *Contracted3Gaussian) VElec(o1, o2, o3 *Contracted3Gaussian) float64 {
	result := 0.
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			for k := 0; k < 3; k++ {
				for l := 0; l < 3; l++ {
					result += cg.coefs[i] * o1.coefs[j] * o2.coefs[k] * o3.coefs[l] * cg.normconst[i] * o1.normconst[j] *
						o2.normconst[k] * o3.normconst[l] * integrals.VElecIJ(cg.ijk[0], cg.ijk[1], cg.ijk[2], o1.ijk[0],
						o1.ijk[1], o1.ijk[2], o2.ijk[0], o2.ijk[1], o2.ijk[2], o3.ijk[0], o3.ijk[1], o3.ijk[2],
						cg.exps[i], o1.exps[j], o2.exps[k], o3.exps[l], cg.shift, o1.shift, o2.shift, o3.shift)
				}
			}
		}
	}
	return result
}
