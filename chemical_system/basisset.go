package chemical_system

import (
	"encoding/json"
	"math"
	"os"
	"strconv"
)

var STO3G = "sto-3g"
var VSTO3G = "vsto-3g"
var cartesianPower = map[int8][][3]int8{
	0: {{0, 0, 0}},
	1: {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
	2: {{1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {2, 0, 0}, {0, 2, 0}, {0, 0, 2}},
}

type BasisSet struct {
	name  *string
	basis map[string][]*Contracted3Gaussian
}

func NewBasisSet(name *string, elements map[int8]*AtomicData) *BasisSet {
	result := &BasisSet{
		name:  name,
		basis: make(map[string][]*Contracted3Gaussian),
	}

	raw, err := os.ReadFile("./" + *name + ".json")
	if err != nil {
		panic(err)
	}
	rBasisSet := &RawBasisSet{}
	err = json.Unmarshal(raw, rBasisSet)
	if err != nil {
		panic(err)
	}
	shift := [3]float64{0., 0., 0.}
	for _, element := range elements {
		for _, eShell := range rBasisSet.EContext[strconv.Itoa(int(element.atnum))].ElectronShells {
			for i, angmom := range eShell.AMoms {
				for _, ikm := range cartesianPower[angmom] {
					bFunc := NewContracted3Gaussian(eShell.Coefs[i], eShell.Exps, shift, ikm, angmom)
					norm := bFunc.S(bFunc)
					normCoefs := make([]float64, 3)
					for j, coef := range bFunc.GetCoefs() {
						normCoefs[j] = coef * math.Sqrt(norm)
					}
					bFunc.SetCoefs([3]float64(normCoefs))
					result.basis[element.symbol] = append(result.basis[element.symbol], bFunc)
				}
			}
		}
	}
	return result
}

func (bs *BasisSet) getBasisFunctionsFor(atom *Atom) []*Contracted3Gaussian {
	return bs.basis[atom.data.symbol]
}

type RawBasisSet struct {
	EContext map[string]*ElementContext `json:"elements"`
}

type ElementContext struct {
	ElectronShells []*ElectronShell `json:"electron_shells"`
}

type ElectronShell struct {
	AMoms []int8       `json:"angular_momentum"`
	Exps  [3]float64   `json:"exponents"`
	Coefs [][3]float64 `json:"coefficients"`
}

type StringToFloat float64
