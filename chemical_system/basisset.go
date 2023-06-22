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
	result := &BasisSet{name: name}

	raw, err := os.ReadFile("./" + *name + ".json")
	if err != nil {
		panic(err)
	}
	data := map[string]interface{}{}
	err = json.Unmarshal(raw, &data)
	if err != nil {
		panic(err)
	}
	shift := [3]float64{0., 0., 0.}
	for _, element := range elements {
		eData := data["elements"].(map[string]interface{})
		aData := eData[strconv.Itoa(int(element.atnum))].(map[string]interface{})
		bData := aData["electron_shells"].([]*ElectronShell)
		for _, bd := range bData {
			rawExps := bd.exponents
			exps := make([]float64, 3)
			for i := 0; i < 3; i++ {
				exps[i] = float64(rawExps[i])
			}
			for i, angmom := range bd.angular_momentum {
				coefs := make([]float64, 3)
				rawCoefs := bd.coefficients[i]
				for i := 0; i < 3; i++ {
					coefs[i] = float64(rawCoefs[i])
				}

				for _, ikm := range cartesianPower[angmom] {
					bFunction := NewContracted3Gaussian([3]float64(coefs), [3]float64(exps), shift, ikm, element)
					norm := bFunction.S(bFunction)
					normCoefs := make([]float64, 3)
					for i, coef := range bFunction.GetCoefs() {
						normCoefs[i] = coef * math.Sqrt(norm)
					}
					bFunction.SetCoefs([3]float64(normCoefs))
					result.basis[element.symbol] = append(result.basis[element.symbol], bFunction)
				}
			}
		}
	}
	return result
}

func (bs *BasisSet) getBasisFunctionsFor(atom *Atom) []*Contracted3Gaussian {
	return bs.basis[atom.data.symbol]
}

type ElectronShell struct {
	function_type    string
	region           string
	angular_momentum []int8
	exponents        [3]float64
	coefficients     [][3]float64
}
