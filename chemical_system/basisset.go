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
	data := map[string]interface{}{}
	err = json.Unmarshal(raw, &data)
	if err != nil {
		panic(err)
	}
	shift := [3]float64{0., 0., 0.}
	for _, element := range elements {
		eData := data["elements"].(map[string]interface{})
		aData := eData[strconv.Itoa(int(element.atnum))].(map[string]interface{})
		//bData := make([]*ElectronShell, 1)
		bData := aData["electron_shells"].([]interface{})
		for _, bdi := range bData {
			bd := bdi.(map[string]interface{})
			rawExps := bd["exponents"].([]interface{})
			exps := make([]float64, 3)
			for i := 0; i < 3; i++ {
				exps[i], _ = strconv.ParseFloat(rawExps[i].(string), 32)
			}
			rawCoefsList := bd["coefficients"].([]interface{})
			for i, angmom := range bd["angular_momentum"].([]interface{}) {
				coefs := make([]float64, 3)
				rawCoefs := rawCoefsList[i].([]interface{})
				for j := 0; j < 3; j++ {
					coefs[j], _ = strconv.ParseFloat(rawCoefs[j].(string), 32)
				}

				for _, ikm := range cartesianPower[int8(angmom.(float64))] {
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
	function_type    interface{}
	region           interface{}
	angular_momentum []interface{}
	exponents        []interface{}
	coefficients     []interface{}
}
