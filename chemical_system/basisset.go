package chemical_system

import (
	"encoding/json"
	"os"
)

var STO3G = "sto-3g"
var VSTO3G = "vsto-3g"

type BasisSet struct {
	name  *string
	basis map[string][]*Contracted3Gaussian
}

func NewBasisSet(name *string) *BasisSet {
	result := &BasisSet{name: name}
	result.initialize()
	return result
}

// todo finish implementation
func (bs *BasisSet) initialize() {
	raw, err := os.ReadFile("./" + *bs.name + ".json")
	if err != nil {
		panic(err)
	}
	data := map[string]string{}
	err = json.Unmarshal(raw, &data)
	if err != nil {
		return
	}
	print(data)
}

func (bs *BasisSet) getBasisFunctionsFor(atom *Atom) []*Contracted3Gaussian {
	return bs.basis[atom.data.symbol]
}
