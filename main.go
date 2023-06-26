package main

import (
	"RHF-go/chemical_system"
	"fmt"
	"gonum.org/v1/gonum/mat"
)

func main() {
	/*h1 := chemical_system.NewAtom("H", [3]float64{1.0, 0.0, 0.0})
	h2 := chemical_system.NewAtom("H", [3]float64{0.0, 1.0, 0.0})
	o1 := chemical_system.NewAtom("O", [3]float64{0.0, 0.0, 0.0})
	h2o := chemical_system.NewMolecule([]*chemical_system.Atom{o1, h1, h2}, &chemical_system.STO3G)*/

	h2o := chemical_system.NewMoleculeFromXYZ("./xyzfiles/Water.xyz", &chemical_system.STO3G)

	h2o.CalcS()
	fmt.Print(mat.Formatted(h2o.GetS()))
}
