package main

import (
	"RHF-go/chemical_system"
	_ "RHF-go/chemical_system"
	"fmt"
)

func main() {
	h1 := chemical_system.NewAtom("H", [3]float32{1.0, 0.0, 0.0})
	h2 := chemical_system.NewAtom("H", [3]float32{0.0, 1.0, 0.0})
	o1 := chemical_system.NewAtom("O", [3]float32{0.0, 0.0, 0.0})
	h2o := chemical_system.NewMolecule([]*chemical_system.Atom{h1, h2, o1})

	fmt.Print(*h2o)
}
