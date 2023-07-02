package main

import (
	"RHF-go/chemical_system"
	"RHF-go/util"
	"errors"
	"fmt"
	"os"
)

func main() {
	args := os.Args
	arglen := len(args)
	if arglen <= 1 {
		panic(errors.New("first argument msut be the basis-set like sto-3g; second argument must be a xyz-file"))
	}
	basisset := args[arglen-2] // basisset should every time be the prae-last argument
	if _, ok := chemical_system.BASIS_SETS[basisset]; ok {
		timeLogger := util.NewTimeLogger()
		xyzFilePath := args[arglen-1] // path to xyz file should every time be the last argumetn

		mol := chemical_system.NewMoleculeFromXYZ(xyzFilePath, &basisset, timeLogger)
		_ = mol.CalcS()

		//fmt.Print(mat.Formatted(S))
		fmt.Print("\n", "-----------------------------", "\n")
		_ = mol.CalcT()
		//fmt.Print(mat.Formatted(T))
	}
}

/*h1 := chemical_system.NewAtom("H", [3]float64{1.0, 0.0, 0.0})
h2 := chemical_system.NewAtom("H", [3]float64{0.0, 1.0, 0.0})
o1 := chemical_system.NewAtom("O", [3]float64{0.0, 0.0, 0.0})
h2o := chemical_system.NewMolecule([]*chemical_system.Atom{o1, h1, h2}, &chemical_system.STO3G)*/

/*h2o := chemical_system.NewMoleculeFromXYZ("./xyzfiles/Water.xyz", &chemical_system.STO3G)

h2o.CalcS()
fmt.Print(mat.Formatted(h2o.GetS()))*/
