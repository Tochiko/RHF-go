package main

import (
	"RHF-go/chemical_system"
	"RHF-go/util"
	"errors"
	"fmt"
	"gonum.org/v1/gonum/mat"
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
		S := mol.CalcS()
		fmt.Print(mat.Formatted(S))
		fmt.Print("\n", "-----------------------------", "\n")
		T := mol.CalcT()
		fmt.Print(mat.Formatted(T))
		fmt.Print("\n", "-----------------------------", "\n")
		VNuc := mol.CalcVNuc()
		fmt.Print(mat.Formatted(VNuc))
	}

	/*args := os.Args
	arglen := len(args)
	basisset := args[arglen-2] // basisset should every time be the prae-last argument
	timeLogger := util.NewTimeLogger()
	h1 := chemical_system.NewAtom("H", [3]float64{1.0, 0.0, 0.0})
	h2 := chemical_system.NewAtom("H", [3]float64{0.0, 1.0, 0.0})
	o1 := chemical_system.NewAtom("O", [3]float64{0.0, 0.0, 0.0})
	mol := chemical_system.NewMolecule([]*chemical_system.Atom{o1, h1, h2}, &basisset, timeLogger)
	S := mol.CalcS()
	fmt.Print(mat.Formatted(S))
	fmt.Print("\n", "-----------------------------", "\n")
	T := mol.CalcT()
	fmt.Print(mat.Formatted(T))
	fmt.Print("\n", "-----------------------------", "\n")
	VNuc := mol.CalcVNuc()
	fmt.Print(mat.Formatted(VNuc))*/

}
