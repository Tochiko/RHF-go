package basisfunctions

import "RHF-go/chemical_system"

type BasisFunction interface {
	Normalize()
	Copy() *BasisFunction
	SetAtom(at *chemical_system.Atom)
}
