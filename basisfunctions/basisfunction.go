package basisfunctions

import "RHF-go/chemical_system"

// TODO: Replace the hard-coded Contracted3Gaussians in chemical system with BasisFunction interface
type BasisFunction interface {
	Normalize()
	Copy() *BasisFunction
	SetAtom(at *chemical_system.Atom)
}
