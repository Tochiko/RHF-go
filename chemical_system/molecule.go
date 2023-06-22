package chemical_system

type Molecule struct {
	atoms []*Atom
}

func NewMolecule(atoms []*Atom) *Molecule {
	return &Molecule{
		atoms: atoms,
	}
}
