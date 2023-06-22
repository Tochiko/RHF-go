package chemical_system

type Molecule struct {
	atoms       []*Atom
	basisSet    *BasisSet
	aoBasis     []*Contracted3Gaussian
	nElectrons  int8
	nVElectrons int8
}

func NewMolecule(atoms []*Atom, nameBs *string) *Molecule {
	result := &Molecule{
		atoms:       atoms,
		nElectrons:  0,
		nVElectrons: 0,
	}
	result.basisSet = NewBasisSet(nameBs)
	result.initialize()
	return result
}

func (m *Molecule) GetCountElectrons() (int8, int8) {
	for _, atom := range m.atoms {
		m.nElectrons += atom.data.atnum
		m.nVElectrons += atom.data.velectrons
	}
	return m.nElectrons, m.nVElectrons
}

func (m *Molecule) initialize() {
	for _, atom := range m.atoms {
		bFunctions := m.basisSet.getBasisFunctionsFor(atom)
		for _, bf := range bFunctions {
			cp := bf.Copy()
			cp.SetLocation(atom.coord)
			m.aoBasis = append(m.aoBasis, cp)
		}
	}
}
