package chemical_system

type Molecule struct {
	atoms       []*Atom
	elements    map[int8]*AtomicData
	basisSet    *BasisSet
	aoBasis     []*Contracted3Gaussian
	nElectrons  int32
	nVElectrons int32
}

func NewMolecule(atoms []*Atom, nameBs *string) *Molecule {
	result := &Molecule{
		atoms:       atoms,
		nElectrons:  0,
		nVElectrons: 0,
	}
	for _, atom := range atoms {
		_, ok := result.elements[atom.data.atnum]
		if !ok {
			result.elements[atom.data.atnum] = atom.data
		}

	}
	result.basisSet = NewBasisSet(nameBs, result.elements)
	result.initialize()
	return result
}

func (m *Molecule) GetCountElectrons() (int32, int32) {
	for _, atom := range m.atoms {
		m.nElectrons += int32(atom.data.atnum)
		m.nVElectrons += int32(atom.data.velectrons)
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
