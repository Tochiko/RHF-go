package chemical_system

import (
	"gonum.org/v1/gonum/mat"
)

type Molecule struct {
	atoms           []*Atom
	elements        map[int8]*AtomicData
	basisSet        *BasisSet
	aoBasis         []*Contracted3Gaussian
	nElectrons      int32
	nVElectrons     int32
	nBasisFunctions int
	S               *mat.Dense
}

func NewMolecule(atoms []*Atom, nameBs *string) *Molecule {
	result := &Molecule{
		atoms:       atoms,
		nElectrons:  0,
		nVElectrons: 0,
		elements:    make(map[int8]*AtomicData),
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
	seqNum := 1
	for _, atom := range m.atoms {
		bFunctions := m.basisSet.getBasisFunctionsFor(atom)
		for _, bf := range bFunctions {
			cp := bf.Copy()
			cp.SetAtom(atom)
			cp.SetSeqNum(seqNum)
			m.aoBasis = append(m.aoBasis, cp)
			seqNum++
		}
	}
	m.nBasisFunctions = len(m.aoBasis)
}

func (m *Molecule) GetS() *mat.Dense {
	m.S = mat.NewDense(m.nBasisFunctions, m.nBasisFunctions, nil)
	for i := 0; i < m.nBasisFunctions; i++ {
		for j := i; j < m.nBasisFunctions; j++ {
			if i == j {
				m.S.Set(i, j, 1.)
				continue
			}
			bfi := m.aoBasis[i]
			bfj := m.aoBasis[j]
			if bfi.atom == bfj.atom && bfi.angMom != bfj.angMom {
				m.S.Set(i, j, 0.)
				m.S.Set(j, i, 0.)
				continue
			}

			m.S.Set(i, j, bfi.S(bfj))
			m.S.Set(j, i, m.S.At(i, j))
		}
	}
	return m.S
}
