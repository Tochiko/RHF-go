package chemical_system

import (
	"RHF-go/util"
	"gonum.org/v1/gonum/mat"
	"os"
	"strconv"
	"strings"
)

type Molecule struct {
	logger          *util.TimeLogger
	atoms           []*Atom
	elements        map[int8]*AtomicData
	basisSet        *BasisSet
	aoBasis         []*Contracted3Gaussian
	nElectrons      int32
	nVElectrons     int32
	nBasisFunctions int
	s               *mat.Dense
}

func NewMolecule(atoms []*Atom, nameBs *string, logger *util.TimeLogger) *Molecule {
	result := &Molecule{
		atoms:       atoms,
		nElectrons:  0,
		nVElectrons: 0,
		elements:    make(map[int8]*AtomicData),
		logger:      logger,
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

func NewMoleculeFromXYZ(path string, nameBs *string, logger *util.TimeLogger) *Molecule {
	data, err := os.ReadFile(path)
	if err != nil {
		panic(err)
	}
	raw := string(data)
	lines := strings.Split(strings.TrimSuffix(raw, "\n"), "\n")
	if len(lines) == 0 {
		return nil
	}
	countAtoms, err := strconv.Atoi(lines[0])
	if err != nil {
		panic(err)
	}
	atoms := make([]*Atom, countAtoms)
	for j, line := range lines[2:] { //First two lines from xyz-file are atom-number and a blank line
		properties := strings.Fields(line)
		var coords [3]float64
		for i, coord := range properties[1:] {
			coords[i], err = strconv.ParseFloat(coord, 64)
			if err != nil {
				panic(err)
			}
		}
		atoms[j] = NewAtom(properties[0], coords)
	}
	return NewMolecule(atoms, nameBs, logger)
}

func (m *Molecule) GetCountElectrons() (int32, int32) {
	for _, atom := range m.atoms {
		m.nElectrons += int32(atom.data.atnum)
		m.nVElectrons += int32(atom.data.velectrons)
	}
	return m.nElectrons, m.nVElectrons
}

func (m *Molecule) initialize() {
	seqNum := 0
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

func (m *Molecule) CalcS() *mat.Dense {
	defer m.logger.LogTimeAfterCompletion("Molecule.CalcS")()
	m.s = mat.NewDense(m.nBasisFunctions, m.nBasisFunctions, nil)
	for i := 0; i < m.nBasisFunctions; i++ {
		for j := i; j < m.nBasisFunctions; j++ {
			if i == j {
				m.s.Set(i, j, 1.)
				continue
			}
			bfi := m.aoBasis[i]
			bfj := m.aoBasis[j]

			m.s.Set(i, j, bfi.S(bfj)) // i, j are exactly the SeqNum of the basis functions
			m.s.Set(j, i, m.s.At(i, j))
		}
	}
	return m.s
}

func (m *Molecule) GetS() *mat.Dense {
	return m.s
}
