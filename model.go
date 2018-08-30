package hd

import (
	"errors"
	"fmt"
	"os"

	"gonum.org/v1/gonum/mat"
)

// Model is structural calculation model
type Model struct {
	// Points is slice of point coordinate
	//
	// [0] - X coordinate
	// [1] - Y coordinate
	Points [][2]float64

	// Beams is slice of point index and beam property
	Beams []BeamProp

	// Supports is slice of fixed supports.
	// Len of support must be same amount of Points
	//
	// [0] - X
	// [1] - Y
	// [2] - M
	Supports [][3]bool

	// LoadCases is slice of load cases
	LoadCases []LoadCase

	// output file
	out *os.File
}

// BeamProp is beam property
type BeamProp struct {
	// Start and end point index
	//
	// [0] - start of beam
	// [1] - end of beam
	N [2]int

	// A cross-section area
	// Unit : sq. meter.
	A float64

	// J is moment inertia
	// Unit : meter^4
	J float64

	// E is modulus of elasticity
	// Unit : Pa
	E float64

	// G is shear modulus
	// Unit : Pa
	G float64

	// Density
	// Unit : N/m^3
	Density float64

	// Coefficient of thermal expansion
	// Unit: 1/deg.C
	Cte float64
}

// LoadCase is summary combination of loads
type LoadCase struct {
	// LoadNodes is nodal loads
	LoadNodes []LoadNode
}

// LoadNode is node load on specific point
type LoadNode struct {
	// N is point index
	N int

	// Forces is node load on each direction
	//
	// [0] - X
	// [1] - Y
	// [2] - M
	Forces [3]float64
}

// Run is run calculation of model
func (m *Model) Run(out *os.File) (err error) {

	// by default output in standart stdio
	if out == nil {
		out = os.Stdout
	}
	m.out = out

	// check
	err = m.checkInputData()
	if err != nil {
		return
	}

	// calculation by load cases
	for i, lc := range m.LoadCases {
		fmt.Fprintf(m.out, "Calculate load case %d of %d\n",
			i,
			len(m.LoadCases))
		err = m.run(&lc)
		if err != nil {
			e := fmt.Sprintf("Error in load case %d: %v", i, err)
			fmt.Fprintf(m.out, e)
			return errors.New(e)
		}
	}
	return nil
}

func (m *Model) run(lc *LoadCase) (err error) {
	fmt.Fprintf(m.out, "Linear Elastic Analysis\n")

	// assembly matrix of stiffiner
	k := m.assemblyK()

	// assembly nodal load
	p := m.assemblyNodalLoad(lc)

	// calculate nodal displacament
	z, err := solve(k, p)
	if err != nil {
		return err
	}

	// TODO
	_ = z

	return nil
}

func (m *Model) assemblyK() *mat.Dense {
	dof := 3 * len(m.Points)
	data := make([]float64, dof*dof)
	k := mat.NewDense(dof, dof, data)
	for i := range m.Beams {
		kr := m.getStiffBeam2d(i)
		for p1 := 0; p1 < 2; p1++ {
			for p2 := 0; p2 < 2; p2++ {
				for r1 := 0; r1 < 3; r1++ {
					for r2 := 0; r2 < 3; r2++ {
						x := m.Beams[i].N[p1]*3 + r1
						y := m.Beams[i].N[p2]*3 + r2
						k.Set(x, y, k.At(x, y)+kr.At(p1*3+r1, p2*3+r2))
					}
				}
			}
		}
	}

	return k
}

func (m *Model) assemblyNodalLoad(lc *LoadCase) (p []float64) {
	return
}

func solve(k *mat.Dense, p []float64) (z []float64, err error) {
	return
}
