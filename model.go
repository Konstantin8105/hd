package hd

import (
	"fmt"
	"io"
	"math"
	"os"
	"sort"

	"github.com/Konstantin8105/errors"
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

	// Pins is slice of pins for beams in local system coordinate.
	// Len of support must be same amount of Points.
	// Or if len is zero, then all DoF(degree of freedom) is rigid.
	//
	// first index is point index
	//
	// [0] - X on start point
	// [1] - Y on start point
	// [2] - M on start point
	// [3] - X on end point
	// [4] - Y on end point
	// [5] - M on end point
	//
	// if `true` then free degree of freedom
	Pins [][6]bool

	// Supports is slice of fixed supports.
	// Len of support must be same amount of Points
	//
	// first index is point index
	//
	// [0] - X
	// [1] - Y
	// [2] - M
	Supports [][3]bool

	// LoadCases is slice of load cases
	LoadCases []LoadCase

	// ModalCases is slice of modal cases
	ModalCases []ModalCase

	// output file
	out io.Writer
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
}

// LoadCase is summary combination of loads
type LoadCase struct {
	// LoadNodes is node loads in global system coordinate
	LoadNodes []LoadNode

	// Point displacement in global system coordinate.
	// Return data.
	//
	// first index is point index
	//
	// [0] - X
	// [1] - Y
	// [2] - M
	// Unit: meter
	PointDisplacementGlobal [][3]float64

	// BeamForces is beam forces in local system coordinate.
	// Return data.
	//
	// first index is beam index
	//
	// [0] - Fx on start point
	// [1] - Fy on start point
	// [2] - M  on start point
	// [3] - Fx on end point
	// [4] - Fy on end point
	// [5] - M  on end point
	// Unit: N
	BeamForces [][6]float64

	// Reactions is reaction loads in support points.
	// Return data.
	//
	// first index is point index
	//
	// [0] - Fx
	// [1] - Fy
	// [2] - M
	// Unit: N
	Reactions [][3]float64
}

// ModalCase is modal calculation case
type ModalCase struct {
	// ModalMasses is modal masses
	ModalMasses []ModalMass

	// Result of modal calculation
	Result []ModalResult
}

type ModalResult struct {
	// Natural frequency
	Hz float64

	// Modal displacament in global system coordinate
	//
	// first index is point index
	//
	// [0] - X direction
	// [1] - Y direction
	// [2] - M direction
	// Unit: Dimensionless
	ModalDisplacement [][3]float64
}

// ModalMass is mass of point
type ModalMass struct {
	// N is point index
	N int

	// Mass of point
	// Unit: N
	Mass float64
}

// LoadNode is node load on specific point in global system coordinate
type LoadNode struct {
	// N is point index
	N int

	// Forces is node loads on each direction
	//
	// [0] - X , Unit: N
	//
	// [1] - Y , Unit: N
	//
	// [2] - M , Unit: N*m
	Forces [3]float64
}

// Run is run calculation of model
func (m *Model) Run(out io.Writer) (err error) {

	// by default output in standart stdio
	if out == nil {
		out = os.Stdout
	}
	m.out = out

	// if pins is empty , then all rigid. So, create with all false
	if len(m.Pins) == 0 {
		m.Pins = make([][6]bool, len(m.Beams))
	}

	// fix pin bug for truss elements
	for beam := 0; beam < len(m.Pins); beam++ {
		if m.Pins[beam][2] && m.Pins[beam][5] {
			// add free by Y direction
			m.Pins[beam][1] = true
			m.Pins[beam][4] = true
		}
	}

	// check
	err = m.checkInputData()
	if err != nil {
		return
	}

	var eCalc errors.Tree
	eCalc.Name = "Calculation errors"

	// calculation by load cases
	for i := range m.LoadCases {
		fmt.Fprintf(m.out, "Calculate load case %d of %d\n",
			i,
			len(m.LoadCases))
		if err := m.runLinearElastic(&m.LoadCases[i]); err != nil {
			eCalc.Add(fmt.Errorf("Error in load case %d: %v", i, err))
		}
	}

	// calculation by modal cases
	for i := range m.ModalCases {
		fmt.Fprintf(m.out, "Calculate modal case %d of %d\n",
			i,
			len(m.ModalCases))
		if err := m.runModal(&m.ModalCases[i]); err != nil {
			eCalc.Add(fmt.Errorf("Error in modal case %d: %v", i, err))
		}
	}

	if eCalc.IsError() {
		fmt.Fprintf(m.out, eCalc.Error())
		return eCalc
	}

	return nil
}

func (m *Model) runLinearElastic(lc *LoadCase) (err error) {
	fmt.Fprintf(m.out, "Linear Elastic Analysis\n")

	// assembly matrix of stiffiner
	k := m.assemblyK()

	// assembly node load
	p := m.assemblyNodeLoad(lc)

	// add support
	m.addSupport(k)

	// calculate node displacament
	dof := 3 * len(m.Points)
	dataDisp := make([]float64, dof)
	d := mat.NewDense(dof, 1, dataDisp)
	err = d.Solve(k, p)
	if err != nil {
		return fmt.Errorf("Linear Elastic calculation error: %v", err)
	}

	// create result information
	lc.PointDisplacementGlobal = make([][3]float64, len(m.Points))
	for p := 0; p < len(m.Points); p++ {
		for i := 0; i < 3; i++ {
			lc.PointDisplacementGlobal[p][i] = d.At(3*p+i, 0)
		}
	}

	// TODO: add model with many load cases
	// TODO: add benchmark for many load cases
	// TODO: try LU decomposition for optimize calc time

	data := make([]float64, 6)
	dataS := make([]float64, 6)
	Zo := mat.NewDense(6, 1, data)
	lc.BeamForces = make([][6]float64, len(m.Beams))
	lc.Reactions = make([][3]float64, len(m.Points))
	for bi, b := range m.Beams {
		for i := 0; i < 3; i++ {
			for j := 0; j < 2; j++ {
				Zo.Set(j*3+i, 0, d.At(b.N[j]*3+i, 0))
			}
		}
		tr := m.getCoordTransStiffBeam2d(bi)
		var z mat.Dense
		z.Mul(tr, Zo)
		kr := m.getStiffBeam2d(bi)
		s := mat.NewDense(6, 1, dataS)
		s.Mul(kr, &z)
		// calculate beam forces
		for i := 0; i < 6; i++ {
			lc.BeamForces[bi][i] = s.At(i, 0)
		}
	}
	// calculate reactions
	k = m.assemblyK()
	for pt := 0; pt < len(m.Points); pt++ {
		for i := 0; i < 3; i++ {
			if !m.Supports[pt][i] {
				// free support
				continue
			}
			// fix support
			react := -p.At(3*pt+i, 0)
			for j := 0; j < dof; j++ {
				react += k.At(3*pt+i, j) * d.At(j, 0)
			}
			lc.Reactions[pt][i] = react
		}
	}

	return nil
}

func (m *Model) assemblyK() *mat.Dense {
	dof := 3 * len(m.Points)
	data := make([]float64, dof*dof)
	k := mat.NewDense(dof, dof, data)

	dataTrT := make([]float64, dof*dof)
	trt := mat.NewDense(dof, dof, dataTrT)

	for i := range m.Beams {
		kr := m.getStiffBeam2d(i)
		tr := m.getCoordTransStiffBeam2d(i)

		trt.Clone(tr)
		trt.T()

		kr.Mul(trt, kr)
		kr.Mul(kr, tr)
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

	// It is happen if we have pin nodes
	// add 1 for diagonal with zero
	for i := 0; i < dof; i++ {
		isZero := true
		for j := 0; j < dof; j++ {
			if k.At(i, j) != 0.0 {
				isZero = false
				break
			}
		}
		if !isZero {
			continue
		}
		k.Set(i, i, 1)
	}

	return k
}

func (m *Model) assemblyNodeLoad(lc *LoadCase) (p *mat.Dense) {
	dof := 3 * len(m.Points)
	data := make([]float64, dof)
	p = mat.NewDense(dof, 1, data)
	// node loads
	for _, ln := range lc.LoadNodes {
		for i := 0; i < 3; i++ {
			p.Set(ln.N*3+i, 0, p.At(ln.N*3+i, 0)+ln.Forces[i])
		}
	}
	return
}

func (m *Model) addSupport(k *mat.Dense) {
	dof := 3 * len(m.Points)
	for n, s := range m.Supports {
		for i := 0; i < len(s); i++ {
			if s[i] {
				for j := 0; j < dof; j++ {
					k.Set(j, n*3+i, 0)
					k.Set(n*3+i, j, 0)
				}
				k.Set(n*3+i, n*3+i, 1)
			}
		}
	}
}

// Earth gravity, m/sq.sec.
const Gravity float64 = 9.80665

func (m *Model) runModal(mc *ModalCase) (err error) {
	mcCases := []struct {
		name      string
		direction int
	}{
		{
			name:      "Modal Analysis by X\n",
			direction: 0,
		}, {
			name:      "Modal Analysis by Y\n",
			direction: 1,
		},
	}

	// memory initialization
	dof := 3 * len(m.Points)

	for _, mcCase := range mcCases {
		// TODO: memory optimize for modal calc
		dataM := make([]float64, dof*dof)
		dataH := make([]float64, dof*dof)

		fmt.Fprintf(m.out, "%s", mcCase.name)

		// assembly matrix of stiffiner
		k := m.assemblyK()

		// add support
		m.addSupport(k)

		M := mat.NewDense(dof, dof, dataM)

		for _, mm := range mc.ModalMasses {
			index := mm.N*3 + mcCase.direction
			M.Set(index, index, M.At(index, index)+mm.Mass/Gravity)
		}

		h := mat.NewDense(dof, dof, dataH)

		for col := 0; col < dof; col++ {
			dataMS := make([]float64, dof)
			MS := mat.NewDense(dof, 1, dataMS)
			isZero := true
			for i := 0; i < dof; i++ {
				if M.At(i, col) != 0 {
					isZero = false
				}
				MS.Set(i, 0, M.At(i, col))
			}

			if isZero {
				continue
			}

			datahh := make([]float64, dof)
			hh := mat.NewDense(dof, 1, datahh)
			err = hh.Solve(k, MS)
			if err != nil {
				return err
			}

			for i := 0; i < dof; i++ {
				h.Set(i, col, hh.At(i, 0))
			}
		}

		var e mat.Eigen
		ok := e.Factorize(h, true, true)
		if !ok {
			return fmt.Errorf("Eigen factorization is not ok")
		}

		// add results
		v := e.Values(nil)
		eVector := e.Vectors()
		for i := 0; i < len(v); i++ {
			if math.Abs(imag(v[i])) > 0 || real(v[i]) == 0 {
				continue
			}
			if real(v[i]) < 0 {
				// TODO: change sign. Check is it important
				v[i] = complex(math.Abs(real(v[i])), 0)
			}

			var mr ModalResult
			mr.Hz = 1. / (math.Sqrt(real(v[i])) * 2.0 * math.Pi)
			var isFound bool
			for j := range mc.Result {
				if math.Abs((mc.Result[j].Hz-mr.Hz)/mr.Hz) < 1e-10 {
					isFound = true
				}
			}
			if isFound {
				continue
			}

			mr.ModalDisplacement = make([][3]float64, len(m.Points))
			for p := 0; p < len(m.Points); p++ {
				mr.ModalDisplacement[p][0] = eVector.At(3*p+0, i)
				mr.ModalDisplacement[p][1] = eVector.At(3*p+1, i)
				mr.ModalDisplacement[p][2] = eVector.At(3*p+2, i)
			}
			mc.Result = append(mc.Result, mr)
		}
	}

	// Sort by frequency
	sort.SliceStable(mc.Result, func(i, j int) bool {
		return mc.Result[i].Hz < mc.Result[j].Hz
	})

	return
}

// TODO: add picture draw. See https://github.com/fogleman/gg

func (m Model) String() (out string) {
	// TODO: use github.com/olekukonko/tablewriter for create table
	out += "\n"
	// points and supports
	out += fmt.Sprintf("Point coordinates:\n")
	out += fmt.Sprintf("%37s %17s\n", "", "   SUPPORT DIRECTION")
	out += fmt.Sprintf("%5s %15s %15s ", "Index", "X, m", "Y, m")
	out += fmt.Sprintf("%5s %5s %5s (0 - free, 1 - fixed)\n", "SX", "SY", "SM")
	for i := 0; i < len(m.Points); i++ {
		out += fmt.Sprintf("%5d %15.5f %15.5f ",
			i, m.Points[i][0], m.Points[i][1])
		for j := 0; j < 3; j++ {
			if m.Supports[i][j] {
				// fixed
				out += fmt.Sprintf("%5d", 1)
			} else {
				// free
				out += fmt.Sprintf("%5d", 0)
			}
			if j < 5 {
				out += " "
			}
		}
		out += fmt.Sprintf("\n")
	}
	// beams
	out += fmt.Sprintf("Beam property:\n")
	out += fmt.Sprintf("%5s %15s %15s ", "Index", "Start point", "End point")
	out += fmt.Sprintf("%15s %15s %15s\n",
		"Area,sq.m", "Moment inertia,m4", "Elasticity,Pa")
	for i := 0; i < len(m.Beams); i++ {
		out += fmt.Sprintf("%5d %15d %15d ",
			i, m.Beams[i].N[0], m.Beams[i].N[1])
		out += fmt.Sprintf("%15.5e %15.5e %15.5e\n",
			m.Beams[i].A, m.Beams[i].J, m.Beams[i].E)
	}
	// pins
	var pinHeader bool
	for pin := 0; pin < len(m.Pins); pin++ {
		var isPin bool
		for i := 0; i < 6; i++ {
			if m.Pins[pin][i] {
				isPin = true
			}
		}
		if !isPin {
			continue
		}
		if !pinHeader {
			out += fmt.Sprintf("Pins of beam in local system coordinate:\n")
			out += fmt.Sprintf("%5s %23s %23s\n",
				"", "START OF BEAM     ", "END OF BEAM       ")
			out += fmt.Sprintf("%5s %23s %23s\n",
				"", "------------------", "------------------")
			out += fmt.Sprintf("%5s %7s %7s %7s %7s %7s %7s\n",
				"Index", "X", "Y", "M", "X", "Y", "M")
			pinHeader = true
		}
		out += fmt.Sprintf("%5d ", pin)
		for i := 0; i < 6; i++ {
			out += fmt.Sprintf("%7v",
				m.Pins[pin][i])
			if i < 5 {
				out += " "
			}
		}
		out += fmt.Sprintf("\n")
	}
	if !pinHeader {
		out += fmt.Sprintf("All beams haven`t pins\n")
	}
	// loads
	for lc := 0; lc < len(m.LoadCases); lc++ {
		out += fmt.Sprintf("\nLoad case #%3d\n", lc)
		out += fmt.Sprintf("%5s %15s %15s %15s\n",
			"Point", "Fx, N", "Fy, N", "M, N*m")
		for _, ln := range m.LoadCases[lc].LoadNodes {
			out += fmt.Sprintf("%5d %15.5f %15.5f %15.5f\n",
				ln.N, ln.Forces[0], ln.Forces[1], ln.Forces[2])
		}
		l := m.LoadCases[lc]
		if len(l.PointDisplacementGlobal) > 0 {
			out += fmt.Sprintf("Point displacament in global system coordinate:\n")
			out += fmt.Sprintf("%5s %15s %15s\n", "Point", "DX, m", "DY, m")
			for i := 0; i < len(l.PointDisplacementGlobal); i++ {
				out += fmt.Sprintf("%5d %15.5e %15.5e\n",
					i, l.PointDisplacementGlobal[i][0], l.PointDisplacementGlobal[i][1])
			}
		}
		// results
		if len(l.BeamForces) > 0 {
			out += fmt.Sprintf("Local force in beam:\n")
			out += fmt.Sprintf("%5s %45s %45s\n",
				"", "START OF BEAM     ", "END OF BEAM       ")
			out += fmt.Sprintf("%5s %45s %45s\n",
				"",
				"------------------------------------",
				"------------------------------------")
			out += fmt.Sprintf("%5s %15s %15s %15s %15s %15s %15s\n",
				"Index", "Fx, N", "Fy, N", "M, N*m", "Fx, N", "Fy, N", "M, N*m")
			for i := 0; i < len(l.BeamForces); i++ {
				out += fmt.Sprintf("%5d ", i)
				for j := 0; j < 6; j++ {
					out += fmt.Sprintf("%15.5e", l.BeamForces[i][j])
					if j < 5 {
						out += " "
					}
				}
				out += fmt.Sprintf("\n")
			}
		}
		if len(l.Reactions) > 0 {
			out += fmt.Sprintf("Reaction on support:\n")
			out += fmt.Sprintf("%5s %15s %15s %15s\n",
				"Index", "Fx, N", "Fy, N", "M, N*m")
			for i := 0; i < len(l.Reactions); i++ {
				if l.Reactions[i][0] == 0 &&
					l.Reactions[i][1] == 0 &&
					l.Reactions[i][2] == 0 {
					continue
				}
				out += fmt.Sprintf("%5d %15.5e %15.5e %15.5e\n",
					i, l.Reactions[i][0], l.Reactions[i][1], l.Reactions[i][2])
			}
		}
	}
	// modal cases
	for mc := 0; mc < len(m.ModalCases); mc++ {
		out += fmt.Sprintf("\nModal case #%3d\n", mc)
		out += fmt.Sprintf("%5s %15s\n",
			"Point", "Mass, N")
		for _, mn := range m.ModalCases[mc].ModalMasses {
			out += fmt.Sprintf("%5d %15.5f\n", mn.N, mn.Mass)
		}
		for _, mr := range m.ModalCases[mc].Result {
			out += fmt.Sprintf("Natural frequency : %15.5f Hz\n", mr.Hz)
			out += fmt.Sprintf("%5s %15s %15s %15s\n",
				"Point", "X", "Y", "M")
			for i := 0; i < len(mr.ModalDisplacement); i++ {
				out += fmt.Sprintf("%5d %15.5e %15.5e %15.5e\n",
					i,
					mr.ModalDisplacement[i][0],
					mr.ModalDisplacement[i][1],
					mr.ModalDisplacement[i][2])
			}
		}
	}

	return
}
