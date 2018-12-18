package hd

import (
	"fmt"
	"io"
	"math"
	"os"
	"sort"

	"github.com/Konstantin8105/errors"
	"github.com/Konstantin8105/sparse"
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
	// Unit: meter and rad.
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
	// Unit: N and N*m
	BeamForces [][6]float64

	// Reactions is reaction loads in support points.
	// Return data.
	//
	// first index is point index
	//
	// [0] - Fx
	// [1] - Fy
	// [2] - M
	// Unit: N and N*m
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
	// remove result data
	for ind := 0; ind < len(m.LoadCases); ind++ {
		m.LoadCases[ind].PointDisplacementGlobal = nil
		m.LoadCases[ind].BeamForces = nil
		m.LoadCases[ind].Reactions = nil
	}
	for ind := 0; ind < len(m.ModalCases); ind++ {
		m.ModalCases[ind].Result = nil
	}

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
	if err := m.runLinearElastic(); err != nil {
		eCalc.Add(fmt.Errorf("Error in load case :%v", err))
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

func (m *Model) runLinearElastic() (err error) {
	fmt.Fprintf(m.out, "Linear Elastic Analysis\n")

	// assembly matrix of stiffiner
	k, ignore := m.assemblyK()

	// add support
	ignore2 := m.addSupport()
	ignore = append(ignore, ignore2...)

	// LU decomposition
	var lu sparse.LU
	err = lu.Factorize(k, ignore...)
	if err != nil {
		return fmt.Errorf("LU error factorization: %v", err)
	}
	// TODO : need concurency solver

	// repair stiffiner matrix
	// k = m.assemblyK()

	// calculate node displacament
	dof := 3 * len(m.Points)

	d := make([]float64, dof)

	// templorary data for displacement in global system coordinate
	data := make([]float64, 6)
	Zo := mat.NewDense(6, 1, data)

	for ilc := range m.LoadCases {
		lc := &m.LoadCases[ilc]

		fmt.Fprintf(m.out, "Calculate load case %d of %d\n", ilc, len(m.LoadCases))

		// assembly node load
		p := m.assemblyNodeLoad(lc)

		// solve by LU decomposition
		d, err = lu.Solve(p)
		if err != nil {
			return fmt.Errorf("Linear Elastic calculation error: %v", err)
		}

		// create result information
		lc.PointDisplacementGlobal = make([][3]float64, len(m.Points))
		for p := 0; p < len(m.Points); p++ {
			for i := 0; i < 3; i++ {
				lc.PointDisplacementGlobal[p][i] = d[3*p+i]
			}
		}

		lc.BeamForces = make([][6]float64, len(m.Beams))
		lc.Reactions = make([][3]float64, len(m.Points))
		for bi, b := range m.Beams {
			for i := 0; i < 3; i++ {
				for j := 0; j < 2; j++ {
					Zo.Set(j*3+i, 0, d[b.N[j]*3+i])
				}
			}
			tr := m.getCoordTransStiffBeam2d(bi)
			var z mat.Dense
			z.Mul(tr, Zo)
			kr := m.getStiffBeam2d(bi)
			s := mat.NewDense(6, 1, lc.BeamForces[bi][:])
			// calculate beam forces
			s.Mul(kr, &z)
		}

		// calculate reactions
		for pt := 0; pt < len(m.Points); pt++ {
			for i := 0; i < 3; i++ {
				if !m.Supports[pt][i] {
					// free support
					continue
				}
				// fix support
				react := -p[3*pt+i]

				row := 3*pt + i
				_, _ = sparse.Fkeep(k, func(i, j int, x float64) bool {
					if i == row {
						react += x * d[j]
					}
					return true
				})
				lc.Reactions[pt][i] = react
			}
		}
	}

	return nil
}

func (m *Model) assemblyK() (k *sparse.Matrix, ignore []int) {
	T, err := sparse.NewTriplet()
	if err != nil {
		panic(err)
	}

	for i := range m.Beams {
		kr := m.getStiffBeam2d(i)
		tr := m.getCoordTransStiffBeam2d(i)

		kr.Mul(tr.T(), kr)
		kr.Mul(kr, tr)

		for p1 := 0; p1 < 2; p1++ {
			for p2 := 0; p2 < 2; p2++ {
				for r1 := 0; r1 < 3; r1++ {
					for r2 := 0; r2 < 3; r2++ {
						x := m.Beams[i].N[p1]*3 + r1
						y := m.Beams[i].N[p2]*3 + r2
						val := kr.At(p1*3+r1, p2*3+r2)
						// relocate zero checking to sparse package
						// if val == 0.0 {
						// 	continue
						// }
						if err := sparse.Entry(T, x, y, val); err != nil {
							panic(err)
						}
					}
				}
			}
		}
	}

	// from triplet to sparse matrix
	k, err = sparse.Compress(T)
	if err != nil {
		panic(err)
	}

	// remove zero elements
	_, _ = sparse.Fkeep(k, func(i, j int, x float64) bool {
		return x != 0.0
	})

	// remove duplicate matrix
	err = sparse.Dupl(k)
	if err != nil {
		panic(err)
	}

	// singinal check
	min, max := math.MaxFloat64, 0.0
	_, err = sparse.Fkeep(k, func(i, j int, x float64) bool {
		if i == j { // diagonal
			if math.Abs(x) > max {
				max = math.Abs(x)
			}
			if math.Abs(x) < min {
				min = math.Abs(x)
			}
		}
		// keep entry
		return true
	})
	if err != nil {
		panic(err)
	}
	if min == 0 {
		panic("singular: zero entry on diagonal")
	}
	if max/min > 1e18 {
		panic(fmt.Sprintf("singular: max/min diagonal entry: %v", max/min))
	}

	// zero diagonals
	r, c := k.Dims()
	if r != c {
		panic("is not symmetric")
	}
	bz := make([]bool, r)
	_, err = sparse.Fkeep(k, func(i, j int, x float64) bool {
		if i == j { // diagonal
			bz[i] = true
		}
		// keep entry
		return true
	})

	for i := range bz {
		if !bz[i] {
			ignore = append(ignore, i)
			// panic(fmt.Errorf("singular %d : %v", i, bz))
		}
	}

	return k, ignore
}

func (m *Model) assemblyNodeLoad(lc *LoadCase) []float64 {
	dof := 3 * len(m.Points)
	p := make([]float64, dof)
	// node loads
	for _, ln := range lc.LoadNodes {
		for i := 0; i < 3; i++ {
			p[ln.N*3+i] += ln.Forces[i]
		}
	}
	return p
}

// TODO: need recode part of solving (Ax=b,eigen) for avoid that function
// For avoid matrix singular value of support is must be:
// 1) not zero
// 2) like absolute value of k
func getAverageValueOfK(k mat.Matrix) (value float64) {
	c, _ := k.Dims()
	for i := 0; i < c; i++ {
		if value < math.Abs(k.At(i, i)) {
			value = math.Abs(k.At(i, i))
		}
	}
	return
}

func (m *Model) addSupport() (ignore []int) {
	// choose value for support
	for n := range m.Supports {
		for i := 0; i < 3; i++ {
			if m.Supports[n][i] {
				ignore = append(ignore, n*3+i)
			}
		}
	}
	return
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

	// LU decomposition
	var lu sparse.LU
	{
		// assembly matrix of stiffiner
		k, ignore2 := m.assemblyK()

		// add support
		ignore := m.addSupport()
		ignore = append(ignore, ignore2...)

		// LU factorization
		err = lu.Factorize(k, ignore...)
		if err != nil {
			return fmt.Errorf("LU error factorization: %v", err)
		}
	}

	// templorary data for calc matrix H
	// hh := make([]float64, dof)

	// templorary data for mass preparing
	MS := make([]float64, dof)

	var e mat.Eigen

	for _, mcCase := range mcCases {
		// TODO: memory optimize for modal calc

		fmt.Fprintf(m.out, "%s", mcCase.name)

		dataM := make([]float64, dof*dof)
		M := mat.NewDense(dof, dof, dataM)

		for _, mm := range mc.ModalMasses {
			index := mm.N*3 + mcCase.direction
			M.Set(index, index, M.At(index, index)+mm.Mass/Gravity)
		}

		dataH := make([]float64, dof*dof)
		h := mat.NewDense(dof, dof, dataH)

		for col := 0; col < dof; col++ {
			for i := 0; i < dof; i++ {
				// initialization by 0.0
				MS[i] = 0.0
			}
			isZero := true
			for i := 0; i < dof; i++ {
				if M.At(i, col) != 0 {
					isZero = false
				}
				MS[i] = M.At(i, col)
			}

			if isZero {
				continue
			}

			// LU decomposition
			hh, err := lu.Solve(MS)
			if err != nil {
				return err
			}

			for i := 0; i < dof; i++ {
				h.Set(i, col, hh[i])
			}
		}

		ok := e.Factorize(h, true, true)
		if !ok {
			return fmt.Errorf("Eigen factorization is not ok")
		}

		// create result report
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
			// TODO: G-beam test have 2 frequency on 2 different direction. I think, that checking must be removed.
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
