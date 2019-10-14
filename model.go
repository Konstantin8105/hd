package hd

import (
	"bytes"
	"fmt"
	"io"
	"math"
	"os"
	"runtime/debug"
	"sort"

	"github.com/Konstantin8105/errors"
	"github.com/Konstantin8105/sparse"
	"gonum.org/v1/gonum/mat"
)

// Model is structural calculation model
type Model struct {
	// Points is slice of point coordinate
	//
	//	[0] - X coordinate
	//	[1] - Y coordinate
	//
	Points [][2]float64

	// Beams is slice of point index and beam property
	Beams []BeamProp

	// Pins is slice of pins for beams in local system coordinate.
	// Len of support must be same amount of beam.
	// Or if len is zero, then all DoF(degree of freedom) is rigid.
	//
	// first index is point index
	//
	//	[0] - X on start point
	//	[1] - Y on start point
	//	[2] - M on start point
	//	[3] - X on end point
	//	[4] - Y on end point
	//	[5] - M on end point
	//
	// if `true` then free degree of freedom
	Pins [][6]bool

	// Supports is slice of fixed supports.
	// Len of support must be same amount of Points
	//
	// first index is point index
	//
	//	[0] - X
	//	[1] - Y
	//	[2] - M
	//
	Supports [][3]bool
}

// BeamProp is beam property
type BeamProp struct {
	// Start and end point index
	//
	//	[0] - start of beam
	//	[1] - end of beam
	//
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
	//	[0] - X
	//	[1] - Y
	//	[2] - M
	//
	// Unit: meter and rad.
	PointDisplacementGlobal [][3]float64

	// BeamForces is beam forces in local system coordinate.
	// Return data.
	//
	// first index is beam index
	//
	//	[0] - Fx on start point
	//	[1] - Fy on start point
	//	[2] - M  on start point
	//	[3] - Fx on end point
	//	[4] - Fy on end point
	//	[5] - M  on end point
	//
	// Unit: N and N*m
	BeamForces [][6]float64

	// Reactions is reaction loads in support points.
	// Return data.
	//
	// first index is point index
	//
	//	[0] - Fx
	//	[1] - Fy
	//	[2] - M
	//
	// Unit: N and N*m
	Reactions [][3]float64

	// TODO : replace to specific struct
	// Amount of calculated forms.
	// If Amount is zero or less zero, then no calculate.
	// LinearBuckling is linear buckling calculation
	AmountLinearBuckling uint16

	// Result of linear buckling calculation
	LinearBucklingResult []BucklingResult
}

func (lc LoadCase) String() (out string) {
	out += fmt.Sprintf("\nLoad case:\n")
	out += fmt.Sprintf("%5s %15s %15s %15s\n",
		"Point", "Fx, N", "Fy, N", "M, N*m")
	for _, ln := range lc.LoadNodes {
		out += fmt.Sprintf("%5d %15.5f %15.5f %15.5f\n",
			ln.N, ln.Forces[0], ln.Forces[1], ln.Forces[2])
	}
	if len(lc.PointDisplacementGlobal) > 0 {
		out += fmt.Sprintf("Point displacament in global system coordinate:\n")
		out += fmt.Sprintf("%5s %15s %15s\n", "Point", "DX, m", "DY, m")
		for i := 0; i < len(lc.PointDisplacementGlobal); i++ {
			out += fmt.Sprintf("%5d %15.5e %15.5e\n",
				i, lc.PointDisplacementGlobal[i][0], lc.PointDisplacementGlobal[i][1])
		}
	}
	// results
	if len(lc.BeamForces) > 0 {
		out += fmt.Sprintf("Local force in beam:\n")
		out += fmt.Sprintf("%5s %45s %45s\n",
			"", "START OF BEAM     ", "END OF BEAM       ")
		out += fmt.Sprintf("%5s %45s %45s\n",
			"",
			"------------------------------------",
			"------------------------------------")
		out += fmt.Sprintf("%5s %15s %15s %15s %15s %15s %15s\n",
			"Index", "Fx, N", "Fy, N", "M, N*m", "Fx, N", "Fy, N", "M, N*m")
		for i := 0; i < len(lc.BeamForces); i++ {
			out += fmt.Sprintf("%5d ", i)
			for j := 0; j < 6; j++ {
				out += fmt.Sprintf("%15.5e", lc.BeamForces[i][j])
				if j < 5 {
					out += " "
				}
			}
			out += fmt.Sprintf("\n")
		}
	}
	if len(lc.Reactions) > 0 {
		out += fmt.Sprintf("Reaction on support:\n")
		out += fmt.Sprintf("%5s %15s %15s %15s\n",
			"Index", "Fx, N", "Fy, N", "M, N*m")
		for i := 0; i < len(lc.Reactions); i++ {
			if lc.Reactions[i][0] == 0 &&
				lc.Reactions[i][1] == 0 &&
				lc.Reactions[i][2] == 0 {
				continue
			}
			out += fmt.Sprintf("%5d %15.5e %15.5e %15.5e\n",
				i, lc.Reactions[i][0], lc.Reactions[i][1], lc.Reactions[i][2])
		}
	}

	// output linear buckling data
	if len(lc.LinearBucklingResult) > 0 {
		out += fmt.Sprintf("\nLinear buckling result:\n")
	} else if lc.AmountLinearBuckling > 0 {
		out += fmt.Sprintf("\nLinear buckling result: is not calculated\n")
	}
	for _, lbr := range lc.LinearBucklingResult {
		out += fmt.Sprintf("Linear buckling factor: %15.5f\n", lbr.Factor)
		out += fmt.Sprintf("%5s %15s %15s %15s\n",
			"Point", "X", "Y", "M")
		for i := 0; i < len(lbr.PointDisplacementGlobal); i++ {
			out += fmt.Sprintf("%5d %15.5e %15.5e %15.5e\n",
				i,
				lbr.PointDisplacementGlobal[i][0],
				lbr.PointDisplacementGlobal[i][1],
				lbr.PointDisplacementGlobal[i][2])
		}
	}

	return
}

// reset - set nil to old results
func (lc *LoadCase) reset() {
	lc.PointDisplacementGlobal = nil
	lc.BeamForces = nil
	lc.Reactions = nil
	lc.LinearBucklingResult = nil
}

// ModalCase is modal calculation case
type ModalCase struct {
	// ModalMasses is modal masses
	ModalMasses []ModalMass

	// Result of modal calculation
	Result []ModalResult
}

// reset - set nil to old results
func (mc *ModalCase) reset() {
	mc.Result = nil
}

func (m ModalCase) String() (out string) {
	out += fmt.Sprintf("\nModal case:\n")
	out += fmt.Sprintf("%5s %15s\n",
		"Point", "Mass, N")
	for _, mn := range m.ModalMasses {
		out += fmt.Sprintf("%5d %15.5f\n", mn.N, mn.Mass)
	}
	for _, mr := range m.Result {
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
	return
}

// ModalResult is result of modal calculation
type ModalResult struct {
	// Natural frequency
	Hz float64

	// Modal displacament in global system coordinate
	//
	// first index is point index
	//
	//	[0] - X direction
	//	[1] - Y direction
	//	[2] - M direction
	//
	// Unit: Dimensionless
	ModalDisplacement [][3]float64
}

// BucklingResult is result of buckling calculation
type BucklingResult struct {
	// Buckling factor
	// Return data.
	Factor float64

	// Point displacement in global system coordinate.
	// Return data.
	//
	// first index is point index
	//
	//	[0] - X
	//	[1] - Y
	//	[2] - M
	//
	// Unit: relative
	PointDisplacementGlobal [][3]float64
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
	//	[0] - X , Unit: N. Positive direction from left to right.
	//	[1] - Y , Unit: N. Positive direction from down to top.
	//	[2] - M , Unit: N*m. Positive direction is counter-clockwise direction.
	//
	Forces [3]float64
}

// getK return stiffiner matrix, LU decomposition, ignore list
func getK(m *Model) (k *sparse.Matrix, lu sparse.LU, ignore []int, err error) {
	// assembly matrix of stiffiner
	k, ignore, err = m.assemblyK(m.getStiffBeam2d)
	if err != nil {
		return
	}

	// LU decomposition
	err = lu.Factorize(k,
		// add support
		append(ignore, m.addSupport()...)...)
	if err != nil {
		err = fmt.Errorf("LU error factorization: %v", err)
	}
	return
}

// prepare input Model
func prepare(in io.Writer, m *Model) (out io.Writer, err error) {
	// by default output in standart stdio
	out = in
	if out == nil {
		out = os.Stdout
	}

	// if pins is empty , then all rigid. So, create with all false
	if len(m.Pins) == 0 {
		m.Pins = make([][6]bool, len(m.Beams))
	}

	// check
	err = m.checkInputData()
	if err != nil {
		return
	}

	return
}

// LinearStatic run linear static analysis.
func LinearStatic(out io.Writer, m *Model, lc *LoadCase) (err error) {
	if out == nil {
		var buf bytes.Buffer
		out = &buf
		_ = buf
	}
	// remove result data
	lc.reset()

	// error handling
	name := "Linear Elastic Analysis"
	fmt.Fprintf(out, "%s\n", name)
	defer func() {
		if err != nil {
			et := errors.New(name)
			et.Add(err)
			err = et
		}
	}()

	// panic free. replace to stacktrace
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("stacktrace from panic: %s\n", debug.Stack())
		}
	}()

	// prepare input data
	{
		et := errors.New("")
		out, err = prepare(out, m)
		et.Add(err)
		et.Add(lc.checkInputData(m))
		if et.IsError() {
			return et
		}
	}

	// calculate node displacament
	dof := 3 * len(m.Points)
	d := make([]float64, dof)

	// templorary data for displacement in global system coordinate
	data := make([]float64, 6)
	Zo := mat.NewDense(6, 1, data)

	// assembly node load
	p, err := m.assemblyNodeLoad(lc)
	if err != nil {
		return fmt.Errorf("Assembly node load: %v", err)
	}

	// generate stiffiner matrix and ignore list
	k, lu, ignore, err := getK(m)
	if err != nil {
		return fmt.Errorf("Assembly node load: %v", err)
	}

	{
		// check loads on ignore free directions
		// ignore load in free direction. usually for pin connection
		var et errors.Tree
		for _, i := range ignore {
			if p[i] != 0.0 {
				et.Add(fmt.Errorf("on direction %d load is not zero : %f", i, p[i]))
			}
		}
		if et.IsError() {
			et.Name = "Warning: List loads on not valid directions"
			return et
		}
	}

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

	// TODO : split to specific function
	if lc.AmountLinearBuckling > 0 {
		// assembly matrix of stiffiner
		g, _, err := m.assemblyK(func(pos int) *mat.Dense {
			return m.getGeometricBeam2d(pos, lc)
		})
		if err != nil {
			return err
		}

		// memory initialization
		dof := 3 * len(m.Points)

		// templorary data for mass preparing
		MS := make([]float64, dof)

		dataM := make([]float64, dof*dof)
		M := mat.NewDense(dof, dof, dataM)

		if _, err = sparse.Fkeep(g, func(i, j int, x float64) bool {
			M.Set(i, j, M.At(i, j)+x)
			// kept value
			return true
		}); err != nil {
			return err
		}

		// for _, mm := range mc.ModalMasses {
		// index := mm.N*3 + mcCase.direction
		// M.Set(index, index, M.At(index, index)+mm.Mass/Gravity)
		// }

		dataH := make([]float64, dof*dof)
		h := mat.NewSymDense(dof, dataH)

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
				h.SetSym(i, col, hh[i])
			}
		}

		var e mat.EigenSym

		ok := e.Factorize(h, true)
		if !ok {
			return fmt.Errorf("Eigen factorization is not ok")
		}

		// create result report
		v := e.Values(nil)
		eVector := mat.NewDense(len(v), len(v), nil)
		e.VectorsTo(eVector)
		for i := 0; i < len(v); i++ {
			if v[i] == 0 {
				continue
			}

			// store the result
			var lbr BucklingResult
			if val := math.Abs(v[i]); val != 0.0 {
				// use only possitive value
				lbr.Factor = 1. / val
			} else {
				// ignore that value
				continue
			}
			for p := 0; p < len(m.Points); p++ {
				lbr.PointDisplacementGlobal = append(lbr.PointDisplacementGlobal, [3]float64{
					eVector.At(3*p+0, i),
					eVector.At(3*p+1, i),
					eVector.At(3*p+2, i),
				})
			}
			lc.LinearBucklingResult = append(lc.LinearBucklingResult, lbr)
		}
		// Sort by factors
		sort.SliceStable(lc.LinearBucklingResult, func(i, j int) bool {
			return lc.LinearBucklingResult[i].Factor < lc.LinearBucklingResult[j].Factor
		})

		// Cut result slice
		if len(lc.LinearBucklingResult) > int(lc.AmountLinearBuckling) {
			lc.LinearBucklingResult = lc.LinearBucklingResult[:lc.AmountLinearBuckling]
		}
	}

	return nil
}

func (m *Model) assemblyK(elementMatrix func(int) *mat.Dense) (
	k *sparse.Matrix, ignore []int, err error) {

	var T *sparse.Triplet
	if T, err = sparse.NewTriplet(); err != nil {
		return
	}

	for i := range m.Beams {
		kr := elementMatrix(i) // m.getStiffBeam2d(i)
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
						if err = sparse.Entry(T, x, y, val); err != nil {
							return
						}
					}
				}
			}
		}
	}

	// from triplet to sparse matrix
	k, err = sparse.Compress(T)
	if err != nil {
		return
	}

	// remove zero elements
	if _, err = sparse.Fkeep(k, func(i, j int, x float64) bool {
		return x != 0.0
	}); err != nil {
		return
	}

	// remove duplicate matrix
	err = sparse.Dupl(k)
	if err != nil {
		return
	}

	// zero diagonals
	r, c := k.Dims()
	if r != c {
		err = fmt.Errorf("matrix is not symmetric")
		return
	}
	bz := make([]bool, r)
	if _, err = sparse.Fkeep(k, func(i, j int, x float64) bool {
		if i == j { // diagonal
			bz[i] = true
		}
		// keep entry
		return true
	}); err != nil {
		return
	}

	for i := range bz {
		if !bz[i] {
			ignore = append(ignore, i)
		}
	}

	return k, ignore, nil
}

func (m *Model) assemblyNodeLoad(lc *LoadCase) (p []float64, err error) {
	dof := 3 * len(m.Points)
	p = make([]float64, dof)
	// node loads
	for _, ln := range lc.LoadNodes {
		for i := 0; i < 3; i++ {
			p[ln.N*3+i] += ln.Forces[i]
		}
	}
	for i := range p {
		if math.IsNaN(p[i]) || math.IsInf(p[i], 0) {
			return nil, fmt.Errorf("not valid node load %e", p[i])
		}
	}
	return p, nil
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

// Gravity is Earth gravity constant, m/sq.sec.
const Gravity float64 = 9.80665

func Modal(out io.Writer, m *Model, mc *ModalCase) (err error) {
	if out == nil {
		var buf bytes.Buffer
		out = &buf
		_ = buf
	}
	// remove result data
	mc.reset()

	// error handling
	name := "Modal Analysis"
	fmt.Fprintf(out, "%s\n", name)
	defer func() {
		if err != nil {
			et := errors.New(name)
			et.Add(err)
			err = et
		}
	}()

	// panic free. replace to stacktrace
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("stacktrace from panic: %s\n", debug.Stack())
		}
	}()

	// prepare input data
	{
		et := errors.New("")
		out, err = prepare(out, m)
		et.Add(err)
		et.Add(mc.checkInputData(m))
		if et.IsError() {
			return et
		}
	}

	// get LU decomposition of stiffiner matrix
	_, lu, _, err := getK(m)
	if err != nil {
		return fmt.Errorf("Assembly node load: %v", err)
	}

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

	// templorary data for mass preparing
	MS := make([]float64, dof)

	var e mat.EigenSym

	for _, mcCase := range mcCases {
		fmt.Fprintf(out, "%s", mcCase.name)

		dataM := make([]float64, dof*dof)
		M := mat.NewDense(dof, dof, dataM)

		for _, mm := range mc.ModalMasses {
			index := mm.N*3 + mcCase.direction
			M.Set(index, index, M.At(index, index)+mm.Mass/Gravity)
		}

		dataH := make([]float64, dof*dof)
		h := mat.NewSymDense(dof, dataH)

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
				h.SetSym(i, col, hh[i])
			}
		}

		ok := e.Factorize(h, true)
		if !ok {
			return fmt.Errorf("Eigen factorization is not ok")
		}

		// create result report
		v := e.Values(nil)
		eVector := mat.NewDense(len(v), len(v), nil)
		e.VectorsTo(eVector)
		for i := 0; i < len(v); i++ {
			if v[i] == 0 {
				continue
			}

			var mr ModalResult

			if val := math.Abs(v[i]); val != 0.0 {
				// use only possitive value
				mr.Hz = 1. / (math.Sqrt(val) * 2.0 * math.Pi)
			} else {
				// ignore that value
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
	return
}

func Run(out io.Writer, m *Model, lcs []LoadCase, mcs []ModalCase) (err error) {
	et := errors.New("Run function")

	if out != nil {
		fmt.Fprintf(out, "%s\n", *m)
	}
	for i := range lcs {
		if err = LinearStatic(out, m, &(lcs[i])); err != nil {
			et.Add(err)
			continue
		}
		if out != nil {
			fmt.Fprintf(out, "\n\n")
			fmt.Fprintf(out, "%s\n", lcs[i])
		}
	}
	if out != nil {
		fmt.Fprintf(out, "\n\n")
	}
	for i := range mcs {
		if err = Modal(out, m, &(mcs[i])); err != nil {
			et.Add(err)
			continue
		}
		if out != nil {
			fmt.Fprintf(out, "\n\n")
			fmt.Fprintf(out, "%s\n", mcs[i])
		}
	}

	if out != nil {
		fmt.Fprintf(out, "\n\n")
	}

	if et.IsError() {
		return et
	}
	return
}
