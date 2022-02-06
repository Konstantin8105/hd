package verif

import (
	"fmt"
	"math"

	"github.com/Konstantin8105/hd"
	"github.com/Konstantin8105/pow"
)

func compare(actual, expect float64) string {
	var tol float64
	if expect == actual {
		tol = 0
	} else if expect != 0.0 {
		tol = (actual - expect) / expect
	} else {
		tol = (actual - expect) / actual
	}
	tol *= 100.0
	return fmt.Sprintf("%9.4e %9.4e : %7.2f%%\n", actual, expect, tol)
}

func compareD(actual, expect [3]float64) (str string) {
	for i := range actual {
		str += compare(actual[i], expect[i])
	}
	return
}

func MSA21() (model hd.Model, lc hd.LoadCase, name string, isOk func(lc *hd.LoadCase) string) {
	name = `
Book:
William McGuire, Richard H.Gallagher, Ronald D.Ziemian
Matrix Structural Analysis

EXAMPLE 2.1
Page 21
`
	model = hd.Model{
		Points: [][2]float64{
			{0, 0},
			{6, 4},
			{9, 0},
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: 6000e-6, J: 1, E: 200000e6},
			{N: [2]int{1, 2}, A: 8000e-6, J: 1, E: 200000e6},
		},
		Pins: [][6]bool{
			{false, true, true, false, true, true},
			{false, true, true, false, true, true},
		},
		Supports: [][3]bool{
			{true, true, false},
			{false, false, false},
			{true, true, false},
		},
	}
	P := 500e3
	lc = hd.LoadCase{
		LoadNodes: []hd.LoadNode{
			{N: 1, Forces: [3]float64{P, 0, 0}},
		},
	}

	isOk = func(lc *hd.LoadCase) (str string) {
		var (
			expectX  = 2.41 / 1000.0
			expectY  = 0.72 / 1000.0
			expectFa = 400.6e3
			expectFb = -277.8e3
		)
		x := lc.PointDisplacementGlobal[1][0]
		y := lc.PointDisplacementGlobal[1][1]
		str += fmt.Sprintf("def x: %s", compare(x, expectX))
		str += fmt.Sprintf("def y: %s", compare(y, expectY))
		Fa := lc.BeamForces[0][3]
		Fb := lc.BeamForces[1][3]
		str += fmt.Sprintf("fa : %s", compare(Fa, expectFa))
		str += fmt.Sprintf("fa : %s", compare(Fb, expectFb))
		return
	}
	return
}

func MSA49() (model hd.Model, lc hd.LoadCase, name string, isOk func(lc *hd.LoadCase) string) {
	name = `
Book:
William McGuire, Richard H.Gallagher, Ronald D.Ziemian
Matrix Structural Analysis

EXAMPLE 4.9
Page 80
`
	model = hd.Model{
		Points: [][2]float64{
			{0, 0},
			{8, 0},
			{13, 0},
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: 6000e-6, J: 200e6 * 1e-12, E: 200000e6},
			{N: [2]int{1, 2}, A: 4000e-6, J: 50.e6 * 1e-12, E: 200000e6},
		},
		Pins: [][6]bool{
			{false, false, false, false, false},
			{false, false, false, false, false},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, true, false},
			{false, false, false},
		},
	}
	P := 5e3
	a := 45.0 * math.Pi / 180.0
	lc = hd.LoadCase{
		LoadNodes: []hd.LoadNode{
			{N: 2, Forces: [3]float64{P * math.Sin(a), -P * math.Cos(a), 0}},
		},
	}

	isOk = func(lc *hd.LoadCase) (str string) {
		var (
			expectReactXa = -3.6e3
			expectReactYa = -3.3e3
			expectReactMa = -8.8e3
			expectReactYb = 6.85e3
			expectDXc     = -19.15e-3
			expectDRc     = -0.00530
		)
		str += fmt.Sprintf("react: %s", compare(lc.Reactions[0][0], expectReactXa))
		str += fmt.Sprintf("react: %s", compare(lc.Reactions[0][1], expectReactYa))
		str += fmt.Sprintf("react: %s", compare(lc.Reactions[0][2], expectReactMa))
		str += fmt.Sprintf("react: %s", compare(lc.Reactions[1][1], expectReactYb))
		str += fmt.Sprintf("defor: %s", compare(lc.PointDisplacementGlobal[2][1], expectDXc))
		str += fmt.Sprintf("defor: %s", compare(lc.PointDisplacementGlobal[2][2], expectDRc))
		return
	}
	return
}

func MSA413() (model hd.Model, lc hd.LoadCase, name string, isOk func(lc *hd.LoadCase) string) {
	name = `
Book:
William McGuire, Richard H.Gallagher, Ronald D.Ziemian
Matrix Structural Analysis

EXAMPLE 4.13
Page 84
`
	model = hd.Model{
		Points: [][2]float64{
			{0, 5},
			{8, 5},
			{8, 0},
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: 6000e-6, J: 200e6 * 1e-12, E: 200000e6},
			{N: [2]int{1, 2}, A: 4000e-6, J: 50.e6 * 1e-12, E: 200000e6},
		},
		Pins: [][6]bool{
			{false, false, false, false, false},
			{false, false, false, false, false},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{true, true, true},
		},
	}
	P := 100e3
	a := 45.0 * math.Pi / 180.0
	M := +50e3
	lc = hd.LoadCase{
		LoadNodes: []hd.LoadNode{
			{N: 1, Forces: [3]float64{P * math.Sin(a), -P * math.Cos(a), M}},
		},
	}

	isOk = func(lc *hd.LoadCase) string {
		var (
			expectDa = [3]float64{
				+0.4414 / 1000.0,
				-0.3998 / 1000.0,
				+0.00169,
			}
		)
		return compareD(lc.PointDisplacementGlobal[1], expectDa)
	}
	return
}

func MSA67() (model hd.Model, lc hd.LoadCase, name string, isOk func(lc *hd.LoadCase) string) {
	name = `
Book:
William McGuire, Richard H.Gallagher, Ronald D.Ziemian
Matrix Structural Analysis

EXAMPLE 6.7
Page 168
`
	// random values
	L := 3 * math.Pi
	J := 50e-6
	A, E := 1.0, 2e11
	model = hd.Model{
		Points: [][2]float64{
			{0, 0},
			{L / 3.0, 0},
			{L/3.0 + L/6.0, 0},
			{L, 0},
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: A, J: J, E: E},
			{N: [2]int{1, 2}, A: A, J: J, E: E},
			{N: [2]int{2, 3}, A: A, J: J, E: E},
		},
		Pins: [][6]bool{
			{false, false, false, false, false},
			{false, false, false, false, false},
			{false, false, false, false, false},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{false, false, false},
			{true, true, true},
		},
	}
	P := 100e3
	lc = hd.LoadCase{
		LoadNodes: []hd.LoadNode{
			{N: 1, Forces: [3]float64{0, -P, 0}},
		},
	}

	isOk = func(lc *hd.LoadCase) string {
		var (
			v = -(5.0 * P * pow.E3(L)) / (1296.0 * E * J)
		)
		return compare(lc.PointDisplacementGlobal[2][1], v)
	}
	return
}

// TODO: add MSA511
// TODO: add Problem 6.12
// TODO: add MSA74, MSA75, MSA76, MSA710

func MSA81() (model hd.Model, lc hd.LoadCase, name string, isOk func(lc *hd.LoadCase) string) {
	name = `
Book:
William McGuire, Richard H.Gallagher, Ronald D.Ziemian
Matrix Structural Analysis

EXAMPLE 8.1
Page 220
`
	L := 4.0
	// k = E*A/L
	k := 1e3
	Area := 1e-3
	E := k * L / Area
	factor := 1e3
	model = hd.Model{
		Points: [][2]float64{
			{L, 0}, // a
			{L, L}, // b
			{0, L}, // c
		},
		Beams: []hd.BeamProp{
			// ab
			{
				N: [2]int{0, 1},
				A: factor * Area,
				J: factor * 1,
				E: factor * E,
			},
			// bc
			{N: [2]int{1, 2},
				A: Area,
				J: 1,
				E: E,
			},
		},
		Pins: [][6]bool{
			{false, true, true, false, true, true},
			{false, true, true, false, true, true},
			// {false, false, false, false, false, false},
		},
		Supports: [][3]bool{
			{true, true, false},
			{false, false, false},
			{true, false, true},
		},
	}
	a := 0.1
	var (
		expectQ = math.Atan(math.Pow(a, 1.0/3.0))
		expectD = L * math.Sin(expectQ)
		expectP = k * L * math.Cos(expectQ) / (1.0 + a*1.0/math.Tan(expectQ))
	)
	P := 0.99 * expectP // k * L

	lc = hd.LoadCase{
		LoadNodes: []hd.LoadNode{
			{N: 1, Forces: [3]float64{a * P, -P, 0}},
		},
	}
	lc.LinearBuckling.Amount = 1

	lc.NonlinearNR.MaxIterations = 29000
	lc.NonlinearNR.Substep = 50

	lc.NonlinearNK.MaxIterations = 29000
	lc.NonlinearNK.Substep = 50

	// TODO: add arc-method

	isOk = func(lc *hd.LoadCase) (str string) {
		var (
			expectLinearBucklingFactor = 1.0
		)
		str += fmt.Sprintf("Linear buckl: %s",
			compare(lc.LinearBuckling.Results[0].Factor, expectLinearBucklingFactor))
		str += fmt.Sprintf("expect Q : %.4e\n", expectQ)
		str += fmt.Sprintf("expect D : %.4e\n", expectD)
		str += fmt.Sprintf("expect P : %.4e  %.4e\n", expectP, expectP/(k*L))
		// compare with linear result
		str += fmt.Sprintf("Deformation: %s",
			compare(lc.PointDisplacementGlobal[1][0], expectD))
		// 		for Q := 0.0; Q < math.Pi/2.0; Q += 0.01 {
		// 			expectP = k * L * math.Cos(Q) / (1.0 + a*1.0/math.Tan(Q))
		// 			fmt.Println(Q*180/3.1415, Q, expectP)
		// 		}
		// TODO for each steps
		// TODO for arc method in point with Yc == 0
		for _, res := range lc.NonlinearNK.Results {
			d := res.PointDisplacementGlobal[1][0]
			p := res.Reactions[0][1]
			Q := math.Asin(d / L)
			expectP := k * L * math.Cos(Q) / (1.0 + a*1.0/math.Tan(Q))
			str += fmt.Sprintf("Steps NK: %s", compare(p, expectP))
		}
		for _, res := range lc.NonlinearNR.Results {
			d := res.PointDisplacementGlobal[1][0]
			p := res.Reactions[0][1]
			Q := math.Asin(d / L)
			expectP := k * L * math.Cos(Q) / (1.0 + a*1.0/math.Tan(Q))
			str += fmt.Sprintf("Steps NR: %s", compare(p, expectP))
		}
		return
	}
	return
}

func MSA91() (model hd.Model, lc hd.LoadCase, name string, isOk func(lc *hd.LoadCase) string) {
	name = `
Book:
William McGuire, Richard H.Gallagher, Ronald D.Ziemian
Matrix Structural Analysis

EXAMPLE 9.1
Page 247
`
	model = hd.Model{
		Points: [][2]float64{
			{0.0, 4.0}, // a
			{4.0, 4.0}, // b
			{4.0, 0.0}, // c
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: 2e-6, J: 1, E: 200000e6},    // ab
			{N: [2]int{1, 2}, A: 5000e-6, J: 1, E: 200000e6}, // bc
		},
		Pins: [][6]bool{
			{false, true, true, false, true, true},
			{false, true, true, false, true, true},
		},
		Supports: [][3]bool{
			{true, true, false},
			{false, false, false},
			{true, true, false},
		},
	}
	P := 500e3
	a := 0.05
	lc = hd.LoadCase{
		LoadNodes: []hd.LoadNode{
			{N: 1, Forces: [3]float64{a * P, -P, 0}},
		},
	}
	lc.NonlinearNR.MaxIterations = 29000
	lc.NonlinearNR.Substep = 50

	lc.NonlinearNK.MaxIterations = 29000
	lc.NonlinearNK.Substep = 50.0

	// TODO: add arc-method

	isOk = func(lc *hd.LoadCase) (str string) {
		var (
			expectD = 1692.0 / 1000.0
			expectP = 340.0 * 1000.0
		)
		for i := range lc.NonlinearNK.Results {
			if i == 0 {
				continue
			}
			lastD := lc.NonlinearNK.Results[i-1].PointDisplacementGlobal[1][0]
			presD := lc.NonlinearNK.Results[i].PointDisplacementGlobal[1][0]
			if lastD <= expectD && expectD <= presD {
				lastP := lc.NonlinearNK.Results[i-1].Reactions[2][1]
				presP := lc.NonlinearNK.Results[i].Reactions[2][1]
				p := lastP + (presP-lastP)*(expectD-lastD)/(presD-lastD)
				str += fmt.Sprintf("P NK: %s", compare(p, expectP))
			}
		}

		for i := range lc.NonlinearNR.Results {
			if i == 0 {
				continue
			}
			lastD := lc.NonlinearNR.Results[i-1].PointDisplacementGlobal[1][0]
			presD := lc.NonlinearNR.Results[i].PointDisplacementGlobal[1][0]
			if lastD <= expectD && expectD <= presD {
				lastP := lc.NonlinearNR.Results[i-1].Reactions[2][1]
				presP := lc.NonlinearNR.Results[i].Reactions[2][1]
				p := lastP + (presP-lastP)*(expectD-lastD)/(presD-lastD)
				str += fmt.Sprintf("P NK: %s", compare(p, expectP))
			}
		}
		return
	}
	return
}
