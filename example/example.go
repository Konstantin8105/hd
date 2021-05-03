package example

import (
	"math"

	"github.com/Konstantin8105/hd"
)

// ConsoleBeam - example of `fd` Model
//
//	||                   🡱
//	||==================== 🡲
//	||
//
func ConsoleBeam() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
	m = hd.Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{1.0, 0.0},
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
		},
	}
	lc = []hd.LoadCase{
		{
			LoadNodes: []hd.LoadNode{
				{N: 1, Forces: [3]float64{0, 2.3, 0}},
				{N: 1, Forces: [3]float64{10, 0, 0}},
			},
		},
		{ // test for 2 cases with different positions
			LoadNodes: []hd.LoadNode{
				{N: 1, Forces: [3]float64{10, 0, 0}},
				{N: 1, Forces: [3]float64{0, 2.3, 0}},
			},
		},
	}
	mc = []hd.ModalCase{
		{
			ModalMasses: []hd.ModalMass{{N: 1, Mass: 10000}},
		},
	}

	return
}

// BucklingBeam - example of `fd` Model
//
//	||
//	||====================  ⟽
//	||
//
func BucklingBeam() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
	m = hd.Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{0.5, 0.0},
			{1.0, 0.0},
			{1.5, 0.0},
			{2.0, 0.0},
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
			{N: [2]int{1, 2}, A: 12e-4, J: 120e-6, E: 2.0e11},
			{N: [2]int{2, 3}, A: 12e-4, J: 120e-6, E: 2.0e11},
			{N: [2]int{3, 4}, A: 12e-4, J: 120e-6, E: 2.0e11},
		},
		Supports: [][3]bool{
			{true, true, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, true, false},
		},
	}
	lc = []hd.LoadCase{
		{
			LoadNodes: []hd.LoadNode{
				{N: 4, Forces: [3]float64{-1.0, 0, 0}},
			},
			LinearBuckling: struct {
				Amount  uint16
				Results []hd.BucklingResult
			}{
				Amount: 2,
			},
		},
	}

	return
}

// GBeam - example of `fd` Model
//
//	            🡱🡲
//	0===========.============0
//
func GBeam() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
	m = hd.Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{1.0, 0.0},
			{1.0, 1.0},
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: 2.5e-3, J: 0.5208e-6, E: 2.0e11},
			{N: [2]int{1, 2}, A: 2.5e-3, J: 0.5208e-6, E: 2.0e11},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{true, true, true},
		},
	}
	lc = []hd.LoadCase{
		{
			LoadNodes: []hd.LoadNode{
				{N: 1, Forces: [3]float64{0, 13e3, 0}},
				{N: 1, Forces: [3]float64{13e3, 0, 0}},
			},
		},
	}
	mc = []hd.ModalCase{
		{
			ModalMasses: []hd.ModalMass{{N: 1, Mass: 10000}},
		},
	}
	return
}

// Truss - example of `fd` Model
//
//	0
//	|\
//	| \
//	|  0
//	| /|\
//	|/ | \
//	0--0--0
//	-  -  -
//
func Truss() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
	m = hd.Model{
		Points: [][2]float64{
			{0.0, 0.0}, // 1
			{0.0, 12.}, // 2
			{4.0, 0.0}, // 3
			{4.0, 6.0}, // 4
			{8.0, 0.0}, // 5
		},
		Beams: []hd.BeamProp{
			{ // 1
				N: [2]int{0, 1}, A: 40e-4, J: 1, E: 2.0e11,
			}, { // 2
				N: [2]int{0, 2}, A: 64e-4, J: 1, E: 2.0e11,
			}, { // 3
				N: [2]int{0, 3}, A: 60e-4, J: 1, E: 2.0e11,
			}, { // 4
				N: [2]int{1, 3}, A: 60e-4, J: 1, E: 2.0e11,
			}, { // 5
				N: [2]int{2, 3}, A: 40e-4, J: 1, E: 2.0e11,
			}, { // 6
				N: [2]int{2, 4}, A: 64e-4, J: 1, E: 2.0e11,
			}, { // 7
				N: [2]int{3, 4}, A: 60e-4, J: 1, E: 2.0e11,
			},
		},
		Supports: [][3]bool{
			{true, true, false},   // 1
			{false, false, false}, // 2
			{false, true, false},  // 3
			{false, false, false}, // 4
			{false, true, false},  // 5
		},
		Pins: [][6]bool{
			{false, true, true, false, true, true}, // 1
			{false, true, true, false, true, true}, // 2
			{false, true, true, false, true, true}, // 3
			{false, true, true, false, true, true}, // 4
			{false, true, true, false, true, true}, // 5
			{false, true, true, false, true, true}, // 6
			{false, true, true, false, true, true}, // 7
		},
	}
	lc = []hd.LoadCase{
		hd.LoadCase{
			LoadNodes: []hd.LoadNode{
				{
					N:      1,
					Forces: [3]float64{-70000, 0, 0},
				}, {
					N:      3,
					Forces: [3]float64{42000, 0, 0},
				},
			},
			LinearBuckling: struct {
				Amount  uint16
				Results []hd.BucklingResult
			}{
				Amount: 1,
			},
		},
	}
	mc = []hd.ModalCase{
		hd.ModalCase{
			ModalMasses: []hd.ModalMass{
				{
					N:    1,
					Mass: 10000,
				},
			},
		},
	}
	return
}

func BeamWithBuckling() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
	A := math.Pi * math.Pow(0.050, 2) / 4.0
	J := math.Pi * math.Pow(0.050, 4) / 64.0
	m = hd.Model{
		Points: [][2]float64{
			{0.000, 0.0}, // 0
			{0.300, 0.0}, // 1
			{0.700, 0.0}, // 2
			{1.000, 0.0}, // 3
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: A, J: J, E: 2e11},
			{N: [2]int{1, 2}, A: A, J: J, E: 2e11},
			{N: [2]int{2, 3}, A: A, J: J, E: 2e11},
		},
		Supports: [][3]bool{
			{true, true, true},    // 0
			{false, false, false}, // 1
			{false, false, false}, // 2
			{false, true, false},  // 3
		},
		Pins: [][6]bool{
			{false, false, false, false, false, false}, // 0
			{false, false, false, false, false, false}, // 1
			{false, false, false, false, false, false}, // 2
		},
	}
	lc = []hd.LoadCase{
		hd.LoadCase{
			LoadNodes: []hd.LoadNode{
				{N: 3, Forces: [3]float64{+70000, 0, 0}},
			},
			LinearBuckling: struct {
				Amount  uint16
				Results []hd.BucklingResult
			}{
				Amount: 2,
			},
		},
		hd.LoadCase{
			LoadNodes: []hd.LoadNode{
				{N: 3, Forces: [3]float64{-70000, 0, 0}},
			},
			LinearBuckling: struct {
				Amount  uint16
				Results []hd.BucklingResult
			}{
				Amount: 2,
			},
			NonlinearNR: struct {
				MaxIterations uint64
				Substep       uint64
				Results       []*hd.LoadCase
			}{
				MaxIterations: 1000,
				Substep:       1,
			},
			NonlinearNK: struct {
				MaxIterations uint64
				Substep       uint64
				Results       []*hd.LoadCase
			}{
				MaxIterations: 1000,
				Substep:       1,
			},
		},
	}
	return
}

func ModalBeam2() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
	A := math.Pi * math.Pow(0.050, 2) / 4.0
	J := math.Pi * math.Pow(0.050, 4) / 64.0
	m = hd.Model{
		Points: [][2]float64{
			{0.000, 0.0}, // 0
			{0.400, 0.0}, // 1
			{1.000, 0.0}, // 2
		},
		Beams: []hd.BeamProp{
			{ // 0
				N: [2]int{0, 1}, A: A, J: J, E: 2e11,
			}, { // 1
				N: [2]int{1, 2}, A: A, J: J, E: 2e11,
			},
		},
		Supports: [][3]bool{
			{true, true, true},    // 0
			{false, false, false}, // 1
			{true, true, true},    // 2
		},
		Pins: [][6]bool{
			{false, false, true, false, false, false}, // 0
			{false, false, false, true, false, true},  // 1
		},
	}
	mc = []hd.ModalCase{
		{
			ModalMasses: []hd.ModalMass{
				{
					N:    1,
					Mass: 100 * hd.Gravity,
				},
			},
		},
	}
	return
}

func ModalBeamRotate2() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
	A := math.Pi * math.Pow(0.050, 2) / 4.0
	J := math.Pi * math.Pow(0.050, 4) / 64.0
	m = hd.Model{
		Points: [][2]float64{
			{0.0, 0.000}, // 0
			{0.0, 0.400}, // 1
			{0.0, 1.000}, // 2
		},
		Beams: []hd.BeamProp{
			{ // 0
				N: [2]int{0, 1}, A: A, J: J, E: 2e11,
			}, { // 1
				N: [2]int{1, 2}, A: A, J: J, E: 2e11,
			},
		},
		Supports: [][3]bool{
			{true, true, true},    // 0
			{false, false, false}, // 1
			{true, true, true},    // 2
		},
		Pins: [][6]bool{
			{false, false, true, false, false, false}, // 0
			{false, false, false, true, false, true},  // 1
		},
	}
	mc = []hd.ModalCase{
		{
			ModalMasses: []hd.ModalMass{
				{
					N:    1,
					Mass: 100 * hd.Gravity,
				},
			},
		},
	}
	return
}

func ModalBeam() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
	A := math.Pi * math.Pow(0.050, 2) / 4.0
	J := math.Pi * math.Pow(0.050, 4) / 64.0
	m = hd.Model{
		Points: [][2]float64{
			{0.000, 0.0}, // 0
			{0.400, 0.0}, // 1
			{1.000, 0.0}, // 2
		},
		Beams: []hd.BeamProp{
			{ // 0
				N: [2]int{0, 1}, A: A, J: J, E: 2e11,
			}, { // 1
				N: [2]int{1, 2}, A: A, J: J, E: 2e11,
			},
		},
		Supports: [][3]bool{
			{true, true, false},   // 0
			{false, false, false}, // 1
			{false, true, false},  // 2
		},
		Pins: [][6]bool{
			{false, false, false, false, false, false}, // 0
			{false, false, false, false, false, false}, // 1
		},
	}
	mc = []hd.ModalCase{
		{
			ModalMasses: []hd.ModalMass{
				{
					N:    1,
					Mass: 100 * hd.Gravity,
				},
			},
		},
	}
	return
}

func ModalBeamRotate() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
	A := math.Pi * math.Pow(0.050, 2) / 4.0
	J := math.Pi * math.Pow(0.050, 4) / 64.0
	m = hd.Model{
		Points: [][2]float64{
			{0.0, 0.000}, // 0
			{0.0, 0.400}, // 1
			{0.0, 1.000}, // 2
		},
		Beams: []hd.BeamProp{
			{ // 0
				N: [2]int{0, 1}, A: A, J: J, E: 2e11,
			}, { // 1
				N: [2]int{1, 2}, A: A, J: J, E: 2e11,
			},
		},
		Supports: [][3]bool{
			{true, true, false},   // 0
			{false, false, false}, // 1
			{true, false, false},  // 2
		},
		Pins: [][6]bool{
			{false, false, false, false, false, false}, // 0
			{false, false, false, false, false, false}, // 1
		},
	}
	mc = []hd.ModalCase{
		{
			ModalMasses: []hd.ModalMass{
				{
					N:    1,
					Mass: 100 * hd.Gravity,
				},
			},
		},
	}
	return
}

func ModalBeam3mass() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
	E := 2e11
	J := 15e-4
	mass := 250.0
	l := 4.0
	A := 12e-2
	m = hd.Model{
		Points: [][2]float64{
			{0.0, 0.000},     // 0
			{0.0, l / 6.0},   // 1
			{0.0, l / 2.0},   // 2
			{0.0, l - l/6.0}, // 3
			{0.0, l},         // 4
		},
		Beams: []hd.BeamProp{
			{ // 0
				N: [2]int{0, 1}, A: A, J: J, E: E,
			}, { // 1
				N: [2]int{1, 2}, A: A, J: J, E: E,
			}, { // 2
				N: [2]int{2, 3}, A: A, J: J, E: E,
			}, { // 3
				N: [2]int{3, 4}, A: A, J: J, E: E,
			},
		},
		Supports: [][3]bool{
			{true, true, false},   // 0
			{false, false, false}, // 1
			{false, false, false}, // 2
			{false, false, false}, // 3
			{true, false, false},  // 4
		},
	}
	mc = []hd.ModalCase{
		{
			ModalMasses: []hd.ModalMass{
				{N: 1, Mass: mass * hd.Gravity},
				{N: 2, Mass: mass * hd.Gravity},
				{N: 3, Mass: mass * hd.Gravity},
			},
		},
	}
	return
}

// baseBeamDc return 2 models for compare maximal displacement
//
//	model 1:
//	        |         |
//	        V         V
//	0-------0---------0------0
//
//	model 2:
//	        |         |
//	        V         V
//	0-------0----0----0------0
//
func BeamDc() (m1, m2 hd.Model, lc1, lc2 []hd.LoadCase) {
	m1 = hd.Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{1.0, 0.0},
			{2.0, 0.0},
			{3.0, 0.0},
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
			{N: [2]int{1, 2}, A: 12e-4, J: 120e-6, E: 2.0e11},
			{N: [2]int{2, 3}, A: 12e-4, J: 120e-6, E: 2.0e11},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{false, false, false},
			{true, true, true},
		},
	}
	lc1 = []hd.LoadCase{
		{
			LoadNodes: []hd.LoadNode{
				{N: 1, Forces: [3]float64{0.0, 10.0, 0.0}},
				{N: 2, Forces: [3]float64{0.0, 10.0, 0.0}},
			},
		},
	}
	m2 = hd.Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{1.0, 0.0},
			{1.5, 0.0},
			{2.0, 0.0},
			{3.0, 0.0},
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
			{N: [2]int{1, 2}, A: 12e-4, J: 120e-6, E: 2.0e11},
			{N: [2]int{2, 3}, A: 12e-4, J: 120e-6, E: 2.0e11},
			{N: [2]int{3, 4}, A: 12e-4, J: 120e-6, E: 2.0e11},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{true, true, true},
		},
	}
	lc2 = []hd.LoadCase{
		{
			LoadNodes: []hd.LoadNode{
				{N: 1, Forces: [3]float64{0.0, 10.0, 0.0}},
				{N: 3, Forces: [3]float64{0.0, 10.0, 0.0}},
			},
		},
	}
	return
}

// Modal analysis of frame
func FrameModal() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
	var (
		A    = 24.0e-4
		J    = 72.0e-8
		E    = 2.00e10
		mass = 100.0 * 10 //  hd.Gravity
	)
	m = hd.Model{
		Points: [][2]float64{
			{0.00, 0.40}, // 1
			{0.20, 0.40}, // 2
			{0.40, 0.00}, // 3
			{0.40, 0.20}, // 4
			{0.40, 0.40}, // 5
			{0.70, 0.40}, // 6
			{1.00, 0.40}, // 7
		},
		Beams: []hd.BeamProp{
			{ // 1
				N: [2]int{0, 1}, A: A, J: J, E: E,
			}, { // 2
				N: [2]int{1, 4}, A: A, J: J, E: E,
			}, { // 3
				N: [2]int{2, 3}, A: A, J: J, E: E,
			}, { // 4
				N: [2]int{3, 4}, A: A, J: J, E: E,
			}, { // 5
				N: [2]int{4, 5}, A: A, J: J, E: E,
			}, { // 6
				N: [2]int{5, 6}, A: A, J: J, E: E,
			},
		},
		Supports: [][3]bool{
			{true, true, true},    // 1
			{false, false, false}, // 2
			{true, true, true},    // 3
			{false, false, false}, // 4
			{false, false, false}, // 5
			{false, false, false}, // 6
			{true, true, true},    // 7
		},
	}
	mc = []hd.ModalCase{
		{
			ModalMasses: []hd.ModalMass{
				{N: 1, Mass: mass}, // 2
				{N: 5, Mass: mass}, // 6
			},
		},
	}
	return
}

// nonlinear example
func G(isLinear bool) (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
	d := 0.14
	E := 2.05e11
	J := math.Pow(d, 4) / 12
	A := d * d
	m = hd.Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{0.0, 1.0},
			{0.0, 2.0},
			{0.0, 3.0},
			{0.0, 4.0},
			{0.0, 5.0},
			{0.1, 5.0},
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: A, J: J, E: E},
			{N: [2]int{1, 2}, A: A, J: J, E: E},
			{N: [2]int{2, 3}, A: A, J: J, E: E},
			{N: [2]int{3, 4}, A: A, J: J, E: E},
			{N: [2]int{4, 5}, A: A, J: J, E: E},
			{N: [2]int{5, 6}, A: A, J: J, E: E},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
		},
	}
	P := 600000.0
	l := hd.LoadCase{
		LoadNodes: []hd.LoadNode{
			{N: 6, Forces: [3]float64{1, -P, 0}},
		},
	}
	if isLinear {
		l.LinearBuckling.Amount = 1
	} else {
		l.NonlinearNR.MaxIterations = 5000
		l.NonlinearNR.Substep = 10
		l.NonlinearNK.MaxIterations = 5000
		l.NonlinearNK.Substep = 10
	}
	lc = append([]hd.LoadCase{}, l)
	return
}

// linear buckling
func Gframe() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
	E := 2e11
	J := 120e-6
	A := 12e-4
	m = hd.Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{0.0, 1.0},
			{0.0, 2.0},
			{0.0, 3.0},
			{0.0, 4.0},
			{3.0, 4.0},
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: A, J: J, E: E},
			{N: [2]int{1, 2}, A: A, J: J, E: E},
			{N: [2]int{2, 3}, A: A, J: J, E: E},
			{N: [2]int{3, 4}, A: A, J: J, E: E},
			{N: [2]int{4, 5}, A: A, J: 2 * J, E: E},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{true, true, false},
		},
	}
	P := 2.04039 * E * J
	l := hd.LoadCase{
		LoadNodes: []hd.LoadNode{
			{N: 4, Forces: [3]float64{0, -P, 0}},
		},
		LinearBuckling: struct {
			Amount  uint16
			Results []hd.BucklingResult
		}{
			Amount: 1,
		},
	}
	lc = append([]hd.LoadCase{}, l)
	return
}

// book:
//	Клейн Г.К., Рекач В.Г., Розенблат Г.И. Руководство к практическим занятиям
//	по курсу строительной механики.–М.: Высш. шк.,1972.–320с
// Example 16. Page 66.
//
// book:
//	П.П. Гайджуров РАСЧЕТ СТЕРЖНЕВЫХ СИСТЕМ НА УСТОЙЧИВОСТЬ И КОЛЕБАНИЯ
//	Учебное пособие
// Example 1. Page 78.
//
// Figure:
//	     | 50P             | 50P
//	 P   V                 V
//	---> *-----------------*
//	     |                 |
//	     |                 |
//	     |                 |
//	     |                 |
//	     |                 |
//	     |                 |
//	    ===               ===
func Pframe() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
	J := 518e-8
	A := 61.2e-4
	E := 2.1e11
	m = hd.Model{
		Points: [][2]float64{
			{0.0, 0.0}, // 0   ===
			{0.0, 1.0}, // 1
			{0.0, 2.0}, // 2
			{0.0, 3.0}, // 3
			{0.0, 4.0}, // 4
			{0.0, 5.0}, // 5
			{0.0, 6.0}, // 6
			{0.0, 7.0}, // 7
			{0.0, 8.0}, // 8   ***
			{1.0, 8.0}, // 9
			{2.0, 8.0}, // 10
			{3.0, 8.0}, // 11
			{4.0, 8.0}, // 12  ***
			{4.0, 7.0}, // 13
			{4.0, 6.0}, // 14
			{4.0, 5.0}, // 15
			{4.0, 4.0}, // 16
			{4.0, 3.0}, // 17
			{4.0, 2.0}, // 18
			{4.0, 1.0}, // 19
			{4.0, 0.0}, // 20  ===
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: A, J: J, E: E},
			{N: [2]int{1, 2}, A: A, J: J, E: E},
			{N: [2]int{2, 3}, A: A, J: J, E: E},
			{N: [2]int{3, 4}, A: A, J: J, E: E},
			{N: [2]int{4, 5}, A: A, J: J, E: E},
			{N: [2]int{5, 6}, A: A, J: J, E: E},
			{N: [2]int{6, 7}, A: A, J: J, E: E},
			{N: [2]int{7, 8}, A: A, J: J, E: E},
			{N: [2]int{8, 9}, A: A, J: J, E: E},
			{N: [2]int{9, 10}, A: A, J: J, E: E},
			{N: [2]int{10, 11}, A: A, J: J, E: E},
			{N: [2]int{11, 12}, A: A, J: J, E: E},
			{N: [2]int{12, 13}, A: A, J: J, E: E},
			{N: [2]int{13, 14}, A: A, J: J, E: E},
			{N: [2]int{14, 15}, A: A, J: J, E: E},
			{N: [2]int{15, 16}, A: A, J: J, E: E},
			{N: [2]int{16, 17}, A: A, J: J, E: E},
			{N: [2]int{17, 18}, A: A, J: J, E: E},
			{N: [2]int{18, 19}, A: A, J: J, E: E},
			{N: [2]int{19, 20}, A: A, J: J, E: E},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{true, true, true},
		},
	}
	P := 1000.0
	l := hd.LoadCase{
		LoadNodes: []hd.LoadNode{
			{N: 8, Forces: [3]float64{P, 0, 0}},
			{N: 8, Forces: [3]float64{0, -50 * P, 0}},
			{N: 12, Forces: [3]float64{0, -50 * P, 0}},
		},
		LinearBuckling: struct {
			Amount  uint16
			Results []hd.BucklingResult
		}{
			Amount: 1,
		},
	}
	l.NonlinearNR.MaxIterations = 50000
	l.NonlinearNR.Substep = 5
	lc = append([]hd.LoadCase{}, l)

	// todo  --- result is not same
	// reactions:
	// V           N      M    in kN, kN*m
	// -0.511	 48.452	 3.460
	// -0.489	 51.549	 3.443
	//
	// deformation of top by X in meter:
	// 0.052

	return
}

// book :
//	И. Ф. Дьяков, С. А. Чернов, А. Н. Черный
//	МЕТОД КОНЕЧНЫХ ЭЛЕМЕНТОВ
//	В РАСЧЁТАХ СТЕРЖНЕВЫХ СИСТЕМ
// page 121
//	            | 0.8P       |P
//	            V            V
//	*-----------*------------*
//	            |            |
//	            |            |
//	            |            |
//	            |            |
//	           ===          ===
// TODO: check in 3D
// func Frame() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
// 	E := 2e11
// 	m = hd.Model{
// 		Points: [][2]float64{
// 			{0.0, 4.0}, // 0
// 			{2.0, 4.0}, // 1
// 			{4.0, 4.0}, // 2
// 			{4.0, 2.0}, // 3
// 			{4.0, 0.0}, // 4
// 			{6.0, 4.0}, // 5
// 			{8.0, 4.0}, // 6
// 			{8.0, 2.0}, // 7
// 			{8.0, 0.0}, // 8
// 		},
// 		Beams: []hd.BeamProp{
// 			{N: [2]int{0, 1}, A: 24e-4, J: 32e-8, E: E},
// 			{N: [2]int{1, 2}, A: 24e-4, J: 32e-8, E: E},
// 			{N: [2]int{2, 3}, A: 24e-4, J: 32e-8, E: E},
// 			{N: [2]int{3, 4}, A: 24e-4, J: 32e-8, E: E},
// 			{N: [2]int{2, 5}, A: 48e-4, J: 64e-8, E: E},
// 			{N: [2]int{5, 6}, A: 48e-4, J: 64e-8, E: E},
// 			{N: [2]int{6, 7}, A: 24e-4, J: 32e-8, E: E},
// 			{N: [2]int{7, 8}, A: 24e-4, J: 32e-8, E: E},
// 		},
// 		Supports: [][3]bool{
// 			{true, true, false},
// 			{false, false, false},
// 			{false, false, false},
// 			{false, false, false},
// 			{true, true, true},
// 			{false, false, false},
// 			{false, false, false},
// 			{false, false, false},
// 			{true, true, true},
// 		},
// 		Pins: [][6]bool{
// 			{false, false, false, false, false, false},
// 			{false, false, false, false, false, false},
// 			{false, false, false, false, false, false},
// 			{false, false, false, false, false, false},
// 			{false, false, false, false, false, false},
// 			{false, false, false, false, false, false},
// 			{false, false, true, false, false, false},
// 			{false, false, false, false, false, false},
// 		},
// 	}
// 	P := 1.0
// 	l := hd.LoadCase{
// 		LoadNodes: []hd.LoadNode{
// 			{N: 2, Forces: [3]float64{0, -0.8 * P, 0}},
// 			{N: 6, Forces: [3]float64{0, -1.0 * P, 0}},
// 		},
// 		LinearBuckling: struct {
// 			Amount  uint16
// 			Results []hd.BucklingResult
// 		}{
// 			Amount: 1,
// 		},
// 	}
// 	lc = append([]hd.LoadCase{}, l)
//  // TODO: critical force ~ 8860 by book
// 	return
// }

// created in according:
//	EFESTS: Educational finite element software for truss structure –
//	Part 3: Geometrically nonlinear static analysis
//	Wenjie Zuo, Ke Huang and Fei Cheng
//
func EFESTS10bar() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
	L := 360e-3
	A := 10e-6
	E := 10e10
	P := 1000.0
	J := 0.1
	m = hd.Model{
		Points: [][2]float64{
			{0.00000, 0.}, // 0
			{1.0 * L, 0.}, // 1
			{2.0 * L, 0.}, // 2
			{2.0 * L, -L}, // 3
			{1.0 * L, -L}, // 4
			{0.00000, -L}, // 5
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: A, J: J, E: E},
			{N: [2]int{4, 5}, A: A, J: J, E: E},
			{N: [2]int{1, 2}, A: A, J: J, E: E},
			{N: [2]int{3, 4}, A: A, J: J, E: E},
			{N: [2]int{0, 4}, A: A, J: J, E: E},
			{N: [2]int{1, 5}, A: A, J: J, E: E},
			{N: [2]int{1, 3}, A: A, J: J, E: E},
			{N: [2]int{2, 4}, A: A, J: J, E: E},
			{N: [2]int{1, 4}, A: A, J: J, E: E},
			{N: [2]int{2, 3}, A: A, J: J, E: E},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{true, true, true},
		},
		Pins: [][6]bool{
			{false, true, true, false, true, true},
			{false, true, true, false, true, true},
			{false, true, true, false, true, true},
			{false, true, true, false, true, true},
			{false, true, true, false, true, true},
			{false, true, true, false, true, true},
			{false, true, true, false, true, true},
			{false, true, true, false, true, true},
			{false, true, true, false, true, true},
			{false, true, true, false, true, true},
		},
	}
	l := hd.LoadCase{
		LoadNodes: []hd.LoadNode{
			{N: 3, Forces: [3]float64{0, -P, 0}},
			{N: 4, Forces: [3]float64{0, -P, 0}},
		},
		LinearBuckling: struct {
			Amount  uint16
			Results []hd.BucklingResult
		}{
			Amount: 1,
		},
	}
	l.NonlinearNR.MaxIterations = 50000
	l.NonlinearNR.Substep = 5
	lc = append([]hd.LoadCase{}, l)

	// todo  --- result is not same
	// deformation on point 2:
	// linear      X  0.848     Y -3.795    - ok
	// nonlinear   X  0.722     Y -3.757
	// nonlinear   X  0.743     Y -3.757
	// deformation on point 3:
	// linear      X -0.952     Y -3.940    - ok
	// nonlinear   X -1.035     Y -3.856
	// nonlinear   X -1.011     Y -3.858

	return
}

// ConsoleBeam - example from research K.J.BATHE and S.BOLOURCHI
//
//	||
//	||==================== Moment
//	||
//
func ConsoleBathe() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
	m = hd.Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{1.0, 0.0},
			{2.0, 0.0},
			{3.0, 0.0},
			{4.0, 0.0},
			{5.0, 0.0},
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
			{N: [2]int{1, 2}, A: 12e-4, J: 120e-6, E: 2.0e11},
			{N: [2]int{2, 3}, A: 12e-4, J: 120e-6, E: 2.0e11},
			{N: [2]int{3, 4}, A: 12e-4, J: 120e-6, E: 2.0e11},
			{N: [2]int{4, 5}, A: 12e-4, J: 120e-6, E: 2.0e11},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
			{false, false, false},
		},
	}
	l := hd.LoadCase{
		LoadNodes: []hd.LoadNode{
			{N: 5, Forces: [3]float64{0, 0, 12063713}},
		},
		LinearBuckling: struct {
			Amount  uint16
			Results []hd.BucklingResult
		}{
			Amount: 1,
		},
	}
	l.NonlinearNR.MaxIterations = 50000
	l.NonlinearNR.Substep = 5
	lc = append([]hd.LoadCase{}, l)

	mc = []hd.ModalCase{}

	return
}
