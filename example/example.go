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
			AmountLinearBuckling: 2,
		},
	}

	return
}

// GBeam - example of `fd` Model
//
//	            🡱🡲
//	0===========.============0
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
		{
			LoadNodes: []hd.LoadNode{
				{
					N:      1,
					Forces: [3]float64{-70000, 0, 0},
				}, {
					N:      3,
					Forces: [3]float64{42000, 0, 0},
				},
			},
		},
	}
	mc = []hd.ModalCase{
		{
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

func TrussWithBuckling() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
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
		{
			LoadNodes: []hd.LoadNode{
				{N: 3, Forces: [3]float64{+70000, 0, 0}},
			},
			AmountLinearBuckling: 2,
		},
		{
			LoadNodes: []hd.LoadNode{
				{N: 3, Forces: [3]float64{-70000, 0, 0}},
			},
			AmountLinearBuckling: 2,
		},
	}
	return
}

func ModalTruss() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
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

func ModalTrussRotate() (m hd.Model, lc []hd.LoadCase, mc []hd.ModalCase) {
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
