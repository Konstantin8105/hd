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
func ConsoleBeam() hd.Model {
	return hd.Model{
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
		LoadCases: []hd.LoadCase{
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
		},
		ModalCases: []hd.ModalCase{
			{
				ModalMasses: []hd.ModalMass{{N: 1, Mass: 10000}},
			},
		},
	}
}

// GBeam - example of `fd` Model
//
//	            🡱🡲
//	0===========.============0
//
func GBeam() hd.Model {
	return hd.Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{1.0, 0.0},
			{1.0, 1.0},
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
			{N: [2]int{1, 2}, A: 12e-4, J: 120e-6, E: 2.0e11},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{true, true, true},
		},
		LoadCases: []hd.LoadCase{
			{
				LoadNodes: []hd.LoadNode{
					{N: 1, Forces: [3]float64{0, 13, 0}},
					{N: 1, Forces: [3]float64{13, 0, 0}},
				},
			},
		},
		ModalCases: []hd.ModalCase{
			{
				ModalMasses: []hd.ModalMass{{N: 1, Mass: 10000}},
			},
		},
	}
}

// DoubleBeam - example of `fd` Model
//
//	             |
//	             |
//	             |
//	             |
//	0=========== 0
//
func DoubleBeam() hd.Model {
	return hd.Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{1.0, 0.0},
			{0.0, 1.0},
			{2.0, 1.0},
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
			{N: [2]int{2, 3}, A: 12e-4, J: 120e-6, E: 2.0e11},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{true, true, true},
			{false, false, false},
		},
		LoadCases: []hd.LoadCase{
			{
				LoadNodes: []hd.LoadNode{
					{N: 1, Forces: [3]float64{0, 2.3, 0}},
					{N: 1, Forces: [3]float64{10, 0, 0}},
					{N: 3, Forces: [3]float64{0, 2.3, 0}},
					{N: 3, Forces: [3]float64{10, 0, 0}},
				},
			},
		},
		ModalCases: []hd.ModalCase{
			{
				ModalMasses: []hd.ModalMass{
					{N: 1, Mass: 10000},
					{N: 3, Mass: 10000},
				},
			},
		},
	}
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
func Truss() hd.Model {
	return hd.Model{
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
			{false, false, true, false, false, true}, // 1
			{false, false, true, false, false, true}, // 2
			{false, false, true, false, false, true}, // 3
			{false, false, true, false, false, true}, // 4
			{false, false, true, false, false, true}, // 5
			{false, false, true, false, false, true}, // 6
			{false, false, true, false, false, true}, // 7
		},
		LoadCases: []hd.LoadCase{
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
		},
		ModalCases: []hd.ModalCase{
			{
				ModalMasses: []hd.ModalMass{
					{
						N:    1,
						Mass: 10000,
					},
				},
			},
		},
	}
}

func ModalTruss() hd.Model {
	A := math.Pi * math.Pow(0.050, 2) / 4.0
	J := math.Pi * math.Pow(0.050, 4) / 64.0
	m := hd.Model{
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
		ModalCases: []hd.ModalCase{
			{
				ModalMasses: []hd.ModalMass{
					{
						N:    1,
						Mass: 100 * hd.Gravity,
					},
				},
			},
		},
	}
	return m
}

func ModalTrussRotate() hd.Model {
	A := math.Pi * math.Pow(0.050, 2) / 4.0
	J := math.Pi * math.Pow(0.050, 4) / 64.0
	m := hd.Model{
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
		ModalCases: []hd.ModalCase{
			{
				ModalMasses: []hd.ModalMass{
					{
						N:    1,
						Mass: 100 * hd.Gravity,
					},
				},
			},
		},
	}
	return m
}

func ModalBeam() hd.Model {
	A := math.Pi * math.Pow(0.050, 2) / 4.0
	J := math.Pi * math.Pow(0.050, 4) / 64.0
	m := hd.Model{
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
		ModalCases: []hd.ModalCase{
			{
				ModalMasses: []hd.ModalMass{
					{
						N:    1,
						Mass: 100 * hd.Gravity,
					},
				},
			},
		},
	}
	return m
}

func ModalBeamRotate() hd.Model {
	A := math.Pi * math.Pow(0.050, 2) / 4.0
	J := math.Pi * math.Pow(0.050, 4) / 64.0
	m := hd.Model{
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
		ModalCases: []hd.ModalCase{
			{
				ModalMasses: []hd.ModalMass{
					{
						N:    1,
						Mass: 100 * hd.Gravity,
					},
				},
			},
		},
	}
	return m
}

func ModalBeam3mass() hd.Model {
	E := 2e11
	J := 15e-4
	m := 250.0
	l := 4.0
	A := 12e-2
	return hd.Model{
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
		ModalCases: []hd.ModalCase{
			{
				ModalMasses: []hd.ModalMass{
					{N: 1, Mass: m * hd.Gravity},
					{N: 2, Mass: m * hd.Gravity},
					{N: 3, Mass: m * hd.Gravity},
				},
			},
		},
	}
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
func BeamDc() (m1, m2 hd.Model) {
	return hd.Model{
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
			LoadCases: []hd.LoadCase{
				{
					LoadNodes: []hd.LoadNode{
						{N: 1, Forces: [3]float64{0.0, 10.0, 0.0}},
						{N: 2, Forces: [3]float64{0.0, 10.0, 0.0}},
					},
				},
			},
		}, hd.Model{
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
			LoadCases: []hd.LoadCase{
				{
					LoadNodes: []hd.LoadNode{
						{N: 1, Forces: [3]float64{0.0, 10.0, 0.0}},
						{N: 3, Forces: [3]float64{0.0, 10.0, 0.0}},
					},
				},
			},
		}
}
