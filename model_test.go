package hd

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"testing"

	"github.com/Konstantin8105/cs"
	"github.com/pmezard/go-difflib/difflib"
)

func baseBeam() Model {
	return Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{1.0, 0.0},
		},
		Beams: []BeamProp{
			{
				N: [2]int{0, 1},
				A: 12e-4,
				J: 120e-6,
				E: 2.0e11,
			},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
		},
		LoadCases: []LoadCase{
			{
				LoadNodes: []LoadNode{
					{N: 1, Forces: [3]float64{0, 2.3, 0}},
					{N: 1, Forces: [3]float64{10, 0, 0}},
				},
			},
			{ // test for 2 cases with different positions
				LoadNodes: []LoadNode{
					{N: 1, Forces: [3]float64{10, 0, 0}},
					{N: 1, Forces: [3]float64{0, 2.3, 0}},
				},
			},
		},
		ModalCases: []ModalCase{
			{
				ModalMasses: []ModalMass{{N: 1, Mass: 10000}},
			},
		},
	}
}

func baseGBeam() Model {
	return Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{1.0, 0.0},
			{1.0, 1.0},
		},
		Beams: []BeamProp{
			{
				N: [2]int{0, 1},
				A: 12e-4,
				J: 120e-6,
				E: 2.0e11,
			},
			{
				N: [2]int{1, 2},
				A: 12e-4,
				J: 120e-6,
				E: 2.0e11,
			},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{true, true, true},
		},
		LoadCases: []LoadCase{
			{
				LoadNodes: []LoadNode{
					{N: 1, Forces: [3]float64{0, 13, 0}},
					{N: 1, Forces: [3]float64{13, 0, 0}},
				},
			},
		},
		ModalCases: []ModalCase{
			{
				ModalMasses: []ModalMass{{N: 1, Mass: 10000}},
			},
		},
	}
}

func baseDoubleBeam() Model {
	return Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{1.0, 0.0},
			{0.0, 1.0},
			{2.0, 1.0},
		},
		Beams: []BeamProp{
			{
				N: [2]int{0, 1},
				A: 12e-4,
				J: 120e-6,
				E: 2.0e11,
			},
			{
				N: [2]int{2, 3},
				A: 12e-4,
				J: 120e-6,
				E: 2.0e11,
			},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{true, true, true},
			{false, false, false},
		},
		LoadCases: []LoadCase{
			{
				LoadNodes: []LoadNode{
					{N: 1, Forces: [3]float64{0, 2.3, 0}},
					{N: 1, Forces: [3]float64{10, 0, 0}},
					{N: 3, Forces: [3]float64{0, 2.3, 0}},
					{N: 3, Forces: [3]float64{10, 0, 0}},
				},
			},
		},
		ModalCases: []ModalCase{
			{
				ModalMasses: []ModalMass{
					{N: 1, Mass: 10000},
					{N: 3, Mass: 10000},
				},
			},
		},
	}
}

func baseTruss() Model {
	return Model{
		Points: [][2]float64{
			{0.0, 0.0}, // 1
			{0.0, 12.}, // 2
			{4.0, 0.0}, // 3
			{4.0, 6.0}, // 4
			{8.0, 0.0}, // 5
		},
		Beams: []BeamProp{
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
		LoadCases: []LoadCase{
			{
				LoadNodes: []LoadNode{
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
		ModalCases: []ModalCase{
			{
				ModalMasses: []ModalMass{
					{
						N:    1,
						Mass: 10000,
					},
				},
			},
		},
	}
}

func TestWrongLoad(t *testing.T) {
	m := baseTruss()
	m.LoadCases[0].LoadNodes[0].Forces[2] = 42.0 // not acceptable moment on pin
	var buf bytes.Buffer
	if err := m.Run(&buf); err == nil {
		t.Fatalf("not acceptable moment on pin")
	}
}

func baseModalTruss() Model {
	A := math.Pi * math.Pow(0.050, 2) / 4.0
	J := math.Pi * math.Pow(0.050, 4) / 64.0
	m := Model{
		Points: [][2]float64{
			{0.000, 0.0}, // 0
			{0.400, 0.0}, // 1
			{1.000, 0.0}, // 2
		},
		Beams: []BeamProp{
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
		ModalCases: []ModalCase{
			{
				ModalMasses: []ModalMass{
					{
						N:    1,
						Mass: 100 * Gravity,
					},
				},
			},
		},
	}
	return m
}

func baseModalTrussRotate() Model {
	A := math.Pi * math.Pow(0.050, 2) / 4.0
	J := math.Pi * math.Pow(0.050, 4) / 64.0
	m := Model{
		Points: [][2]float64{
			{0.0, 0.000}, // 0
			{0.0, 0.400}, // 1
			{0.0, 1.000}, // 2
		},
		Beams: []BeamProp{
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
		ModalCases: []ModalCase{
			{
				ModalMasses: []ModalMass{
					{
						N:    1,
						Mass: 100 * Gravity,
					},
				},
			},
		},
	}
	return m
}

func baseModalBeam() Model {
	A := math.Pi * math.Pow(0.050, 2) / 4.0
	J := math.Pi * math.Pow(0.050, 4) / 64.0
	m := Model{
		Points: [][2]float64{
			{0.000, 0.0}, // 0
			{0.400, 0.0}, // 1
			{1.000, 0.0}, // 2
		},
		Beams: []BeamProp{
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
		ModalCases: []ModalCase{
			{
				ModalMasses: []ModalMass{
					{
						N:    1,
						Mass: 100 * Gravity,
					},
				},
			},
		},
	}
	return m
}

func baseModalBeamRotate() Model {
	A := math.Pi * math.Pow(0.050, 2) / 4.0
	J := math.Pi * math.Pow(0.050, 4) / 64.0
	m := Model{
		Points: [][2]float64{
			{0.0, 0.000}, // 0
			{0.0, 0.400}, // 1
			{0.0, 1.000}, // 2
		},
		Beams: []BeamProp{
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
		ModalCases: []ModalCase{
			{
				ModalMasses: []ModalMass{
					{
						N:    1,
						Mass: 100 * Gravity,
					},
				},
			},
		},
	}
	return m
}

func baseModalBeam3mass() Model {
	E := 2e11
	J := 15e-4
	m := 250.0
	l := 4.0
	A := 12e-2
	return Model{
		Points: [][2]float64{
			{0.0, 0.000},     // 0
			{0.0, l / 6.0},   // 1
			{0.0, l / 2.0},   // 2
			{0.0, l - l/6.0}, // 3
			{0.0, l},         // 4
		},
		Beams: []BeamProp{
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
		ModalCases: []ModalCase{
			{
				ModalMasses: []ModalMass{
					{N: 1, Mass: m * Gravity},
					{N: 2, Mass: m * Gravity},
					{N: 3, Mass: m * Gravity},
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
func baseBeamDc() (m1, m2 Model) {
	return Model{
			Points: [][2]float64{
				{0.0, 0.0},
				{1.0, 0.0},
				{2.0, 0.0},
				{3.0, 0.0},
			},
			Beams: []BeamProp{
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
			LoadCases: []LoadCase{
				{
					LoadNodes: []LoadNode{
						{N: 1, Forces: [3]float64{0.0, 10.0, 0.0}},
						{N: 2, Forces: [3]float64{0.0, 10.0, 0.0}},
					},
				},
			},
		}, Model{
			Points: [][2]float64{
				{0.0, 0.0},
				{1.0, 0.0},
				{1.5, 0.0},
				{2.0, 0.0},
				{3.0, 0.0},
			},
			Beams: []BeamProp{
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
			LoadCases: []LoadCase{
				{
					LoadNodes: []LoadNode{
						{N: 1, Forces: [3]float64{0.0, 10.0, 0.0}},
						{N: 3, Forces: [3]float64{0.0, 10.0, 0.0}},
					},
				},
			},
		}
}

func TestJsonModel(t *testing.T) {
	m := baseBeam()
	b, err := json.Marshal(m)
	if err != nil {
		t.Fatal(err)
	}

	var mu Model
	err = json.Unmarshal(b, &mu)
	if err != nil {
		t.Fatal(err)
	}

	bu, err := json.Marshal(m)
	if err != nil {
		t.Fatal(err)
	}

	if !bytes.Equal(b, bu) {
		t.Fatalf("Is not same:\n%s\n%s", b, bu)
	}
}

func TestModelFail(t *testing.T) {
	ms := []Model{
		// error of input data
		{
			Points: [][2]float64{
				{-math.MaxFloat64 - 1, math.NaN()},
				{0.0, 0.0},
				{0.0, math.Inf(0)},
			},
			Beams: []BeamProp{
				{N: [2]int{-1, 1}, A: -12e-4, J: -120e-6, E: -2.0e11},
				{N: [2]int{1, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
				{N: [2]int{1, 20}, A: 12e-4, J: 120e-6, E: math.Inf(-1)},
			},
			Supports: [][3]bool{
				{true, true, true},
				{false, false, false},
				{false, true, false},
				{false, false, false},
				{false, false, false},
				{false, false, false},
				{false, false, false},
				{false, false, false},
				{false, false, false},
			},
			Pins: [][6]bool{
				{true, true, true, true, true, false},
				{true, false, false, true, false, false},
			},
			LoadCases: []LoadCase{
				{
					LoadNodes: []LoadNode{
						{N: -1, Forces: [3]float64{0, 2.3, 0}},
						{N: 5, Forces: [3]float64{math.Inf(1), 0, math.NaN()}},
					},
				},
			},
			ModalCases: []ModalCase{
				{ModalMasses: []ModalMass{{N: 7, Mass: -100}}},
				{ModalMasses: []ModalMass{{N: -1, Mass: math.NaN()}}},
				{ModalMasses: []ModalMass{{N: 0, Mass: math.NaN()}}},
				{}, // Modal mass is empty
			},
		},
		// error of inpossible solving - all points are free. No fix supports
		{
			Points: [][2]float64{
				{0, 0},
				{0, 1},
				{1, 1},
			},
			Beams: []BeamProp{
				{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
				{N: [2]int{1, 2}, A: 12e-4, J: 120e-6, E: 2.0e11},
			},
			Supports: [][3]bool{
				{false, false, false},
				{false, false, false},
				{false, false, false},
			},
			LoadCases: []LoadCase{
				{
					LoadNodes: []LoadNode{
						{N: 1, Forces: [3]float64{0, 2.3, 0}},
					},
				},
			},
			ModalCases: []ModalCase{
				{ModalMasses: []ModalMass{{N: 1, Mass: 100}}},
			},
		},
	}
	for i := range ms {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			var b bytes.Buffer
			err := ms[i].Run(&b)
			if err == nil {
				t.Fatalf("Error : %v", err)
			}
			t.Log(err)
		})
	}
}

func TestSplit(t *testing.T) {
	models := []func() Model{
		baseBeam, baseDoubleBeam,
		baseGBeam,
		baseTruss,
		baseModalBeam, baseModalBeamRotate,
		baseModalTruss, baseModalTrussRotate,
	}
	for mIndex := range models {
		t.Run(fmt.Sprintf("Model%d", mIndex), func(t *testing.T) {
			m := models[mIndex]()
			var b bytes.Buffer
			if err := m.Run(&b); err != nil {
				t.Fatalf("Error : %v", err)
			}
			expectResult := m.LoadCases
			hz := m.ModalCases[0].Result[0].Hz

			for i := 1; i < 10; i++ {
				t.Run(fmt.Sprintf("Split%d", i), func(t *testing.T) {
					mLocal := models[mIndex]()

					// split each beams
					var b bytes.Buffer
					amountBeams := len(mLocal.Beams)
					for j := 0; j < amountBeams; j++ {
						if err := mLocal.SplitBeam(j, i); err != nil {
							t.Fatalf("Cannot split %d: %v", i, err)
						}
					}

					// calculation
					if err := mLocal.Run(&b); err != nil {
						t.Fatalf("Error : %v", err)
					}

					// eps
					eps := 1e-9

					if err := compare(expectResult, mLocal.LoadCases, eps); err != nil {
						t.Log(mLocal.String())
						t.Errorf("Result is not same: %v", err)
					}

					h := mLocal.ModalCases[0].Result[0].Hz
					diff := math.Abs((hz - h) / h)
					if diff > eps {
						t.Logf("Narural frequency: %15.5e != %15.5e", hz, h)
						t.Errorf("Diff in natural frequency is not ok : %15.5e", diff)
					}
				})
			}
		})
	}
}

func compare(exp, act []LoadCase, eps float64) error {

	// compare reactions
	for e := range exp {
		for r := range exp[e].Reactions {
			expReact := exp[e].Reactions[r]
			actReact := act[e].Reactions[r]
			for j := 0; j < 3; j++ {
				if expReact[j] == 0.0 {
					if math.Abs(actReact[j]) > eps {
						return fmt.Errorf("expReact is zero, actReact != 0")
					}
					continue
				}
				diff := math.Abs((expReact[j] - actReact[j]) / expReact[j])
				if diff > eps {
					return fmt.Errorf("Diff[%d] is not ok : %15.5e", j, diff)
				}
			}
		}
	}

	// compare displacement
	for e := range exp {
		for d := range exp[e].PointDisplacementGlobal {
			expDisp := exp[e].PointDisplacementGlobal[d]
			actDisp := act[e].PointDisplacementGlobal[d]
			for j := 0; j < 3; j++ {
				if expDisp[j] == 0.0 {
					if math.Abs(actDisp[j]) > eps {
						return fmt.Errorf("expDisp is zero, actDisp != 0\n"+
							"actDisp = %10.8e", actDisp[j])
					}
					continue
				}
				diff := math.Abs((expDisp[j] - actDisp[j]) / expDisp[j])
				if diff > eps {
					return fmt.Errorf("Diff[%d] is not ok : %15.5e", j, diff)
				}
			}
		}
	}

	return nil
}

func TestSplitFail(t *testing.T) {
	m := baseBeam()
	tcs := []struct {
		beamIndex, amounts int
	}{
		{-1, 5},
		{100, 5},
		{0, -1},
	}
	for i, tc := range tcs {
		t.Run(fmt.Sprintf("Case #%d", i), func(t *testing.T) {
			if err := m.SplitBeam(tc.beamIndex, tc.amounts); err == nil {
				t.Errorf("Error for %v", tc)
			}
		})
	}
}

func TestCodeStyle(t *testing.T) {
	cs.All(t)
}

func TestWriter(t *testing.T) {
	m := baseTruss()
	tf, err := ioutil.TempFile("", "testWriter")
	if err != nil {
		t.Fatalf("Cannot create temp file: %v", err)
	}
	old := os.Stdout
	os.Stdout = tf
	defer func() {
		os.Stdout = old
	}()
	if err := m.Run(nil); err != nil {
		t.Fatalf("Cannot calculate : %v", err)
	}
	if err := tf.Close(); err != nil {
		t.Fatalf("Cannot close file : %v", err)
	}
	b, err := ioutil.ReadFile(tf.Name())
	if err != nil {
		t.Fatalf("Cannot read file : %v", err)
	}
	if len(b) == 0 {
		t.Fatalf("File is empty")
	}
}

func TestModelString(t *testing.T) {
	ms := []struct {
		m        Model
		filename string
	}{{
		m:        baseBeam(),
		filename: "beam",
	}, {
		m:        baseTruss(),
		filename: "truss",
	}, {
		m:        baseGBeam(),
		filename: "g-beam",
	}, {
		m:        baseDoubleBeam(),
		filename: "double-beam",
	}, {
		m:        baseModalTruss(),
		filename: "truss-modal",
	}, {
		m:        baseModalTrussRotate(),
		filename: "truss-modal-rotate",
	}, {
		m:        baseModalBeam(),
		filename: "beam-modal",
	}, {
		m:        baseModalBeamRotate(),
		filename: "beam-modal-rotate",
	}, {
		m:        baseModalBeam3mass(),
		filename: "beam-modal-3mass",
	}, {
		m: func() Model {
			m, _ := baseBeamDc()
			return m
		}(),
		filename: "beam-dc-part1",
	}, {
		m: func() Model {
			_, m := baseBeamDc()
			return m
		}(),
		filename: "beam-dc-part2",
	}}
	for _, m := range ms {
		t.Run(m.filename, func(t *testing.T) {
			var b bytes.Buffer
			if err := m.m.Run(&b); err != nil {
				t.Errorf("Cannot calculate : %v", err)
			}
			b.Reset()

			// compare files
			actual := []byte(m.m.String())

			if os.Getenv("UPDATE") != "" {
				err := ioutil.WriteFile("./testdata/"+m.filename, actual, 0644)
				if err != nil {
					t.Fatalf("Cannot Update: %v", err)
				}
			}

			expect, err := ioutil.ReadFile("./testdata/" + m.filename)
			if err != nil {
				t.Fatalf("Cannot read file : %v", err)
			}

			if !bytes.Equal(expect, actual) {
				// show a diff between files
				diff := difflib.UnifiedDiff{
					A:        difflib.SplitLines(string(expect)),
					B:        difflib.SplitLines(string(actual)),
					FromFile: "Original",
					ToFile:   "Current",
					Context:  30000,
				}
				text, _ := difflib.GetUnifiedDiffString(diff)
				t.Log(text)
				t.Errorf("result is not same")
			}
		})
	}
}

func TestDirectionLoadNode(t *testing.T) {
	m := Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{1.0, 0.0},
		},
		Beams: []BeamProp{
			{
				N: [2]int{0, 1},
				A: 12e-4,
				J: 120e-6,
				E: 2.0e11,
			},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
		},
		LoadCases: []LoadCase{
			{
				LoadNodes: []LoadNode{
					{N: 1, Forces: [3]float64{100.0, 0.0, 0.0}},
				},
			},
			{
				LoadNodes: []LoadNode{
					{N: 1, Forces: [3]float64{0.0, 100.0, 0.0}},
				},
			},
			{
				LoadNodes: []LoadNode{
					{N: 1, Forces: [3]float64{0.0, 0.0, 100.0}},
				},
			},
		},
	}

	var b bytes.Buffer
	if err := m.Run(&b); err != nil {
		t.Fatalf("Cannot calculate : %v", err)
	}
	b.Reset()

	if !(m.LoadCases[0].PointDisplacementGlobal[1][0] > 0) {
		t.Errorf("load in direction X is not ok. See: %v", m.LoadCases[0].PointDisplacementGlobal[1])
	}
	if !(m.LoadCases[1].PointDisplacementGlobal[1][1] > 0) {
		t.Errorf("load in direction Y is not ok. See: %v", m.LoadCases[1].PointDisplacementGlobal[1])
	}
	if !(m.LoadCases[2].PointDisplacementGlobal[1][1] > 0) {
		t.Errorf("load in direction M is not ok. See: %v", m.LoadCases[2].PointDisplacementGlobal[1])
	}
	t.Logf("case 0 : %v", m.LoadCases[0].PointDisplacementGlobal)
	t.Logf("case 1 : %v", m.LoadCases[1].PointDisplacementGlobal)
	t.Logf("case 2 : %v", m.LoadCases[2].PointDisplacementGlobal)
}

func TestLoadUniform(t *testing.T) {
	t.Run("error checking", func(t *testing.T) {
		// error checking
		errs := []error{
			func() error {
				var m *Model
				_, err := m.LoadUniform(0, false, [2]float64{0, 0})
				return err
			}(),
			func() error {
				m := Model{}
				_, err := m.LoadUniform(0, false, [2]float64{0, 0})
				return err
			}(),
			func() error {
				m := Model{}
				_, err := m.LoadUniform(-1, false, [2]float64{0, 0})
				return err
			}(),
			func() error {
				m := Model{}
				_, err := m.LoadUniform(100, false, [2]float64{0, 0})
				return err
			}(),
			func() error {
				m := Model{
					Beams: []BeamProp{
						{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
					},
				}
				_, err := m.LoadUniform(1, false, [2]float64{0, 0})
				return err
			}(),
			func() error {
				m := Model{
					Beams: []BeamProp{
						{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
					},
				}
				_, err := m.LoadUniform(0, false, [2]float64{math.Inf(0), math.NaN()})
				return err
			}(),
		}

		for i := range errs {
			if errs[i] == nil {
				t.Errorf("case %d is not correct", i)
			}
			t.Log(errs[i])
		}
	})

	// calculation checking
	t.Run("check on beam with part load", func(t *testing.T) {
		for _, proj := range []bool{false, true} {
			for _, mirror := range []bool{false, true} {
				for _, d := range []float64{200, -0.5, 0.5, -200, 0} {
					size := 3
					t.Run(fmt.Sprintf("Mirror:%v/Proj:%v/%3.1f/Size:%d", mirror, proj, d, size), func(t *testing.T) {
						m := Model{
							Points: [][2]float64{
								{0.0, 0.0},
								{1.0, d},
								{2.0, 0.0},
							},
							Beams: []BeamProp{
								{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
								{N: [2]int{1, 2}, A: 12e-4, J: 120e-6, E: 2.0e11},
							},
							Supports: [][3]bool{
								{true, true, true},
								{false, false, false},
								{false, true, false},
							},
							LoadCases: []LoadCase{
								{LoadNodes: []LoadNode{{N: 1, Forces: [3]float64{100.0, 0.0, 0.0}}}},
							},
						}

						if mirror {
							m.Points[0][0] = 2.0
							m.Points[2][0] = 0.0
						}

						if err := m.SplitBeam(0, size); err != nil {
							t.Fatal(err)
						}
						if err := m.SplitBeam(1, size); err != nil {
							t.Fatal(err)
						}

						// reset and allocate memory
						m.LoadCases = make([]LoadCase, 1)
						// loads
						ux := -10.0
						uy := -25.0
						// uniform load
						for i := range m.Beams {
							if m.Points[m.Beams[i].N[0]][0] > 1.0 || m.Points[m.Beams[i].N[1]][0] > 1.0 {
								// not add on right part of beam
								continue
							}
							un, err := m.LoadUniform(i, proj, [2]float64{ux, uy})
							if err != nil {
								t.Fatal(err)
							}
							m.LoadCases[0].LoadNodes = append(m.LoadCases[0].LoadNodes, un...)
						}
						// node load
						mn := m
						mn.LoadCases = make([]LoadCase, 1)
						dx := math.Abs(m.Points[0][0] - m.Points[1][0])
						dy := math.Abs(m.Points[0][1] - m.Points[1][1])
						if !proj {
							// not projection
							dx = math.Sqrt(dx*dx + dy*dy)
							dy = dx
						}
						for i := range mn.Points {
							if m.Points[i][0] > 1.0 {
								// not add on right part of beam
								continue
							}
							if m.Points[i][0] == 0.0 || m.Points[i][0] == 1.0 {
								mn.LoadCases[0].LoadNodes = append(mn.LoadCases[0].LoadNodes, LoadNode{
									N: i,
									Forces: [3]float64{
										ux * dy / float64(size+1) / 2.0, // X
										uy * dx / float64(size+1) / 2.0, // Y
										0.0,                             // M
									},
								})
								continue
							}
							mn.LoadCases[0].LoadNodes = append(mn.LoadCases[0].LoadNodes, LoadNode{
								N: i,
								Forces: [3]float64{
									ux * dy / float64(size+1), // X
									uy * dx / float64(size+1), // Y
									0.0,                       // M
								},
							})
						}

						// calculation
						var buf bytes.Buffer
						if err := m.Run(&buf); err != nil {
							t.Fatalf("m model error : %v", err)
						}
						buf.Reset()
						if err := mn.Run(&buf); err != nil {
							t.Fatalf("mn model error : %v", err)
						}
						buf.Reset()

						// comparing displacement
						eps := 0.05 // 5%
						var actual, sum float64
						for i := 0; i < 3; i++ {
							actual += math.Pow(
								m.LoadCases[0].PointDisplacementGlobal[1][i]-
									mn.LoadCases[0].PointDisplacementGlobal[1][i],
								2)
							sum += math.Pow(m.LoadCases[0].PointDisplacementGlobal[1][i], 2)
						}
						actual = math.Sqrt(actual)
						sum = math.Sqrt(sum)
						if actual/sum >= eps {
							t.Log(actual / sum)
							t.Log(m.LoadCases[0].PointDisplacementGlobal[1])
							t.Log(mn.LoadCases[0].PointDisplacementGlobal[1])
							t.Errorf("displacement precision of calculation is not ok: %10.5f >= 0.05", actual/sum)
						}

						// comparing reactions
						actual = 0.0
						sum = 0.0
						for i := 0; i < 3; i++ {
							actual += math.Pow(
								m.LoadCases[0].Reactions[0][i]-
									mn.LoadCases[0].Reactions[0][i],
								2)
							sum += math.Pow(mn.LoadCases[0].Reactions[0][i], 2)
						}
						actual = math.Sqrt(actual)
						sum = math.Sqrt(sum)
						if actual/sum >= eps {
							t.Log(actual / sum)
							t.Log(m.LoadCases[0].Reactions[0])
							t.Log(mn.LoadCases[0].Reactions[0])
							t.Errorf("reaction     precision of calculation is not ok: %10.5f >= 0.05", actual/sum)
						}
					})
				}
			}
		}
	})
}

func BenchmarkRun(b *testing.B) {
	minimalLoadCases := 20
	for ic := 1; ic <= 128; ic *= 2 {
		b.Run(fmt.Sprintf("%5d-cases%d", ic, minimalLoadCases), func(b *testing.B) {
			// prepare model
			m := baseBeam()
			// add more finite elements
			err := m.SplitBeam(0, ic)
			if err != nil {
				panic(err)
			}
			// add more cases
			for i := 0; len(m.LoadCases) < minimalLoadCases; i++ {
				lc := m.LoadCases[0]
				lc.LoadNodes[0].Forces[1] += float64(i)
				m.LoadCases = append(m.LoadCases, lc)
			}
			for i := 2; i < len(m.Points); i += 2 {
				m.ModalCases[0].ModalMasses = append(m.ModalCases[0].ModalMasses,
					ModalMass{N: i, Mass: 100})
			}
			var bb bytes.Buffer
			b.ResetTimer()
			// run benchmark
			for i := 0; i < b.N; i++ {
				err := m.Run(&bb)
				if err != nil {
					panic(err)
				}
				bb.Reset()
			}
		})
	}
}
