package hd

import (
	"bufio"
	"bytes"
	"encoding/json"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func baseModel() Model {
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
					{
						N:      1,
						Forces: [3]float64{0, 2.3, 0},
					},
					{
						N:      1,
						Forces: [3]float64{10, 0, 0},
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
				N: [2]int{0, 1},
				A: 40e-4,
				J: 1,
				E: 2.0e11,
			}, { // 2
				N: [2]int{0, 2},
				A: 64e-4,
				J: 1,
				E: 2.0e11,
			}, { // 3
				N: [2]int{0, 3},
				A: 60e-4,
				J: 1,
				E: 2.0e11,
			}, { // 4
				N: [2]int{1, 3},
				A: 60e-4,
				J: 1,
				E: 2.0e11,
			}, { // 5
				N: [2]int{2, 3},
				A: 40e-4,
				J: 1,
				E: 2.0e11,
			}, { // 6
				N: [2]int{2, 4},
				A: 64e-4,
				J: 1,
				E: 2.0e11,
			}, { // 7
				N: [2]int{3, 4},
				A: 60e-4,
				J: 1,
				E: 2.0e11,
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

func TestJsonModel(t *testing.T) {
	m := baseModel()
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
	models := []Model{baseModel(), baseTruss()}
	for mindex, m := range models {
		var b bytes.Buffer
		if err := m.Run(&b); err != nil {
			t.Fatalf("Error : %v", err)
		}
		reaction := m.LoadCases[0].Reactions[0]
		displacament := m.LoadCases[0].PointDisplacementGlobal[1]
		hz := m.ModalCases[0].Result[0].Hz

		// TODO: if i = 1, then Truss model is error with singular problem
		for i := 2; i < 10; i++ {
			t.Run(fmt.Sprintf("Model%d Split%d", mindex, i), func(t *testing.T) {
				mlocal := models[mindex]
				// remove results from last iteration
				for lc := 0; lc < len(mlocal.LoadCases); lc++ {
					mlocal.LoadCases[lc].Reactions = [][3]float64{}
					mlocal.LoadCases[lc].BeamForces = [][6]float64{}
					mlocal.LoadCases[lc].PointDisplacementGlobal = [][3]float64{}
				}
				for mc := 0; mc < len(mlocal.ModalCases); mc++ {
					mlocal.ModalCases[mc].Result = []ModalResult{}
				}

				// split each beams
				var b bytes.Buffer
				amountBeams := len(mlocal.Beams)
				for j := 0; j < amountBeams; j++ {
					if err := mlocal.SplitBeam(j, i); err != nil {
						t.Fatalf("Cannot split %d: %v", i, err)
					}
				}

				// calculation
				if err := mlocal.Run(&b); err != nil {
					t.Fatalf("Error : %v", err)
				}

				// eps
				eps := 1e-9

				// compare results
				r := mlocal.LoadCases[0].Reactions[0]
				for j := 0; j < 3; j++ {
					diff := math.Abs((reaction[j] - r[j]) / r[j])
					t.Logf("Reactions : %15.5e != %15.5e , %15.5e",
						reaction[j], r[j], diff)
					if diff > eps {
						t.Errorf("Diff[%d] is not ok : %15.5e", j, diff)
					}
				}

				d := mlocal.LoadCases[0].PointDisplacementGlobal[1]
				for j := 0; j < 3; j++ {
					diff := math.Abs((displacament[j] - d[j]) / d[j])
					if diff > eps {
						t.Logf("Displacament: %15.5e != %15.5e", displacament[j], d[j])
						t.Errorf("Diff[%d] is not ok : %15.5e", j, diff)
					}
				}

				h := mlocal.ModalCases[0].Result[0].Hz
				diff := math.Abs((hz - h) / h)
				if diff > eps {
					t.Logf("Narural frequency: %15.5e != %15.5e", hz, h)
					t.Errorf("Diff in natural frequency is not ok : %15.5e", diff)
				}
			})
		}
	}
}

func TestSplitFail(t *testing.T) {
	m := baseModel()
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

func TestTodo(t *testing.T) {
	// Show all to do`s in comment code
	source, err := filepath.Glob(fmt.Sprintf("./%s", "*.go"))
	if err != nil {
		t.Fatal(err)
	}

	var amount int

	for i := range source {
		t.Run(source[i], func(t *testing.T) {
			file, err := os.Open(source[i])
			if err != nil {
				t.Fatal(err)
			}
			defer file.Close()

			pos := 1
			scanner := bufio.NewScanner(file)
			for scanner.Scan() {
				line := scanner.Text()
				pos++
				index := strings.Index(line, "//")
				if index < 0 {
					continue
				}
				if !strings.Contains(strings.ToUpper(line), "TODO") {
					continue
				}
				t.Logf("%d %s", pos, line[index:])
				amount++
			}

			if err := scanner.Err(); err != nil {
				log.Fatal(err)
			}
		})
	}
	if amount > 0 {
		t.Logf("Amount TODO: %d", amount)
	}
}

func TestFmt(t *testing.T) {
	// Show all fmt`s in comments code
	source, err := filepath.Glob(fmt.Sprintf("./%s", "*.go"))
	if err != nil {
		t.Fatal(err)
	}

	var amount int

	for i := range source {
		t.Run(source[i], func(t *testing.T) {
			file, err := os.Open(source[i])
			if err != nil {
				t.Fatal(err)
			}
			defer file.Close()

			pos := 1
			scanner := bufio.NewScanner(file)
			for scanner.Scan() {
				line := scanner.Text()
				pos++
				index := strings.Index(line, "//")
				if index < 0 {
					continue
				}
				if !strings.Contains(line, "fmt.") {
					continue
				}
				t.Logf("%d %s", pos, line[index:])
				amount++
			}

			if err := scanner.Err(); err != nil {
				log.Fatal(err)
			}
		})
	}
	if amount > 0 {
		t.Logf("Amount commented fmt: %d", amount)
	}
}

func TestDebug(t *testing.T) {
	// Show all debug information in code
	source, err := filepath.Glob(fmt.Sprintf("./%s", "*.go"))
	if err != nil {
		t.Fatal(err)
	}

	for i := range source {
		t.Run(source[i], func(t *testing.T) {
			file, err := os.Open(source[i])
			if err != nil {
				t.Fatal(err)
			}
			defer file.Close()

			pos := 1
			scanner := bufio.NewScanner(file)
			for scanner.Scan() {
				line := scanner.Text()
				pos++
				if !strings.Contains(line, "fmt"+"."+"Print") {
					continue
				}
				t.Errorf("Debug line: %d in file %s", pos, source[i])
			}

			if err := scanner.Err(); err != nil {
				log.Fatal(err)
			}
		})
	}
}

func ExampleString() {
	ms := []Model{baseModel(), baseTruss()}
	for _, m := range ms {
		var b bytes.Buffer
		if err := m.Run(&b); err != nil {
			continue
		}
		b.Reset()
		fmt.Fprintln(os.Stdout, m.String())
	}

	// Output:

	// Point coordinates:
	// Index            X, m            Y, m    SX    SY    SM (0 - free, 1 - fixed)
	//     0         0.00000         0.00000     1     1     1
	//     1         1.00000         0.00000     0     0     0
	// Beam property:
	// Index     Start point       End point       Area,sq.m Moment inertia,m4   Elasticity,Pa
	//     0               0               1     1.20000e-03     1.20000e-04     2.00000e+11
	// All beams haven`t pins
	//
	// Load case #  0
	// Point           Fx, N           Fy, N          M, N*m
	//     1         0.00000         2.30000         0.00000
	//     1        10.00000         0.00000         0.00000
	// Point displacament in global system coordinate:
	// Point           DX, m           DY, m
	//     0     0.00000e+00     0.00000e+00
	//     1     4.16667e-08     3.19444e-08
	// Local force in beam:
	// Index           Fx, N           Fy, N          M, N*m           Fx, N           Fy, N          M, N*m
	//     0    -1.00000e+01    -2.30000e+00    -2.30000e+00     1.00000e+01     2.30000e+00     8.88178e-16
	// Reaction on support:
	// Index           Fx, N           Fy, N          M, N*m
	//     0    -1.00000e+01    -2.30000e+00    -2.30000e+00
	//
	// Modal case #  0
	// Point         Mass, N
	//     1     10000.00000
	// Natural frequency :        42.29088 Hz
	// Point               X               Y               M
	//     0     0.00000e+00     0.00000e+00     0.00000e+00
	//     1     0.00000e+00     5.54700e-01     8.32050e-01
	// Natural frequency :        77.21223 Hz
	// Point               X               Y               M
	//     0     0.00000e+00     0.00000e+00     0.00000e+00
	//     1     1.00000e+00     0.00000e+00     0.00000e+00
	//
	//
	// Point coordinates:
	// Index            X, m            Y, m    SX    SY    SM (0 - free, 1 - fixed)
	//     0         0.00000         0.00000     1     1     0
	//     1         0.00000        12.00000     0     0     0
	//     2         4.00000         0.00000     0     1     0
	//     3         4.00000         6.00000     0     0     0
	//     4         8.00000         0.00000     0     1     0
	// Beam property:
	// Index     Start point       End point       Area,sq.m Moment inertia,m4   Elasticity,Pa
	//     0               0               1     4.00000e-03     1.00000e+00     2.00000e+11
	//     1               0               2     6.40000e-03     1.00000e+00     2.00000e+11
	//     2               0               3     6.00000e-03     1.00000e+00     2.00000e+11
	//     3               1               3     6.00000e-03     1.00000e+00     2.00000e+11
	//     4               2               3     4.00000e-03     1.00000e+00     2.00000e+11
	//     5               2               4     6.40000e-03     1.00000e+00     2.00000e+11
	//     6               3               4     6.00000e-03     1.00000e+00     2.00000e+11
	// Pins of beam in local system coordinate:
	// Index       X       Y       M       X       Y       M
	//     0   false   false    true   false   false    true
	//     1   false   false    true   false   false    true
	//     2   false   false    true   false   false    true
	//     3   false   false    true   false   false    true
	//     4   false   false    true   false   false    true
	//     5   false   false    true   false   false    true
	//     6   false   false    true   false   false    true
	//
	// Load case #  0
	// Point           Fx, N           Fy, N          M, N*m
	//     1    -70000.00000         0.00000         0.00000
	//     3     42000.00000         0.00000         0.00000
	// Point displacament in global system coordinate:
	// Point           DX, m           DY, m
	//     0     0.00000e+00     0.00000e+00
	//     1    -4.61042e-03    -1.57500e-03
	//     2    -1.06771e-04     0.00000e+00
	//     3    -3.80192e-04     3.33751e-04
	//     4    -2.13541e-04     0.00000e+00
	// Local force in beam:
	// Index           Fx, N           Fy, N          M, N*m           Fx, N           Fy, N          M, N*m
	//     0     1.05000e+05     0.00000e+00     0.00000e+00    -1.05000e+05     0.00000e+00     0.00000e+00
	//     1     3.41666e+04     0.00000e+00     0.00000e+00    -3.41666e+04     0.00000e+00     0.00000e+00
	//     2    -1.11170e+04     9.56479e-10     0.00000e+00     1.11170e+04    -9.56479e-10     0.00000e+00
	//     3    -1.26194e+05     8.73289e-09     0.00000e+00     1.26194e+05    -8.73289e-09     0.00000e+00
	//     4    -4.45001e+04     0.00000e+00     0.00000e+00     4.45001e+04     0.00000e+00     0.00000e+00
	//     5     3.41666e+04     0.00000e+00     0.00000e+00    -3.41666e+04     0.00000e+00     0.00000e+00
	//     6    -6.15948e+04    -8.86350e-11     0.00000e+00     6.15948e+04     8.86350e-11     0.00000e+00
	// Reaction on support:
	// Index           Fx, N           Fy, N          M, N*m
	//     0     2.80000e+04    -9.57501e+04     0.00000e+00
	//     2     0.00000e+00     4.45001e+04     0.00000e+00
	//     4     0.00000e+00     5.12499e+04     0.00000e+00
	//
	// Modal case #  0
	// Point         Mass, N
	//     1     10000.00000
	// Natural frequency :        18.42543 Hz
	// Point               X               Y               M
	//     0     0.00000e+00     0.00000e+00     0.00000e+00
	//     1     9.39632e-01     2.88945e-01     0.00000e+00
	//     2     3.09602e-02     0.00000e+00     0.00000e+00
	//     3     1.56363e-01    -6.60315e-02     0.00000e+00
	//     4     6.19204e-02     0.00000e+00     0.00000e+00
}

func BenchmarkRun(b *testing.B) {
	tcs := []int{1, 2, 4, 8, 16, 32, 64, 128, 256, 512}
	for _, tc := range tcs {
		b.Run(fmt.Sprintf("%5d", tc), func(b *testing.B) {
			m := baseModel()
			var bb bytes.Buffer
			_ = m.SplitBeam(0, tc)
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				_ = m.Run(&bb)
				bb.Reset()
			}
		})
	}
}
