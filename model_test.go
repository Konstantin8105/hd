package hd

import (
	"bufio"
	"bytes"
	"encoding/json"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/bradleyjkemp/cupaloy"
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
	m := Model{
		Points: [][2]float64{
			{-math.MaxFloat64 - 1, math.NaN()},
			{0.0, 0.0},
			{0.0, math.Inf(0)},
		},
		Beams: []BeamProp{
			{
				N: [2]int{-1, 1},
				A: -12e-4,
				J: -120e-6,
				E: -2.0e11,
			},
			{
				N: [2]int{1, 1},
				A: 12e-4,
				J: 120e-6,
				E: 2.0e11,
			},
			{
				N: [2]int{1, 20},
				A: 12e-4,
				J: 120e-6,
				E: math.Inf(-1),
			},
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
		LoadCases: []LoadCase{
			{
				LoadNodes: []LoadNode{
					{
						N:      -1,
						Forces: [3]float64{0, 2.3, 0},
					},
					{
						N:      5,
						Forces: [3]float64{math.Inf(1), 0, math.NaN()},
					},
				},
			},
		},
		ModalCases: []ModalCase{
			{
				ModalMasses: []ModalMass{
					{
						N:    7,
						Mass: -100,
					},
				},
			},
			{
				ModalMasses: []ModalMass{
					{
						N:    -1,
						Mass: math.NaN(),
					},
				},
			},
			{
				ModalMasses: []ModalMass{
					{
						N:    0,
						Mass: math.NaN(),
					},
				},
			},
		},
	}
	f, err := ioutil.TempFile("", "")
	if err != nil {
		t.Fatalf("cannot create temp file : %v", err)
	}
	stdout := os.Stdout
	os.Stdout = f
	defer func() {
		os.Stdout = stdout
	}()

	err = m.Run(nil)
	if err == nil {
		t.Fatalf("Error : %v", err)
	}
	t.Log(err)
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

		for i := 1; i < 10; i++ {
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

				// compare results
				r := mlocal.LoadCases[0].Reactions[0]
				for j := 0; j < 3; j++ {
					diff := math.Abs((reaction[j] - r[j]) / r[j])
					if diff > 1e-10 {
						t.Fatalf("Diff[%d] is not ok : %15.5e", j, diff)
					}
				}

				d := mlocal.LoadCases[0].PointDisplacementGlobal[1]
				for j := 0; j < 3; j++ {
					diff := math.Abs((displacament[j] - d[j]) / d[j])
					if diff > 1e-10 {
						t.Fatalf("Diff[%d] is not ok : %15.5e", j, diff)
					}
				}

				h := mlocal.ModalCases[0].Result[0].Hz
				diff := math.Abs((hz - h) / h)
				if diff > 1e-10 {
					t.Fatalf("Diff in natural frequency is not ok : %15.5e", diff)
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

func TestModelString(t *testing.T) {
	m := baseModel()
	if err := m.SplitBeam(0, 10); err != nil {
		t.Fatalf("Cannot split : %v", err)
	}
	var b bytes.Buffer
	if err := m.Run(&b); err != nil {
		t.Errorf("error : %v", err)
	}
	t.Log(m.String())
	if err := cupaloy.SnapshotMulti("TestModelString", m.String()); err != nil {
		t.Fatalf("error: %s", err)
	}
}

func TestTodo(t *testing.T) {
	// Show all todos in code
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
				if !strings.Contains(line, "TODO") {
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
	// Show all todos in code
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
				if !strings.Contains(line, "fmt") {
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
	// Show all todos in code
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

func TestTruss(t *testing.T) {
	m := baseTruss()
	var b bytes.Buffer
	if err := m.Run(&b); err != nil {
		t.Fatalf("Error : %v", err)
	}
	t.Log(m.String())
	if err := cupaloy.SnapshotMulti("Truss", m.String()); err != nil {
		t.Fatalf("error: %s", err)
	}
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
