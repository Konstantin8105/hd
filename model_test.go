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
			[2]float64{0.0, 0.0},
			[2]float64{1.0, 0.0},
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
			[3]bool{true, true, true},
			[3]bool{false, false, false},
		},
		LoadCases: []LoadCase{
			LoadCase{
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
			ModalCase{
				ModalMasses: []ModalMass{
					{
						N:    1,
						Mass: 10000,
					},
				},
			},
			// ModalCase{
			// 	ModalMasses: []ModalMass{
			// 		{
			// 			N: 1,
			// 			// Mass is zero
			// 		},
			// 	},
			// },
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
			[2]float64{-math.MaxFloat64 - 1, math.NaN()},
			[2]float64{0.0, 0.0},
			[2]float64{0.0, math.Inf(0)},
		},
		Beams: []BeamProp{
			{
				N: [2]int{-1, 1},
				A: -12e-4,
				J: -120e-6,
				E: -2.0e11,
			},
			{
				N: [2]int{1, 2},
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
			[3]bool{true, true, true},
			[3]bool{false, false, false},
			[3]bool{false, true, false},
			[3]bool{false, false, false},
			[3]bool{false, false, false},
			[3]bool{false, false, false},
			[3]bool{false, false, false},
			[3]bool{false, false, false},
			[3]bool{false, false, false},
		},
		LoadCases: []LoadCase{
			LoadCase{
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
			ModalCase{
				ModalMasses: []ModalMass{
					{
						N:    7,
						Mass: -100,
					},
				},
			},
			ModalCase{
				ModalMasses: []ModalMass{
					{
						N:    -1,
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
	m := baseModel()
	var b bytes.Buffer
	if err := m.Run(&b); err != nil {
		t.Fatalf("Error : %v", err)
	}
	reaction := m.LoadCases[0].Reactions[0]
	displacament := m.LoadCases[0].PointDisplacementGlobal[1]

	for i := 1; i < 10; i++ {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			m := baseModel()
			var b bytes.Buffer
			if err := m.SplitBeam(0, i); err != nil {
				t.Fatalf("Cannot split %d: %v", i, err)
			}
			if err := m.Run(&b); err != nil {
				t.Fatalf("Error : %v", err)
			}

			r := m.LoadCases[0].Reactions[0]
			for j := 0; j < 3; j++ {
				diff := math.Abs((reaction[j] - r[j]) / r[j])
				if diff > 1e-10 {
					t.Fatalf("Diff[%d] is not ok : %15.5e", j, diff)
				}
			}

			d := m.LoadCases[0].PointDisplacementGlobal[1]
			for j := 0; j < 3; j++ {
				diff := math.Abs((displacament[j] - d[j]) / d[j])
				if diff > 1e-10 {
					t.Fatalf("Diff[%d] is not ok : %15.5e", j, diff)
				}
			}
		})
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
