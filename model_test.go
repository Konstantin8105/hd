package hd_test

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"runtime/debug"
	"testing"

	"github.com/Konstantin8105/cs"
	"github.com/Konstantin8105/hd"
	"github.com/Konstantin8105/hd/example"
	"github.com/pmezard/go-difflib/difflib"
	"gonum.org/v1/gonum/mat"
)

func TestWrongLoad(t *testing.T) {
	m, lc, _ := example.Truss()
	lc[0].LoadNodes[0].Forces[2] = 42.0 // not acceptable moment on pin
	var buf bytes.Buffer
	if err := hd.LinearStatic(&buf, &m, &lc[0]); err == nil {
		t.Fatalf("not acceptable moment on pin")
	}
}

func TestJsonModel(t *testing.T) {
	m, _, _ := example.ConsoleBeam()
	b, err := json.Marshal(m)
	if err != nil {
		t.Fatal(err)
	}

	var mu hd.Model
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
	ms := []struct {
		m   hd.Model
		lcs []hd.LoadCase
		mcs []hd.ModalCase
	}{
		// error of input data
		{
			m: hd.Model{
				Points: [][2]float64{
					{-math.MaxFloat64 - 1, math.NaN()},
					{0.0, 0.0},
					{0.0, math.Inf(0)},
				},
				Beams: []hd.BeamProp{
					{N: [2]int{-1, 1}, A: -12e-4, J: -120e-6, E: -2.0e11},
					{N: [2]int{1, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
					{N: [2]int{1, 20}, A: 12e-4, J: math.NaN(), E: math.Inf(-1)},
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
			},
			lcs: []hd.LoadCase{
				{
					LoadNodes: []hd.LoadNode{
						{N: -1, Forces: [3]float64{0, 2.3, 0}},
						{N: 5, Forces: [3]float64{math.Inf(1), 0, math.NaN()}},
					},
					AmountLinearBuckling: 1,
				},
			},
			mcs: []hd.ModalCase{
				{ModalMasses: []hd.ModalMass{{N: 7, Mass: -100}}},
				{ModalMasses: []hd.ModalMass{{N: -1, Mass: math.NaN()}}},
				{ModalMasses: []hd.ModalMass{{N: 0, Mass: math.NaN()}}},
				{}, // Modal mass is empty
			},
		},
		// error of inpossible solving - all points are free. No fix supports
		{
			m: hd.Model{
				Points: [][2]float64{
					{0, 0},
					{0, 1},
					{1, 1},
				},
				Beams: []hd.BeamProp{
					{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
					{N: [2]int{1, 2}, A: 12e-4, J: 120e-6, E: 2.0e11},
				},
				Supports: [][3]bool{
					{false, false, false},
					{false, false, false},
					{false, false, false},
				},
			},
			lcs: []hd.LoadCase{
				{
					LoadNodes: []hd.LoadNode{
						{N: 1, Forces: [3]float64{0, 2.3, 0}},
					},
				},
			},
			mcs: []hd.ModalCase{
				{ModalMasses: []hd.ModalMass{{N: 1, Mass: 100}}},
			},
		},
		// error - too big model
		{
			m: hd.Model{
				Points: [][2]float64{
					{0, 0},
					{-1e300, 1e300},
					{1e300, 1e300},
				},
				Beams: []hd.BeamProp{
					{N: [2]int{0, 1}, A: 12e-300, J: 120e-300, E: 2.0e300},
					{N: [2]int{1, 2}, A: 12e+300, J: 120e+300, E: 2.0e300},
				},
				Supports: [][3]bool{
					{true, true, true},
					{false, false, false},
					{false, false, false},
				},
			},
			lcs: []hd.LoadCase{
				{
					LoadNodes: []hd.LoadNode{
						{N: 1, Forces: [3]float64{0, 2.3, 0}},
					},
				},
			},
			mcs: []hd.ModalCase{
				{ModalMasses: []hd.ModalMass{{N: 1, Mass: 100}}},
			},
		},
		// error - too big load
		{
			m: hd.Model{
				Points: [][2]float64{
					{0, 0},
					{0, 1},
					{1, 1},
				},
				Beams: []hd.BeamProp{
					{N: [2]int{0, 1}, A: 12e-30, J: 120e-30, E: 2.0e30},
					{N: [2]int{1, 2}, A: 12e-30, J: 120e-30, E: 2.0e30},
				},
				Supports: [][3]bool{
					{true, true, true},
					{false, false, false},
					{false, false, false},
				},
			},
			lcs: []hd.LoadCase{
				{
					LoadNodes: []hd.LoadNode{
						{N: 2, Forces: [3]float64{1e300, 1e308, 1e300}},
						{N: 2, Forces: [3]float64{1e300, 1e308, 1e300}},
					},
				},
			},
			mcs: []hd.ModalCase{
				{ModalMasses: []hd.ModalMass{{N: 1, Mass: 1e300}}},
			},
		},
		// error - too small load
		{
			m: hd.Model{
				Points: [][2]float64{
					{0, 0},
					{0, 1},
					{1, 1},
				},
				Beams: []hd.BeamProp{
					{N: [2]int{0, 1}, A: 12e-30, J: math.SmallestNonzeroFloat64, E: math.SmallestNonzeroFloat64},
					{N: [2]int{1, 2}, A: 12e-30, J: 120e-30, E: 2.0e30},
				},
				Pins: [][6]bool{
					{false, false, false, false, false, false},
					{false, false, true, false, false, true},
					{false, false, false, false, false, false},
				},
				Supports: [][3]bool{
					{true, true, true},
					{false, false, false},
					{false, false, false},
				},
			},
			lcs: []hd.LoadCase{
				{
					LoadNodes: []hd.LoadNode{
						{N: 20, Forces: [3]float64{0, 1, 0}},
						{N: -20, Forces: [3]float64{0, 1, 0}},
					},
				},
			},
			mcs: []hd.ModalCase{
				{ModalMasses: []hd.ModalMass{
					{N: 1, Mass: 1},
					{N: -1, Mass: -1},
					{N: 10000, Mass: 0},
				}},
			},
		},
		// error - rigid
		{
			m: hd.Model{
				Points: [][2]float64{
					{0, 0},
					{0, math.SmallestNonzeroFloat64},
					{1, math.SmallestNonzeroFloat64},
				},
				Beams: []hd.BeamProp{
					{N: [2]int{0, 1}, A: 1e200, J: 1e200, E: 1e200},
					{N: [2]int{1, 2}, A: 1e200, J: 1e200, E: 1e200},
				},
				Supports: [][3]bool{
					{true, true, true},
					{false, false, false},
					{false, false, false},
				},
			},
			lcs: []hd.LoadCase{
				{
					LoadNodes: []hd.LoadNode{
						{N: 2, Forces: [3]float64{1e20, 1e20, 0}},
					},
				},
			},
			mcs: []hd.ModalCase{
				{ModalMasses: []hd.ModalMass{
					{N: 1, Mass: 1e20},
					{N: -10, Mass: 1e20},
					{N: 10, Mass: 1e20},
				}},
			},
		},
		// error - two models in one model
		{
			m: hd.Model{
				Points: [][2]float64{
					{0, 0},
					{1, 0},
					{2, 0},
					{0, 1},
					{1, 1},
					{2, 0},
				},
				Beams: []hd.BeamProp{
					{N: [2]int{0, 1}, A: 1e200, J: 1e200, E: 1e200},
					{N: [2]int{1, 2}, A: 1e200, J: 1e200, E: 1e200},
					{N: [2]int{3, 4}, A: 1e200, J: 1e200, E: 1e200},
					{N: [2]int{4, 5}, A: 1e200, J: 1e200, E: 1e200},
				},
				Supports: [][3]bool{
					{true, true, true},
					{false, false, false},
					{false, false, false},
					{true, true, true},
					{false, false, false},
					{false, false, false},
				},
			},
			lcs: []hd.LoadCase{
				{
					LoadNodes: []hd.LoadNode{
						{N: 2, Forces: [3]float64{1e20, 1e20, 0}},
						{N: 5, Forces: [3]float64{1e20, 1e20, 0}},
					},
				},
			},
			mcs: []hd.ModalCase{
				{ModalMasses: []hd.ModalMass{
					{N: 1, Mass: 1},
				}},
			},
		},
	}
	for i := range ms {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			var err error
			defer func() {
				if err == nil {
					if r := recover(); r == nil {
						t.Errorf("panic is not happen and no error")
					} else {
						t.Log(r)
						debug.PrintStack()
					}
				}
			}()
			var b bytes.Buffer
			err = hd.Run(&b, &(ms[i].m), ms[i].lcs, ms[i].mcs)
			if err == nil {
				t.Fatalf("Error : %v", err)
			}
			b.Write([]byte(err.Error()))
			t.Log(err)
		})
	}
}

func TestCodeStyle(t *testing.T) {
	cs.All(t)
}

func TestWriter(t *testing.T) {
	m, lcs, mcs := example.Truss()
	var tf *os.File
	var err error
	if tf, err = ioutil.TempFile("", "testWriter"); err != nil {
		t.Fatalf("Cannot create temp file: %v", err)
	}
	if err = hd.Run(tf, &m, lcs, mcs); err != nil {
		t.Fatalf("Cannot calculate : %v", err)
	}
	if err = tf.Close(); err != nil {
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

func TestDirectionLoadNode(t *testing.T) {
	m := hd.Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{1.0, 0.0},
		},
		Beams: []hd.BeamProp{
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
	}
	lcs := []hd.LoadCase{
		{
			LoadNodes: []hd.LoadNode{
				{N: 1, Forces: [3]float64{100.0, 0.0, 0.0}},
			},
		},
		{
			LoadNodes: []hd.LoadNode{
				{N: 1, Forces: [3]float64{0.0, 100.0, 0.0}},
			},
		},
		{
			LoadNodes: []hd.LoadNode{
				{N: 1, Forces: [3]float64{0.0, 0.0, 100.0}},
			},
		},
	}

	var b bytes.Buffer
	if err := hd.Run(&b, &m, lcs, nil); err != nil {
		t.Fatalf("Cannot calculate : %v", err)
	}
	b.Reset()

	if !(lcs[0].PointDisplacementGlobal[1][0] > 0) {
		t.Errorf("load in direction X is not ok. See: %v", lcs[0].PointDisplacementGlobal[1])
	}
	if !(lcs[1].PointDisplacementGlobal[1][1] > 0) {
		t.Errorf("load in direction Y is not ok. See: %v", lcs[1].PointDisplacementGlobal[1])
	}
	if !(lcs[2].PointDisplacementGlobal[1][1] > 0) {
		t.Errorf("load in direction M is not ok. See: %v", lcs[2].PointDisplacementGlobal[1])
	}
	t.Logf("case 0 : %v", lcs[0].PointDisplacementGlobal)
	t.Logf("case 1 : %v", lcs[1].PointDisplacementGlobal)
	t.Logf("case 2 : %v", lcs[2].PointDisplacementGlobal)
}

// TODO: add test with rotation of beam with load

func Example() {

	m := hd.Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{1.0, 0.0},
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
		},
		Pins: [][6]bool{
			{false, false, false, false, false, true},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
		},
	}
	lcs := []hd.LoadCase{
		{
			LoadNodes: []hd.LoadNode{
				{N: 1, Forces: [3]float64{0, 2.3, 0}},
				{N: 1, Forces: [3]float64{-10, 0, 0}},
			},
			AmountLinearBuckling: 1,
		},
		{ // test for 2 cases with different positions
			LoadNodes: []hd.LoadNode{
				{N: 1, Forces: [3]float64{-10, 0, 0}},
				{N: 1, Forces: [3]float64{0, 2.3, 0}},
			},
			AmountLinearBuckling: 2,
		},
	}
	mcs := []hd.ModalCase{
		{
			ModalMasses: []hd.ModalMass{{N: 1, Mass: 10000}},
		},
	}

	var b bytes.Buffer
	if err := hd.Run(&b, &m, lcs, mcs); err != nil {
		panic(fmt.Errorf("Cannot calculate : %v", err))
	}

	expect, err := ioutil.ReadFile("./example/testdata/model.String")
	if err != nil {
		panic(fmt.Errorf("Cannot read file : %v", err))
	}

	if bytes.Equal(expect, b.Bytes()) {
		fmt.Fprintf(os.Stdout, "same")
	} else {
		fmt.Fprintln(os.Stdout, b.String())

		// show a diff between files
		diff := difflib.UnifiedDiff{
			A:        difflib.SplitLines(string(expect)),
			B:        difflib.SplitLines(string(b.Bytes())),
			FromFile: "Original",
			ToFile:   "Current",
			Context:  30000,
		}
		text, _ := difflib.GetUnifiedDiffString(diff)
		fmt.Fprintf(os.Stdout, text)

		fmt.Fprintf(os.Stdout, "not same")
	}

	// Output:
	// same
}

func ExampleKHM() {
	K := mat.NewDense(3, 3, []float64{1, 2, 3, 2, 5, 2, 3, 2, 1})
	if testing.Verbose() {
		fa := mat.Formatted(K, mat.Prefix("    "), mat.Squeeze())
		fmt.Fprintf(os.Stdout, "K = %.3g\n\n", fa)
	}

	M := mat.NewDense(3, 3, []float64{5, 0, 0, 0, 5, 0, 0, 0, 5})
	if testing.Verbose() {
		fa := mat.Formatted(M, mat.Prefix("    "), mat.Squeeze())
		fmt.Fprintf(os.Stdout, "M = %.3g\n\n", fa)
	}

	H, _ := hd.KHM(K, M)
	if testing.Verbose() {
		fa := mat.Formatted(H, mat.Prefix("    "), mat.Squeeze())
		fmt.Fprintf(os.Stdout, "H = %.3g\n\n", fa)
	}
	// Output:
	// K = ⎡1  2  3⎤
	//     ⎢2  5  2⎥
	//     ⎣3  2  1⎦
	//
	// M = ⎡5  0  0⎤
	//     ⎢0  5  0⎥
	//     ⎣0  0  5⎦
	//
	// H = ⎡-0.208  -0.833    2.29⎤
	//     ⎢-0.833    1.67  -0.833⎥
	//     ⎣  2.29  -0.833  -0.208⎦
}
