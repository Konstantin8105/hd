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

	if err := hd.Run(nil, &m, lcs, nil); err != nil {
		t.Fatalf("Cannot calculate : %v", err)
	}

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

func TestRotateBeamLinear(t *testing.T) {
	for angle := 0.0; angle <= 360.0; angle += 5 {
		const (
			F float64 = 100.0
			L         = 2.200
		)
		var (
			Fx = F * math.Cos(math.Pi/180.0*angle)
			Fy = F * math.Sin(math.Pi/180.0*angle)
			x  = L * math.Cos(math.Pi/180.0*angle)
			y  = L * math.Sin(math.Pi/180.0*angle)
		)
		m := hd.Model{
			Points: [][2]float64{{0.0, 0.0}, {x / 2, y / 2}, {x, y}},
			Beams: []hd.BeamProp{
				{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
				{N: [2]int{1, 2}, A: 12e-4, J: 120e-6, E: 2.0e11},
			},
			Supports: [][3]bool{{true, true, true}, {false, false, false}, {false, false, false}},
		}
		lcs := hd.LoadCase{LoadNodes: []hd.LoadNode{{N: 2, Forces: [3]float64{Fx, Fy, 0.0}}}}
		if err := hd.LinearStatic(nil, &m, &lcs); err != nil {
			t.Fatalf("Cannot calculate : %v", err)
		}
		// tolerance
		tol := 1e-6

		// compare beam force
		{
			react := math.Sqrt(math.Pow(lcs.BeamForces[0][0], 2) + math.Pow(lcs.BeamForces[0][1], 2))
			if diff := math.Abs((react - F) / F); diff > tol {
				t.Errorf("Not valid beam force on angle %5.1f : %8.5e", angle, diff)
			}
		}
		{
			react := math.Sqrt(math.Pow(lcs.BeamForces[0][3], 2) + math.Pow(lcs.BeamForces[0][4], 2))
			if diff := math.Abs((react - F) / F); diff > tol {
				t.Errorf("Not valid beam force on angle %5.1f : %8.5e", angle, diff)
			}
		}
		{
			react := math.Sqrt(math.Pow(lcs.BeamForces[1][0], 2) + math.Pow(lcs.BeamForces[1][1], 2))
			if diff := math.Abs((react - F) / F); diff > tol {
				t.Errorf("Not valid beam force on angle %5.1f : %8.5e", angle, diff)
			}
		}
		{
			react := math.Sqrt(math.Pow(lcs.BeamForces[1][3], 2) + math.Pow(lcs.BeamForces[1][4], 2))
			if diff := math.Abs((react - F) / F); diff > tol {
				t.Errorf("Not valid beam force on angle %5.1f : %8.5e", angle, diff)
			}
		}
		// compare reaction
		{
			react := math.Sqrt(math.Pow(lcs.Reactions[0][0], 2) + math.Pow(lcs.Reactions[0][1], 2))
			if diff := math.Abs((react - F) / F); diff > tol {
				t.Errorf("Not valid reaction on angle %5.1f : %8.5e", angle, diff)
			}
		}
	}
}

func TestRotateBeamModal(t *testing.T) {
	var hz [2]float64
	for angle := 0.0; angle <= 360.0; angle += 5 {
		const (
			L = 2.200
		)
		var (
			x = L * math.Cos(math.Pi/180.0*angle)
			y = L * math.Sin(math.Pi/180.0*angle)
		)
		m := hd.Model{
			Points: [][2]float64{{0.0, 0.0}, {x / 2, y / 2}, {x, y}},
			Beams: []hd.BeamProp{
				{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
				{N: [2]int{1, 2}, A: 12e-4, J: 120e-6, E: 2.0e11},
			},
			Supports: [][3]bool{{true, true, true}, {false, false, false}, {false, false, false}},
		}
		mass := hd.ModalCase{ModalMasses: []hd.ModalMass{{N: 1, Mass: 50.0}, {N: 2, Mass: 100.0}}}
		if err := hd.Modal(nil, &m, &mass); err != nil {
			t.Fatalf("Cannot calculate : %v", err)
		}
		if angle == 0.0 {
			hz[0] = mass.Result[0].Hz
			hz[1] = mass.Result[1].Hz
			continue
		}
		// tolerance
		tol := 1e-6

		// compare
		for i := 0; i < 2; i++ {
			if diff := math.Abs((hz[i] - mass.Result[i].Hz) / hz[i]); diff > tol {
				t.Errorf("Not valid frequency on angle %5.1f: %8.5e", angle, diff)
			}
		}
	}
}

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
		fmt.Fprintf(os.Stdout, "Cannot calculate : %v", err)
		return
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
