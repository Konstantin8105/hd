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
					LinearBuckling: struct {
						Amount  uint16
						Results []hd.BucklingResult
					}{
						Amount: 1,
					},
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
			F = 100.0
			L = 2.200
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
			// section 10B1 STO ASCH
			{N: [2]int{0, 1}, A: 10.32e-4, J: 171e-8, E: 2.05e11},
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
				{N: 1, Forces: [3]float64{0, 2.3e3, 0}},
				{N: 1, Forces: [3]float64{-10e3, 0, 0}},
			},
			LinearBuckling: struct {
				Amount  uint16
				Results []hd.BucklingResult
			}{
				Amount: 1,
			},
		},
		{ // test for 2 cases with different positions
			LoadNodes: []hd.LoadNode{
				{N: 1, Forces: [3]float64{-10e3, 0, 0}},
				{N: 1, Forces: [3]float64{0, 2.3e3, 0}},
			},
			LinearBuckling: struct {
				Amount  uint16
				Results []hd.BucklingResult
			}{
				Amount: 2,
			},
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

	filename := "./example/testdata/model.String"

	if os.Getenv("UPDATE") != "" {
		err := ioutil.WriteFile(filename, b.Bytes(), 0644)
		if err != nil {
			panic(fmt.Errorf("Cannot Update: %v", err))
		}
	}

	expect, err := ioutil.ReadFile(filename)
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
			B:        difflib.SplitLines(b.String()),
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

func ExampleNonlinear() {

	// 	d := 0.14
	// 	E := 2.05e11
	// 	J := math.Pow(d, 4) / 12
	// 	A := d * d
	m := hd.Model{
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
	lc := hd.LoadCase{
		LoadNodes: []hd.LoadNode{
			{N: 1, Forces: [3]float64{a * P, -P, 0}},
		},
	}
	// 	if isLinear {
	//		lc.LinearBuckling.Amount = 1
	// 	} else {
	lc.NonlinearNK.MaxIterations = 29000
	lc.NonlinearNK.Substep = 40
	// 	}

	// fmt.Println(m)

	// var buf bytes.Buffer
	err := hd.LinearStatic(os.Stdout, &m, &lc)
	if err != nil {
		// panic(err)
		fmt.Println(">>>>>>>>", err)
		return
	}

	// fmt.Println(lc)
	fmt.Println("> results :")
	for i := range lc.NonlinearNK.Results {
		fmt.Printf("%12.5f, %10.3f\n",
			lc.NonlinearNK.Results[i].PointDisplacementGlobal[1][0]*1000,
			lc.NonlinearNK.Results[i].Reactions[2][1]/1000)
	}

	// 	lc = append([]hd.LoadCase{}, l)
	// 	return

	// 	// Test based on example from document:
	// 	// Nonlinear Analysis of Structures
	// 	// The Arc Length Method: Formulation, Implementation and Applications
	// 	// Nikolaos Vasios
	//
	// 	var (
	// 		EA = 1.0
	// 		E  = 2.0e11
	// 		Ao = EA / E
	// 		Jo = 1.0e-8
	//
	// 		Lo = 10.0
	// 		dy = 5.0
	// 		dx = math.Sqrt(Lo*Lo - dy*dy)
	// 		θo = math.Atan(dy / dx)
	//
	// 		C  = 0.02 // EAL
	// 		Lc = 10.0
	// 		Ac = C * Lc / E
	// 	)
	//
	// 	m := hd.Model{
	// 		Points: [][2]float64{
	// 			{0.0, 0.0},                       // 0  Support
	// 			{dx * 1.0 / 6.0, dy * 1.0 / 6.0}, // 1
	// 			{dx * 2.0 / 6.0, dy * 2.0 / 6.0}, // 2
	// 			{dx * 3.0 / 6.0, dy * 3.0 / 6.0}, // 3
	// 			{dx * 4.0 / 6.0, dy * 4.0 / 6.0}, // 4
	// 			{dx * 5.0 / 6.0, dy * 5.0 / 6.0}, // 5
	// 			{dx, dy},                         // 6 Load point
	// 			{dx, dy - Lc},
	// 		},
	// 		Beams: []hd.BeamProp{
	// 			{N: [2]int{0, 1}, A: Ao, J: Jo, E: E},
	// 			{N: [2]int{1, 2}, A: Ao, J: Jo, E: E},
	// 			{N: [2]int{2, 3}, A: Ao, J: Jo, E: E},
	// 			{N: [2]int{3, 4}, A: Ao, J: Jo, E: E},
	// 			{N: [2]int{4, 5}, A: Ao, J: Jo, E: E},
	// 			{N: [2]int{5, 6}, A: Ao, J: Jo, E: E},
	// 			{N: [2]int{6, 7}, A: Ac, J: Jo, E: E},
	// 		},
	// 		// Pins: [][6]bool{
	// 		// 	{false, true, true, false, false, false},
	// 		// 	{false, false, false, false, false, false},
	// 		// 	{false, false, false, false, false, false},
	// 		// 	{false, false, false, false, false, false},
	// 		// 	{false, false, false, false, false, false},
	// 		// 	{false, false, false, false, true, true},
	// 		// 	{false, true, true, false, true, true},
	// 		// },
	// 		Supports: [][3]bool{
	// 			{true, true, false},
	// 			{false, false, false},
	// 			{false, false, false},
	// 			{false, false, false},
	// 			{false, false, false},
	// 			{false, false, false},
	// 			{true, false, false},
	// 			{true, false, false},
	// 		},
	// 	}
	//
	// 	P := 1.4
	//
	// 	lc := hd.LoadCase{
	// 		LoadNodes: []hd.LoadNode{
	// 			{N: 7, Forces: [3]float64{0, -P / 2.0, 0}},
	// 		},
	// 		NonlinearNR: struct {
	// 			MaxIterations uint64
	// 			Substep       uint64
	// 			Results       []*hd.LoadCase
	// 		}{
	// 			MaxIterations: 10000,
	// 			Substep:       10,
	// 		},
	// 	}
	//
	// 	var buf bytes.Buffer
	// 	err := hd.LinearStatic(&buf, &m, &lc)
	// 	t.Log(err)
	//
	// 	k := E * Ao / Lo
	// 	t.Logf("dx = %v", dx)
	// 	t.Logf("dy = %v", dy)
	// 	t.Logf("Lo = %v", Lo)
	// 	t.Logf("k  = %v", k)
	// 	t.Logf("θo = %v", θo)
	// 	betta := E * Ao / Lo
	// 	t.Logf("betta = E*A/L   = %.5f", betta)
	// 	t.Logf("w     = betta/k = %.5f", betta/k)
	//
	// 	w := -lc.PointDisplacementGlobal[6][1]
	// 	t.Logf("w    = %.4f displacement of %.4f", w, dy)
	//
	// 	λ := math.Abs(lc.Reactions[0][1] * 2 / (2.0 * k * Lo))
	// 	a := math.Abs((-lc.PointDisplacementGlobal[6][1]) / Lo)
	// 	L2 := (1.0/math.Sqrt(1-2*a*math.Sin(θo)+a*a) - 1) * (math.Sin(θo) - a)
	// 	t.Logf("λ    = %.5f a = %.5f L2 = %.5f", λ, a, L2)
	//
	// 	// t.Logf("%s", m)
	// 	// t.Logf("%s", lc)
	// 	t.Logf("Force N\tDisplacement m")
	// 	for _, v := range lc.NonlinearNR.Results {
	// 		t.Logf(" %.4f\t%.4f\t%.f",
	// 			v.LoadNodes[0].Forces[1],
	// 			-v.PointDisplacementGlobal[6][1],
	// 			-v.PointDisplacementGlobal[7][1],
	// 		)
	// 	}
	//
	// 	// Z := dy - w
	// 	// t.Logf("P by research = %.4e", E*Ao*Z*(2*Z*w-w*w)/(2*Lo*Lo*Lo))
	// 	// t.Logf("z             = %.4e %s", Z, " displacament")
	// 	// t.Logf("P             = %.4e  %.4e", P, 2*k*Lo*L2)
	// 	// t.Logf("L/Lo          = %f %f", math.Sqrt(1.0-2*w/Lo*math.Sin(θo)+math.Pow(w/Lo, 2.0)), ( math.Sqrt(dx*dx+Z*Z) )/Lo)
	//
	// 	// ArcLength(E, Ao, Lo, θo)
	//
	// 	// }
	// 	// t.Log(buf.String())

	// Output:
}
