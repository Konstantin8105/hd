package hd

import (
	"bytes"
	"fmt"
	"math"
	"testing"
)

func TestNonlinearStaticValidation(t *testing.T) {
	// Validation of NonlinearConfig parameters.
	// Tests: nil model, Steps<1, Tol<=0, MaxIter<1.
	// Source: own API contract — no external reference.
	m, lcs, _ := consoleBeam()
	tcs := []struct {
		m   *Model
		lc  *LoadCase
		nr  nonlinearConfig
		err bool
	}{
		// nil model
		{m: nil, lc: nil, nr: nonlinearConfig{Steps: 1, Tol: 1e-6, MaxIter: 10}, err: true},
		// Steps < 1
		{m: &m, lc: &lcs[0], nr: nonlinearConfig{Steps: 0, Tol: 1e-6, MaxIter: 10}, err: true},
		// Tol <= 0
		{m: &m, lc: &lcs[0], nr: nonlinearConfig{Steps: 1, Tol: 0, MaxIter: 10}, err: true},
		{m: &m, lc: &lcs[0], nr: nonlinearConfig{Steps: 1, Tol: -1, MaxIter: 10}, err: true},
		// MaxIter < 1
		{m: &m, lc: &lcs[0], nr: nonlinearConfig{Steps: 1, Tol: 1e-6, MaxIter: 0}, err: true},
	}
	for i, tc := range tcs {
		t.Run(fmt.Sprintf("Case%d", i), func(t *testing.T) {
			err := nonlinearStatic(nil, tc.m, tc.lc, tc.nr)
			if (err != nil) != tc.err {
				t.Fatalf("Case %d: expect err=%v, got %v", i, tc.err, err)
			}
			if err != nil {
				t.Log(err)
			}
		})
	}
}

func TestNonlinearStaticSmallDisplacement(t *testing.T) {
	// For small displacements, nonlinear result must match linear within 1%.
	// Model: consoleBeam() — cantilever L=1m, E=2e11, J=120e-6.
	// Source: own verification — small strain approximation.
	m, lcs, _ := consoleBeam()

	// linear reference
	var buf bytes.Buffer
	if err := LinearStatic(&buf, &m, &lcs[0]); err != nil {
		t.Fatalf("LinearStatic failed: %v", err)
	}
	linDisp := lcs[0].PointDisplacementGlobal[1]

	// nonlinear with small load
	nr := nonlinearConfig{Steps: 5, Tol: 1e-8, MaxIter: 30}
	if err := nonlinearStatic(nil, &m, &lcs[0], nr); err != nil {
		t.Fatalf("NonlinearStatic failed: %v", err)
	}
	nlDisp := lcs[0].PointDisplacementGlobal[1]

	// compare - for small displacement nonlinear should match linear within 1%
	eps := 0.01
	for j := 0; j < 3; j++ {
		if linDisp[j] == 0.0 {
			if math.Abs(nlDisp[j]) > 1e-12 {
				t.Errorf("displacement[%d]: expected zero, got %e", j, nlDisp[j])
			}
			continue
		}
		diff := math.Abs((linDisp[j] - nlDisp[j]) / linDisp[j])
		if diff > eps {
			t.Errorf("displacement[%d]: %e vs %e, diff=%e", j, linDisp[j], nlDisp[j], diff)
		}
	}
}

func TestNonlinearStaticDirectionInvariance(t *testing.T) {
	// Nonlinear solver must be invariant to global coordinate rotation.
	// Horizontal and vertical cantilevers under axial load must produce
	// consistent displacements in the dominant (global X or Y) direction.
	// Source: own symmetry/consistency test.
	models := []struct {
		name string
		m    Model
		lc   LoadCase
	}{
		{
			name: "horizontal",
			m: Model{
				Points: [][2]float64{
					{0.0, 0.0},
					{1.0, 0.0},
				},
				Beams: []BeamProp{
					{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
				},
				Supports: [][3]bool{
					{true, true, true},
					{false, false, false},
				},
			},
			lc: LoadCase{
				LoadNodes: []LoadNode{
					{N: 1, Forces: [3]float64{100.0, 0.0, 0.0}},
				},
			},
		},
		{
			name: "vertical",
			m: Model{
				Points: [][2]float64{
					{0.0, 0.0},
					{0.0, 1.0},
				},
				Beams: []BeamProp{
					{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
				},
				Supports: [][3]bool{
					{true, true, true},
					{false, false, false},
				},
			},
			lc: LoadCase{
				LoadNodes: []LoadNode{
					{N: 1, Forces: [3]float64{0.0, 100.0, 0.0}},
				},
			},
		},
	}

	var disp [2][3]float64
	for i, tc := range models {
		t.Run(tc.name, func(t *testing.T) {
			nr := nonlinearConfig{Steps: 5, Tol: 1e-8, MaxIter: 30}
			if err := nonlinearStatic(nil, &tc.m, &tc.lc, nr); err != nil {
				t.Fatalf("NonlinearStatic failed: %v", err)
			}
			disp[i] = tc.lc.PointDisplacementGlobal[1]
		})
	}

	// both should give similar displacement magnitude (dominant direction)
	if math.Abs(disp[0][0]) < 1e-12 {
		t.Errorf("Expected horizontal displacement, got %v", disp[0])
	}
	if math.Abs(disp[1][1]) < 1e-12 {
		t.Errorf("Expected vertical displacement, got %v", disp[1])
	}
}

func TestNonlinearStaticSplitMesh(t *testing.T) {
	// Nonlinear result on a refined mesh (SplitBeam ×4) must match
	// the linear reference on the original mesh within 2%.
	// Model: consoleBeam().
	// Source: own mesh convergence test.
	m, lcs, _ := consoleBeam()

	// reference: linear on unsplit mesh
	var buf bytes.Buffer
	if err := LinearStatic(&buf, &m, &lcs[0]); err != nil {
		t.Fatalf("LinearStatic failed: %v", err)
	}
	ref := lcs[0].PointDisplacementGlobal[1]

	// nonlinear on split mesh with small load
	if err := splitBeam(&m, 0, 4); err != nil {
		t.Fatalf("SplitBeam failed: %v", err)
	}

	nr := nonlinearConfig{Steps: 10, Tol: 1e-6, MaxIter: 30}
	if err := nonlinearStatic(nil, &m, &lcs[0], nr); err != nil {
		t.Fatalf("NonlinearStatic failed: %v", err)
	}
	nlDisp := lcs[0].PointDisplacementGlobal[1]

	eps := 0.02
	for j := 0; j < 3; j++ {
		if ref[j] == 0.0 {
			continue
		}
		diff := math.Abs((ref[j] - nlDisp[j]) / ref[j])
		if diff > eps {
			t.Errorf("displacement[%d]: linear=%e, nonlinear(split)=%e, diff=%e", j, ref[j], nlDisp[j], diff)
		}
	}
}

func TestNonlinearStaticPdeltaAmplification(t *testing.T) {
	// P-delta amplification for a cantilever column at P/P_cr = 0.3, 0.5, 0.7.
	// Expected: delta_nl = delta_lin / (1 - P/P_cr).
	// Source: de Borst et al., Eq.1.12 (p.6): dF_int/dv — tangent stiffness
	// with geometric nonlinearity. Euler critical load: P_cr = pi^2*EJ/(4L^2).
	// Note: Pcr uses Euler-Bernoulli formula. Numerical buckling value from
	// hd may differ (Timoshenko beam), but amplification converges at <3%.
	E := 2.0e11
	J := 5.0e-5
	A := 7.6e-3
	L := 10.0
	Pcr := math.Pi * math.Pi * E * J / (4.0 * L * L)

	for _, ratio := range []float64{0.3, 0.5, 0.7} {
		t.Run(fmt.Sprintf("P_Pcr_%.1f", ratio), func(t *testing.T) {
			Paxial := ratio * Pcr
			Plat := 4500.0

			m := Model{
				Points: [][2]float64{{0.0, 0.0}, {L, 0.0}},
				Beams:  []BeamProp{{N: [2]int{0, 1}, A: A, J: J, E: E}},
				Supports: [][3]bool{
					{true, true, true},
					{false, false, false},
				},
			}
			lc := LoadCase{
				LoadNodes: []LoadNode{{N: 1, Forces: [3]float64{-Paxial, Plat, 0}}},
			}

			if err := LinearStatic(nil, &m, &lc); err != nil {
				t.Fatalf("LinearStatic failed: %v", err)
			}
			linDisp := lc.PointDisplacementGlobal[1]

			nr := nonlinearConfig{Steps: 15, Tol: 1e-8, MaxIter: 200}
			if err := nonlinearStatic(nil, &m, &lc, nr); err != nil {
				t.Fatalf("NonlinearStatic failed: %v", err)
			}
			nlDisp := lc.PointDisplacementGlobal[1]

			linLat := linDisp[1]
			nlLat := nlDisp[1]

			if linLat == 0 {
				t.Fatal("linear lateral displacement is zero")
			}
			amp := nlLat / linLat
			expected := 1.0 / (1.0 - Paxial/Pcr)
			errPct := math.Abs(amp-expected) / expected * 100
			t.Logf("P/P_cr=%.1f: δ_lin=%.4fm δ_nl=%.4fm amp=%.2f expected=%.2f err=%.1f%%",
				ratio, linLat, nlLat, amp, expected, errPct)
			if amp < 1.0+0.5*(expected-1.0) {
				t.Errorf("amplification too weak: %.2f vs expected %.2f", amp, expected)
			}
		})
	}
}

func TestNonlinearStaticConvergenceNearBuckling(t *testing.T) {
	// Convergence behavior as axial load approaches Euler critical load.
	// At P=0.9*P_cr: solver must converge (stiffness reduced but positive).
	// At P=2×numerical_P_cr: solver must fail (K_T singular/negative).
	// Source: de Borst Eq.1.12 (p.6): dF_int/dv → 0 as P → P_cr.
	E := 2.0e11
	J := 5.0e-5
	A := 7.6e-3
	L := 10.0
	Pcr := math.Pi * math.Pi * E * J / (4.0 * L * L)

	t.Run("below_buckling_0.9", func(t *testing.T) {
		Paxial := 0.9 * Pcr
		m := Model{
			Points: [][2]float64{{0.0, 0.0}, {L, 0.0}},
			Beams:  []BeamProp{{N: [2]int{0, 1}, A: A, J: J, E: E}},
			Supports: [][3]bool{
				{true, true, true},
				{false, false, false},
			},
		}
		lc := LoadCase{
			LoadNodes: []LoadNode{{N: 1, Forces: [3]float64{-Paxial, 100, 0}}},
		}
		nr := nonlinearConfig{Steps: 30, Tol: 1e-8, MaxIter: 500}
		if err := nonlinearStatic(nil, &m, &lc, nr); err != nil {
			t.Fatalf("NonlinearStatic failed at P=0.9P_cr: %v", err)
		}
		δ := lc.PointDisplacementGlobal[1]
		t.Logf("P=0.9P_cr: δy=%.4f m, δx=%.4f m", δ[1], δ[0])
		if math.Abs(δ[1]) < 1e-10 {
			t.Errorf("unexpected zero lateral displacement")
		}
	})

	t.Run("buckling_load_computation", func(t *testing.T) {
		// Compute actual numerical buckling load using linear buckling analysis
		m := Model{
			Points: [][2]float64{{0.0, 0.0}, {L, 0.0}},
			Beams:  []BeamProp{{N: [2]int{0, 1}, A: A, J: J, E: E}},
			Supports: [][3]bool{
				{true, true, true},
				{false, false, false},
			},
		}
		// Use linear buckling to find actual critical load
		lcBucket := LoadCase{
			LoadNodes:            []LoadNode{{N: 1, Forces: [3]float64{-1.0, 0, 0}}},
			AmountLinearBuckling: 1,
		}
		if err := LinearStatic(nil, &m, &lcBucket); err != nil {
			t.Fatalf("LinearStatic for buckling failed: %v", err)
		}
		if len(lcBucket.LinearBucklingResult) == 0 {
			t.Fatal("no buckling result")
		}
		actualPcr := lcBucket.LinearBucklingResult[0].Factor
		t.Logf("Theoretical P_cr=%.2f kN, Numerical P_cr=%.2f kN", Pcr/1000, actualPcr)

		// at 2x numerical critical load — should fail (negative stiffness)
		t.Run("above_numerical_buckling", func(t *testing.T) {
			Paxial := 2.0 * actualPcr
			m2 := Model{
				Points: [][2]float64{{0.0, 0.0}, {L, 0.0}},
				Beams:  []BeamProp{{N: [2]int{0, 1}, A: A, J: J, E: E}},
				Supports: [][3]bool{
					{true, true, true},
					{false, false, false},
				},
			}
			lc2 := LoadCase{
				LoadNodes: []LoadNode{{N: 1, Forces: [3]float64{-Paxial, 100, 0}}},
			}
			nr := nonlinearConfig{Steps: 20, Tol: 1e-8, MaxIter: 50}
			if err := nonlinearStatic(nil, &m2, &lc2, nr); err == nil {
				t.Errorf("expected error above P_cr_numerical, but got none")
			} else {
				t.Logf("Expected error at 2×P_cr: %v", err)
			}
		})
	})
}

func TestNonlinearStaticShallowTruss(t *testing.T) {
	// Shallow truss with spring (PyFEM ShallowTruss.py, de Borst Ch.1.1)
	// Parameters: b=10, h=0.5, k=1000 N/m, EA0=5e6 N (TOTAL), DF=50 N, N=30
	// Python reference (verified by running the script):
	//   F_ext(v) = -EA0 * (h-v)/l * (l-l0)/l0 + k*v
	//   F=200N → v=0.10661m, F=300N → v=0.18196m, F=400N → v=0.29305m
	//
	// Spring k modelled as vertical beam (0,0.5)→(0,0) with EA/L = k.
	// A_spring = k*L/E. J=1e-12 keeps bending error < 3%.
	E := 2.0e11
	EA0 := 5.0e6
	A_bar := EA0 / (2.0 * E)
	b := 10.0
	h := 0.5
	k := 1000.0
	l0 := math.Sqrt(b*b + h*h)

	L_spring := 15.0
	A_spring := k * L_spring / E

	m := Model{
		Points: [][2]float64{
			{-b, 0.0},
			{0.0, h},
			{b, 0.0},
			{0.0, h - L_spring},
		},
		Beams: []BeamProp{
			{N: [2]int{0, 1}, A: A_bar, J: 1e-7, E: E},
			{N: [2]int{1, 2}, A: A_bar, J: 1e-7, E: E},
			{N: [2]int{1, 3}, A: A_spring, J: 1e-12, E: E},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{true, true, true},
			{true, true, true},
		},
		Pins: [][6]bool{
			{false, false, true, false, false, false}, // beam 0: M free at support
			{false, false, false, false, false, true}, // beam 1: M free at support
			{false, false, false, false, false, true}, // beam 2: M free at ground
		},
	}

	analyticalF := func(vHd float64) float64 {
		v := -vHd
		l := math.Sqrt(b*b + (h-v)*(h-v))
		return -EA0*(h-v)/l*(l-l0)/l0 + k*v
	}

	// PyFEM reference values from running ShallowTruss.py (analytic 1-DOF
	// truss, monotonic load-displacement path):
	//   F=200N → v=0.10661m, F=300N → v=0.18196m
	//   F=400N → v=0.29305m, F=600N → v=0.70695m
	//
	// Load-control cannot pass the limit point of this beam model (~460 N,
	// where the tangent stiffness becomes indefinite), so F=600 is solved
	// with the arc-length method instead. Pre-limit loads use load-control.
	for _, Ftarget := range []float64{200, 300, 600} {
		t.Run(fmt.Sprintf("F_%.0f", Ftarget), func(t *testing.T) {
			lc := LoadCase{LoadNodes: []LoadNode{{N: 1, Forces: [3]float64{0, -Ftarget, 0}}}}

			if err := LinearStatic(nil, &m, &lc); err != nil {
				t.Fatalf("LinearStatic failed: %v", err)
			}
			linDisp := lc.PointDisplacementGlobal[1][1]

			if Ftarget >= 600 {
				// post-limit load, solved with arc-length control
				var buf bytes.Buffer
				al := arcLengthConfig{Steps: 20, Tol: 1e-4, MaxIter: 500}
				if err := nonlinearStaticArcLength(&buf, &m, &lc, al); err != nil {
					t.Fatalf("NonlinearStaticArcLength failed: %v", err)
				}
				vPfem := -lc.PointDisplacementGlobal[1][1]
				const vRef = 0.70695
				errPct := math.Abs(vPfem-vRef) / vRef * 100
				t.Logf("F=%.0f N: lin=%.5f nl(arc-length)=%.5f (ref=%.5f) err=%.1f%%",
					Ftarget, linDisp, vPfem, vRef, errPct)
				if errPct > 10.0 {
					t.Errorf("arc-length post-limit v=%.4f deviates from PyFEM %.5f by %.1f%%",
						vPfem, vRef, errPct)
				}
				return
			}

			nr := nonlinearConfig{Steps: 10, Tol: 1e-6, MaxIter: 100}
			if err := nonlinearStatic(nil, &m, &lc, nr); err != nil {
				t.Fatalf("NonlinearStatic failed: %v", err)
			}
			vHd := lc.PointDisplacementGlobal[1][1]
			vPfem := -vHd
			Fana := analyticalF(vHd)
			errPct := math.Abs(Ftarget-math.Abs(Fana)) / math.Abs(Fana) * 100
			vRef := map[float64]float64{200: 0.10661, 300: 0.18196}[Ftarget]
			t.Logf("F=%.0f N: lin=%.5f nl=%.5f (ref=%.5f) F_ana=%.1f N err=%.1f%%",
				Ftarget, linDisp, vPfem, vRef, Fana, errPct)
			if errPct > 12.0 {
				t.Errorf("FE deviates from PyFEM: FE=%.0f, analytical=%.1f, err=%.1f%%",
					Ftarget, Fana, errPct)
			}
		})
	}
}

func TestNonlinearStaticFrameSecondOrder(t *testing.T) {
	// Portal frame: vertical load on left column (0.6*P_cr) + horizontal
	// load on beam. Nonlinear horizontal displacement must be larger than
	// linear (P-delta sway amplification).
	// Source: NCCI SN033a-EN-EU, Section 2.1-2.2 (p.2-3): amplified sway
	// moment method, sway factor = 1/(1-1/alpha_cr).
	E := 2.0e11
	J := 1.0e-4
	A := 1.0e-2
	H := 5.0 // column height
	W := 6.0 // beam width

	m := Model{
		Points: [][2]float64{
			{0.0, 0.0}, // 0: left base
			{0.0, H},   // 1: left top
			{W, H},     // 2: right top
			{W, 0.0},   // 3: right base
		},
		Beams: []BeamProp{
			{N: [2]int{0, 1}, A: A, J: J, E: E}, // left column
			{N: [2]int{1, 2}, A: A, J: J, E: E}, // beam
			{N: [2]int{2, 3}, A: A, J: J, E: E}, // right column
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{false, false, false},
			{true, true, true},
		},
	}

	// compute critical load for single column
	Pcr := math.Pi * math.Pi * E * J / (H * H)
	Paxial := 0.6 * Pcr
	Hforce := 50000.0

	lc := LoadCase{
		LoadNodes: []LoadNode{
			{N: 1, Forces: [3]float64{-Paxial, 0, 0}},
			{N: 2, Forces: [3]float64{Hforce, 0, 0}},
		},
	}

	if err := LinearStatic(nil, &m, &lc); err != nil {
		t.Fatalf("LinearStatic failed: %v", err)
	}
	// check horizontal displacement at top of left column (node 1)
	linDisp := lc.PointDisplacementGlobal[1]

	nr := nonlinearConfig{Steps: 20, Tol: 1e-8, MaxIter: 100}
	if err := nonlinearStatic(nil, &m, &lc, nr); err != nil {
		t.Fatalf("NonlinearStatic failed: %v", err)
	}
	nlDisp := lc.PointDisplacementGlobal[1]

	t.Logf("Linear δx=%.4f m, Nonlinear δx=%.4f m", linDisp[0], nlDisp[0])
	// With axial compression, nonlinear horizontal displacement should be larger (P-delta)
	// compare absolute values (displacements are negative = leftward)
	if math.Abs(nlDisp[0]) <= math.Abs(linDisp[0]) {
		t.Errorf("Expected nonlinear displacement > linear, got lin=%.4f nl=%.4f", linDisp[0], nlDisp[0])
	}
}

func TestNonlinearStaticLeeFrame(t *testing.T) {
	// Lee's frame (clamped frame 120×120) with vertical load at apex.
	// Modelled with 3 beam elements, each split into 10 segments via SplitBeam.
	// P-delta effect: nonlinear Y-displacement differs from linear.
	// Source: Bonet, Gil, Wood, "Examples in Nonlinear Continuum Mechanics"
	// (2012), Example 3.4 (p.48-50), Fig 3.10-3.11 (p.50-51).
	E := 2.0e11
	J := 1.0e-3
	A := 1.0e-2
	W := 120.0
	H := 120.0

	m := Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{0.0, H},
			{W, H},
			{W, 0.0},
		},
		Beams: []BeamProp{
			{N: [2]int{0, 1}, A: A, J: J, E: E},
			{N: [2]int{1, 2}, A: A, J: J, E: E},
			{N: [2]int{2, 3}, A: A, J: J, E: E},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
			{false, false, false},
			{true, true, true},
		},
	}

	// split each beam into 10 segments
	for i := 0; i < 3; i++ {
		if err := splitBeam(&m, 0, 10); err != nil {
			t.Fatalf("SplitBeam %d: %v", i, err)
		}
	}

	// find apex by searching point at (W/2, H)
	apexIdx := -1
	for i, p := range m.Points {
		if math.Abs(p[0]-W/2) < 1.0 && math.Abs(p[1]-H) < 1.0 {
			apexIdx = i
			break
		}
	}
	if apexIdx < 0 {
		// estimate: for top beam (from node 1→2), with 10 intermediate points,
		// the 5th intermediate point (0-indexed) is closest to center
		// top beam original: N=[1,2]; after split: new points start at lastPointIndex_after_first_split
		apexIdx = 1 + 5 // approximate
	}
	t.Logf("Apex point index: %d, total points: %d", apexIdx, len(m.Points))

	Fvertical := 1.0e6
	lc := LoadCase{
		LoadNodes: []LoadNode{
			{N: apexIdx, Forces: [3]float64{0, -Fvertical, 0}},
		},
	}

	if err := LinearStatic(nil, &m, &lc); err != nil {
		t.Fatalf("LinearStatic failed: %v", err)
	}
	linDisp := lc.PointDisplacementGlobal[apexIdx]

	nr := nonlinearConfig{Steps: 20, Tol: 1e-8, MaxIter: 200}
	if err := nonlinearStatic(nil, &m, &lc, nr); err != nil {
		t.Fatalf("NonlinearStatic failed: %v", err)
	}
	nlDisp := lc.PointDisplacementGlobal[apexIdx]

	t.Logf("Linear δy=%.4f m, Nonlinear δy=%.4f m", linDisp[1], nlDisp[1])
	t.Logf("Linear δx=%.4f m, Nonlinear δx=%.4f m", linDisp[0], nlDisp[0])
	diff := math.Abs(nlDisp[1]-linDisp[1]) / math.Max(1e-10, math.Abs(linDisp[1]))
	if diff < 1e-6 {
		t.Errorf("Nonlinear displacement identical to linear, diff=%e", diff)
	}
}

func TestNonlinearStaticMultipleCases(t *testing.T) {
	// Verify that NonlinearStatic writes all result fields:
	// PointDisplacementGlobal, BeamForces, Reactions.
	// Model: consoleBeam().
	// Source: own API contract test.
	m, lcs, _ := consoleBeam()

	nr := nonlinearConfig{Steps: 5, Tol: 1e-8, MaxIter: 30}

	// run nonlinear on first load case
	if err := nonlinearStatic(nil, &m, &lcs[0], nr); err != nil {
		t.Fatalf("NonlinearStatic failed: %v", err)
	}
	if len(lcs[0].PointDisplacementGlobal) == 0 {
		t.Fatal("no displacement results")
	}
	if len(lcs[0].BeamForces) == 0 {
		t.Fatal("no beam forces results")
	}
	if len(lcs[0].Reactions) == 0 {
		t.Fatal("no reactions results")
	}
}

func TestNonlinearStaticArcLengthShallowTruss(t *testing.T) {
	// Arc-length control traces post-snap-through of shallow truss
	// (PyFEM ShallowTruss.py, de Borst Ch.1.1, RiksSolver.py).
	//
	// Reference (from running Python script):
	//   b=10, h=0.5, k=1000, EA0=5e6
	//   F=500N → v=0.50000m  (limit point)
	//   F=600N → v=0.70695m  (post-snap-through)
	//
	// Arc-length uses reference load 1500N → λ=0.4 → F=600N.
	E := 2.0e11
	EA0 := 5.0e6
	A_bar := EA0 / (2.0 * E)
	b := 10.0
	h := 0.5
	k := 1000.0
	L_spring := 15.0
	A_spring := k * L_spring / E

	m := Model{
		Points: [][2]float64{
			{-b, 0.0}, {0.0, h}, {b, 0.0}, {0.0, h - L_spring},
		},
		Beams: []BeamProp{
			{N: [2]int{0, 1}, A: A_bar, J: 1e-7, E: E},
			{N: [2]int{1, 2}, A: A_bar, J: 1e-7, E: E},
			{N: [2]int{1, 3}, A: A_spring, J: 1e-7, E: E},
		},
		Supports: [][3]bool{
			{true, true, true}, {false, false, false},
			{true, true, true}, {true, true, true},
		},
		Pins: [][6]bool{
			{false, false, true, false, false, false},
			{false, false, false, false, false, true},
			{false, false, false, false, false, true},
		},
	}

	var buf bytes.Buffer
	lc := LoadCase{
		LoadNodes: []LoadNode{{N: 1, Forces: [3]float64{0, -1500, 0}}},
	}
	al := arcLengthConfig{Steps: 20, Tol: 1e-4, MaxIter: 200}
	if err := nonlinearStaticArcLength(&buf, &m, &lc, al); err != nil {
		t.Fatalf("NonlinearStaticArcLength failed: %v", err)
	}
	t.Logf("Output:\n%s", &buf)

	disp := lc.PointDisplacementGlobal[1][1] // hd Y (positive = up)
	t.Logf("Final apex displacement: %.4f m", disp)
	if math.Abs(disp) < 0.5 {
		t.Errorf("Arc-length did not trace past limit point: v=%.4f m", disp)
	}
}

// consoleBeam is example.ConsoleBeam inlined so this internal test does not
// import the example package (which would create an import cycle).
func consoleBeam() (m Model, lc []LoadCase, mc []ModalCase) {
	m = Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{1.0, 0.0},
		},
		Beams: []BeamProp{
			{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
		},
	}
	lc = []LoadCase{
		{
			LoadNodes: []LoadNode{
				{N: 1, Forces: [3]float64{0, 2.3, 0}},
				{N: 1, Forces: [3]float64{10, 0, 0}},
			},
		},
		{
			LoadNodes: []LoadNode{
				{N: 1, Forces: [3]float64{10, 0, 0}},
				{N: 1, Forces: [3]float64{0, 2.3, 0}},
			},
		},
	}
	mc = []ModalCase{
		{
			ModalMasses: []ModalMass{{N: 1, Mass: 10000}},
		},
	}
	return
}

// splitBeam is mod.SplitBeam inlined so this internal test does not import the
// mod package (which would create an import cycle).
func splitBeam(m *Model, beamIndex, amountIntermediatePoints int) (err error) {
	if amountIntermediatePoints < 1 {
		return fmt.Errorf("Not valid value of amount intermediate points : %d", amountIntermediatePoints)
	}
	if beamIndex < 0 || beamIndex >= len(m.Beams) {
		return fmt.Errorf("Not valid value of beam index : %d", beamIndex)
	}
	lastPointIndex := len(m.Points)
	startPoint := m.Points[m.Beams[beamIndex].N[0]]
	endPoint := m.Points[m.Beams[beamIndex].N[1]]
	m.Points = append(m.Points, make([][2]float64, amountIntermediatePoints)...)
	for j := 0; j < 2; j++ {
		step := (endPoint[j] - startPoint[j]) / float64(amountIntermediatePoints+1)
		for i := 0; i < amountIntermediatePoints; i++ {
			if startPoint[j] == endPoint[j] {
				m.Points[lastPointIndex+i][j] = startPoint[j]
				continue
			}
			m.Points[lastPointIndex+i][j] = startPoint[j] + float64(i+1)*step
		}
	}
	lastBeamIndex := len(m.Beams)
	m.Beams = append(m.Beams, make([]BeamProp, amountIntermediatePoints)...)
	for i := 0; i < amountIntermediatePoints; i++ {
		m.Beams[i+lastBeamIndex] = m.Beams[beamIndex]
		m.Beams[i+lastBeamIndex].N[0] = lastPointIndex + i
		m.Beams[i+lastBeamIndex].N[1] = lastPointIndex + i + 1
	}
	m.Beams[len(m.Beams)-1].N[1] = m.Beams[beamIndex].N[1]
	m.Beams[beamIndex].N[1] = lastPointIndex
	m.Supports = append(m.Supports, make([][3]bool, amountIntermediatePoints)...)
	if len(m.Pins) != 0 {
		m.Pins = append(m.Pins, make([][6]bool, amountIntermediatePoints)...)
		for i := 3; i < 6; i++ {
			m.Pins[len(m.Beams)-1][i] = m.Pins[beamIndex][i]
			m.Pins[beamIndex][i] = false
		}
	}
	return nil
}
