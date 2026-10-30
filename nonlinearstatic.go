package hd

import (
	"fmt"
	"io"
	"math"
	"runtime/debug"

	"github.com/Konstantin8105/errors"
	"github.com/Konstantin8105/sparse"
	"gonum.org/v1/gonum/mat"
)

// nonlinearConfig contains iteration parameters.
type nonlinearConfig struct {
	Steps   int     // number of load steps (>= 1)
	Tol     float64 // convergence tolerance (relative displacement change)
	MaxIter int     // max iterations per step
}

// nonlinearStatic performs geometrically nonlinear static analysis
// using fixed-point iteration on the tangent stiffness.
//
// At each load step:
//
//	(K_linear - K_geometric(N(u_k))) * u_{k+1} = p
//
// K_geometric is subtracted from K_linear (compression reduces stiffness).
// Convergence is checked on the relative change in displacements.
//
// Note: lc.PointDisplacementGlobal, lc.BeamForces, lc.Reactions are
// overwritten with results.
func nonlinearStatic(out io.Writer, m *Model, lc *LoadCase, nr nonlinearConfig) (err error) {
	if out == nil {
		out = io.Discard
	}

	name := "Nonlinear Static Analysis (Geometric nonlinearity)"
	fmt.Fprintf(out, "%s\n", name)

	defer func() {
		if err != nil {
			et := errors.New(name)
			_ = et.Add(err)
			err = et
		}
	}()

	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("stacktrace from panic: %s\n", debug.Stack())
		}
	}()

	if nr.Steps < 1 {
		return fmt.Errorf("nonlinearConfig.Steps must be >= 1, got %d", nr.Steps)
	}
	if nr.Tol <= 0 {
		return fmt.Errorf("nonlinearConfig.Tol must be > 0, got %e", nr.Tol)
	}
	if nr.MaxIter < 1 {
		return fmt.Errorf("nonlinearConfig.MaxIter must be >= 1, got %d", nr.MaxIter)
	}

	out, err = prepare(out, m)
	if err != nil {
		return fmt.Errorf("prepare model: %v", err)
	}
	if err := lc.checkInputData(m); err != nil {
		return fmt.Errorf("check load case: %v", err)
	}

	dof := 3 * len(m.Points)

	K_linear, ignore, err := m.assemblyK(func(pos int) *mat.Dense {
		return m.getStiffBeam2d(pos)
	})
	if err != nil {
		return fmt.Errorf("assembly linear stiffness: %v", err)
	}

	p_total, err := m.assemblyNodeLoad(lc)
	if err != nil {
		return fmt.Errorf("assembly node load: %v", err)
	}

	for _, i := range ignore {
		if p_total[i] != 0.0 {
			return fmt.Errorf("on direction %d load is not zero: %f", i, p_total[i])
		}
	}

	supportIgnore := m.addSupport()
	allIgnore := make([]int, 0, len(ignore)+len(supportIgnore))
	allIgnore = append(allIgnore, ignore...)
	allIgnore = append(allIgnore, supportIgnore...)

	d_total := make([]float64, dof)
	beamForces := make([][6]float64, len(m.Beams))
	var lastK_T *sparse.Matrix

	for step := 1; step <= nr.Steps; step++ {
		factor := float64(step) / float64(nr.Steps)
		p_step := make([]float64, dof)
		for i := range p_step {
			p_step[i] = p_total[i] * factor
		}

		d_current := make([]float64, dof)
		copy(d_current, d_total)

		converged := false
		for iter := 0; iter < nr.MaxIter; iter++ {
			// compute axial forces from current displacements
			for bi := range m.Beams {
				b := m.Beams[bi]
				n0, n1 := b.N[0], b.N[1]
				u_global := mat.NewDense(6, 1, nil)
				for i := 0; i < 3; i++ {
					u_global.Set(i, 0, d_current[n0*3+i])
					u_global.Set(3+i, 0, d_current[n1*3+i])
				}
				tr := m.getCoordTransStiffBeam2d(bi)
				var u_local mat.Dense
				u_local.Mul(tr, u_global)
				k_local := m.getStiffBeam2d(bi)
				var f_local mat.Dense
				f_local.Mul(k_local, &u_local)
				for j := 0; j < 6; j++ {
					beamForces[bi][j] = f_local.At(j, 0)
				}
			}

			// build geometric stiffness
			lcGeom := &LoadCase{BeamForces: beamForces}
			K_geometric, _, err := m.assemblyK(func(pos int) *mat.Dense {
				return m.getGeometricBeam2d(pos, lcGeom)
			})
			if err != nil {
				return fmt.Errorf("assembly geometric stiffness at step %d, iter %d: %v", step, iter, err)
			}

			// K_T = K_linear - K_geometric
			K_T, err := sparse.Add(K_linear, K_geometric, 1.0, -1.0)
			if err != nil {
				return fmt.Errorf("add stiffness matrices at step %d, iter %d: %v", step, iter, err)
			}
			lastK_T = K_T

			var lu sparse.LU
			if err := lu.Factorize(K_T, allIgnore...); err != nil {
				return fmt.Errorf("LU factorization failed at step %d, iter %d: %v", step, iter, err)
			}

			d_new, err := lu.Solve(p_step)
			if err != nil {
				return fmt.Errorf("LU solve failed at step %d, iter %d: %v", step, iter, err)
			}

			// convergence: relative change in displacements
			normD := 0.0
			normDiff := 0.0
			for i := range d_new {
				diff := d_new[i] - d_current[i]
				normDiff += diff * diff
				normD += d_new[i] * d_new[i]
			}
			normD = math.Sqrt(normD)
			normDiff = math.Sqrt(normDiff)

			if normD == 0 || normDiff < nr.Tol*normD {
				copy(d_current, d_new)
				converged = true
				break
			}

			copy(d_current, d_new)

			if iter == nr.MaxIter-1 {
				return fmt.Errorf("fixed-point iteration did not converge at step %d, relative change: %e", step, normDiff/normD)
			}
		}

		if !converged {
			return fmt.Errorf("Newton-Raphson did not converge at step %d", step)
		}

		copy(d_total, d_current)
	}

	// store displacements
	lc.PointDisplacementGlobal = make([][3]float64, len(m.Points))
	for p := 0; p < len(m.Points); p++ {
		for i := 0; i < 3; i++ {
			lc.PointDisplacementGlobal[p][i] = d_total[3*p+i]
		}
	}

	// beam forces from final displacements
	lc.BeamForces = make([][6]float64, len(m.Beams))
	for bi, b := range m.Beams {
		u_global := mat.NewDense(6, 1, nil)
		for i := 0; i < 3; i++ {
			u_global.Set(i, 0, d_total[b.N[0]*3+i])
			u_global.Set(3+i, 0, d_total[b.N[1]*3+i])
		}
		tr := m.getCoordTransStiffBeam2d(bi)
		var z mat.Dense
		z.Mul(tr, u_global)
		kr := m.getStiffBeam2d(bi)
		s := mat.NewDense(6, 1, lc.BeamForces[bi][:])
		s.Mul(kr, &z)
	}

	// reactions
	lc.Reactions = make([][3]float64, len(m.Points))
	if lastK_T != nil {
		for pt := 0; pt < len(m.Points); pt++ {
			for i := 0; i < 3; i++ {
				if !m.Supports[pt][i] {
					continue
				}
				row := 3*pt + i
				react := 0.0
				sparse.Fkeep(lastK_T, func(ki, kj int, x float64) bool {
					if ki == row {
						react += x * d_total[kj]
					}
					return true
				})
				lc.Reactions[pt][i] = react
			}
		}
	}

	return nil
}

// arcLengthConfig contains parameters for arc-length control method.
type arcLengthConfig struct {
	Steps   int     // number of arc-length steps
	MaxIter int     // max iterations per step
	Tol     float64 // residual tolerance ‖R‖ < tol * ‖λ*f_ref‖
}

// nonlinearStaticArcLength performs geometrically nonlinear static analysis
// using Crisfield's spherical arc-length method.
//
// Reference: de Borst (2012) Ch.1.1, PyFEM RiksSolver.py; Crisfield (1991) Ch.9.
// Constraint: Δu^T * Δu + Δλ² = Δs²   (ψ = 1)
// Corrector:  δλ = (-Δu^T*δu_I - a0) / (Δu^T*δu_II + Δλ)
//
//	where a0 = (Δu^T*Δu + Δλ² - Δs²)/2
func nonlinearStaticArcLength(out io.Writer, m *Model, lc *LoadCase, al arcLengthConfig) (err error) {
	if out == nil {
		out = io.Discard
	}

	fmt.Fprintf(out, "\nNonlinear Static Analysis (Arc-length control)\n")

	defer func() {
		if err != nil {
			et := errors.New("Arc-length static analysis")
			_ = et.Add(err)
			err = et
		}
	}()
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("stacktrace from panic: %s\n", debug.Stack())
		}
	}()

	if al.Steps < 1 {
		return fmt.Errorf("arcLengthConfig.Steps must be >= 1, got %d", al.Steps)
	}
	if al.Tol <= 0 {
		return fmt.Errorf("arcLengthConfig.Tol must be > 0, got %e", al.Tol)
	}
	if al.MaxIter < 1 {
		return fmt.Errorf("arcLengthConfig.MaxIter must be >= 1, got %d", al.MaxIter)
	}

	out, err = prepare(out, m)
	if err != nil {
		return fmt.Errorf("prepare model: %v", err)
	}
	if err := lc.checkInputData(m); err != nil {
		return fmt.Errorf("check load case: %v", err)
	}

	dof := 3 * len(m.Points)

	K_linear, ignore, err := m.assemblyK(func(pos int) *mat.Dense {
		return m.getStiffBeam2d(pos)
	})
	if err != nil {
		return fmt.Errorf("assembly linear stiffness: %v", err)
	}

	f_ref, err := m.assemblyNodeLoad(lc)
	if err != nil {
		return fmt.Errorf("assembly reference load: %v", err)
	}
	for _, i := range ignore {
		if f_ref[i] != 0.0 {
			return fmt.Errorf("on direction %d load is not zero: %f", i, f_ref[i])
		}
	}

	supportIgnore := m.addSupport()
	allIgnore := make([]int, 0, len(ignore)+len(supportIgnore))
	allIgnore = append(allIgnore, ignore...)
	allIgnore = append(allIgnore, supportIgnore...)

	u := make([]float64, dof)
	lam := 0.0
	beamForces := make([][6]float64, len(m.Beams))

	// ignore map for convergence check on free DOFs only
	ignoreMap := make(map[int]struct{}, len(allIgnore))
	for _, idx := range allIgnore {
		ignoreMap[idx] = struct{}{}
	}

	// first step: solve linear problem to get scaling
	ds := 0.0
	{
		K0, _, _ := m.assemblyK(func(pos int) *mat.Dense { return m.getStiffBeam2d(pos) })
		var lu0 sparse.LU
		if lu0.Factorize(K0, allIgnore...) == nil {
			if d0, err := lu0.Solve(f_ref); err == nil {
				nrm := 0.0
				for i := range d0 {
					nrm += d0[i] * d0[i]
				}
				ds = math.Sqrt(nrm+1) / 20.0
			}
		}
	}
	if ds == 0 {
		ds = 0.5
	}

	fmt.Fprintf(out, "  ds = %.6f\n", ds)
	fmt.Fprintf(out, "\n  step  : load factor :  displ. apex :  residual :  iters\n")

	// previous step displacement increment, used for the predictor
	// direction so the equilibrium path is traced continuously.
	var prevΔu []float64

	for step := 1; step <= al.Steps; step++ {
		// state at the beginning of the step, restored on substep retry
		uStepStart := make([]float64, dof)
		copy(uStepStart, u)
		lamStepStart := lam

		var (
			converged  bool
			finalIter  int
			finalNormR float64
			Δu         []float64
		)

		// Adaptive substepping: if the corrector fails (typically when the
		// tangent stiffness becomes indefinite past a limit point) the step
		// is retried with a halved arc length.
		const maxRetry = 30
		for attempt := 0; attempt < maxRetry && !converged; attempt++ {
			copy(u, uStepStart)
			lam = lamStepStart

			// ---- predictor ----
			// compute tangent stiffness at current state
			for bi := range m.Beams {
				b := m.Beams[bi]
				n0, n1 := b.N[0], b.N[1]
				ug := mat.NewDense(6, 1, nil)
				for i := 0; i < 3; i++ {
					ug.Set(i, 0, u[n0*3+i])
					ug.Set(3+i, 0, u[n1*3+i])
				}
				tr := m.getCoordTransStiffBeam2d(bi)
				var ul mat.Dense
				ul.Mul(tr, ug)
				kl := m.getStiffBeam2d(bi)
				var fl mat.Dense
				fl.Mul(kl, &ul)
				for j := 0; j < 6; j++ {
					beamForces[bi][j] = fl.At(j, 0)
				}
			}
			lcGeom := &LoadCase{BeamForces: beamForces}
			Kg, _, _ := m.assemblyK(func(pos int) *mat.Dense { return m.getGeometricBeam2d(pos, lcGeom) })
			K_T, _ := sparse.Add(K_linear, Kg, 1.0, -1.0)
			var lu sparse.LU
			if err := lu.Factorize(K_T, allIgnore...); err != nil {
				return fmt.Errorf("LU at step %d: %v", step, err)
			}

			// predictor: δu_II = K_T^{-1} * f_ref
			δu_II, err := lu.Solve(f_ref)
			if err != nil {
				return fmt.Errorf("solve predictor at step %d: %v", step, err)
			}
			// Predictor direction: continue along the previous displacement
			// increment so the equilibrium path is traced continuously
			// (Crisfield). Without this the predictor may jump to the
			// opposite branch and the solver leaves the physical path.
			sign := 1.0
			if step > 1 {
				dot := 0.0
				for i := range δu_II {
					dot += δu_II[i] * prevΔu[i]
				}
				if dot < 0 {
					sign = -1.0
				}
			}

			// Δu₀ = sign * ds * δu_II / ‖(δu_II, 1)‖ (scaled to arc-length ds)
			nrm := 0.0
			for i := range δu_II {
				nrm += δu_II[i] * δu_II[i]
			}
			nrm = math.Sqrt(nrm + 1)
			Δλ0 := sign * ds / nrm
			Δu0 := make([]float64, dof)
			for i := range Δu0 {
				Δu0[i] = Δλ0 * δu_II[i]
			}

			// predictor update
			for i := range u {
				u[i] += Δu0[i]
			}
			lam += Δλ0

			// store predictor increment for constraint
			Δu = make([]float64, dof)
			copy(Δu, Δu0)
			Δλ := Δλ0

			// ---- corrector ----
			for iter := 0; iter < al.MaxIter; iter++ {

				// rebuild K_T at current u
				for bi := range m.Beams {
					b := m.Beams[bi]
					n0, n1 := b.N[0], b.N[1]
					ug := mat.NewDense(6, 1, nil)
					for i := 0; i < 3; i++ {
						ug.Set(i, 0, u[n0*3+i])
						ug.Set(3+i, 0, u[n1*3+i])
					}
					tr := m.getCoordTransStiffBeam2d(bi)
					var ul mat.Dense
					ul.Mul(tr, ug)
					kl := m.getStiffBeam2d(bi)
					var fl mat.Dense
					fl.Mul(kl, &ul)
					for j := 0; j < 6; j++ {
						beamForces[bi][j] = fl.At(j, 0)
					}
				}
				lcGeom = &LoadCase{BeamForces: beamForces}
				Kg, _, _ = m.assemblyK(func(pos int) *mat.Dense { return m.getGeometricBeam2d(pos, lcGeom) })
				K_T, _ = sparse.Add(K_linear, Kg, 1.0, -1.0)
				if err := lu.Factorize(K_T, allIgnore...); err != nil {
					return fmt.Errorf("LU at step %d, iter %d: %v", step, iter, err)
				}

				// internal forces: f_int = K_T * u
				f_int := make([]float64, dof)
				sparse.Fkeep(K_T, func(i, j int, x float64) bool {
					f_int[i] += x * u[j]
					return true
				})

				// residual
				R := make([]float64, dof)
				normR := 0.0
				normP := 0.0
				for i := range R {
					R[i] = lam*f_ref[i] - f_int[i]
					if _, ignored := ignoreMap[i]; ignored {
						R[i] = 0.0
						continue
					}
					normR += R[i] * R[i]
					normP += lam * f_ref[i] * lam * f_ref[i]
				}
				normR = math.Sqrt(normR)
				normP = math.Sqrt(normP)

				if normP > 0 && normR < al.Tol*normP {
					finalNormR = normR
					finalIter = iter
					converged = true
					break
				}
				finalNormR = normR

				// d1 = K_T^{-1} * f_ref
				d1, err := lu.Solve(f_ref)
				if err != nil {
					return fmt.Errorf("solve d1 at step %d, iter %d: %v", step, iter, err)
				}
				// d2 = K_T^{-1} * R
				d2, err := lu.Solve(R)
				if err != nil {
					return fmt.Errorf("solve d2 at step %d, iter %d: %v", step, iter, err)
				}

				// Crisfield spherical arc-length constraint:
				//   (Δu + δλ*d1 + d2)^T (Δu + δλ*d1 + d2) + (Δλ + δλ)² = ds²
				// which is a quadratic in δλ:
				//   a*δλ² + b*δλ + c = 0
				//   a = d1·d1 + 1
				//   b = 2*(Δu·d1 + d1·d2 + Δλ)
				//   c = Δu·Δu + 2*Δu·d2 + d2·d2 + Δλ² - ds²
				var d1d1, d1d2, d2d2 float64
				var ud1, ud2, uu float64
				for i := range d1 {
					d1d1 += d1[i] * d1[i]
					d1d2 += d1[i] * d2[i]
					d2d2 += d2[i] * d2[i]
					ud1 += Δu[i] * d1[i]
					ud2 += Δu[i] * d2[i]
					uu += Δu[i] * Δu[i]
				}
				a := d1d1 + 1.0
				b := 2.0 * (ud1 + d1d2 + Δλ)
				c := uu + 2.0*ud2 + d2d2 + Δλ*Δλ - ds*ds

				disc := b*b - 4.0*a*c
				if disc < 0 {
					disc = 0
				}
				sq := math.Sqrt(disc)
				δλ1 := (-b + sq) / (2.0 * a)
				δλ2 := (-b - sq) / (2.0 * a)

				// Select the root that keeps tracing the same branch: the one
				// whose resulting increment has the largest projection on the
				// previous step increment (or on the current step increment at
				// the first step).
				ref := Δu
				if step > 1 {
					ref = prevΔu
				}
				proj := func(dl float64) float64 {
					var s float64
					for i := range Δu {
						s += (Δu[i] + dl*d1[i] + d2[i]) * ref[i]
					}
					return s
				}
				δλ := δλ1
				if proj(δλ2) > proj(δλ1) {
					δλ = δλ2
				}

				if math.IsNaN(δλ) || math.IsInf(δλ, 0) {
					return fmt.Errorf("δλ is NaN/Inf at step %d, iter %d", step, iter)
				}

				// update increments
				for i := range Δu {
					Δu[i] += δλ*d1[i] + d2[i]
				}
				Δλ += δλ

				// update total
				for i := range u {
					u[i] += δλ*d1[i] + d2[i]
				}
				lam += δλ

				if iter == al.MaxIter-1 {
					// corrector did not converge: leave the loop and retry the
					// step with a smaller arc length.
					break
				}
			}

			// on failure, halve the arc length and retry the step
			if !converged {
				ds *= 0.5
				if ds < 1e-12 {
					return fmt.Errorf("arc-length step %d: arc length underflow", step)
				}
				continue
			}

			// step size control
			optIter := 6
			factor := math.Sqrt(float64(optIter) / float64(finalIter+1))
			if factor > 1.5 {
				factor = 1.5
			}
			if factor < 0.5 {
				factor = 0.5
			}
			ds *= factor

			fmt.Fprintf(out, "  %4d : %12.6f : %12.6f : %12.6e : %4d\n", step, lam, u[4], finalNormR, finalIter+1)

			// store increment of this step for the next predictor direction
			prevΔu = Δu
		}
	}

	// store results
	lc.PointDisplacementGlobal = make([][3]float64, len(m.Points))
	for p := 0; p < len(m.Points); p++ {
		for i := 0; i < 3; i++ {
			lc.PointDisplacementGlobal[p][i] = u[3*p+i]
		}
	}
	lc.BeamForces = make([][6]float64, len(m.Beams))
	for bi, b := range m.Beams {
		ug := mat.NewDense(6, 1, nil)
		for i := 0; i < 3; i++ {
			ug.Set(i, 0, u[b.N[0]*3+i])
			ug.Set(3+i, 0, u[b.N[1]*3+i])
		}
		tr := m.getCoordTransStiffBeam2d(bi)
		var z mat.Dense
		z.Mul(tr, ug)
		kr := m.getStiffBeam2d(bi)
		s := mat.NewDense(6, 1, lc.BeamForces[bi][:])
		s.Mul(kr, &z)
	}
	lc.Reactions = make([][3]float64, len(m.Points))
	sparse.Fkeep(K_linear, func(i, j int, x float64) bool {
		for pt := 0; pt < len(m.Points); pt++ {
			for k := 0; k < 3; k++ {
				if m.Supports[pt][k] && i == pt*3+k {
					lc.Reactions[pt][k] += x * u[j]
				}
			}
		}
		return true
	})

	return nil
}
