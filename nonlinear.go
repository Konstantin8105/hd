package hd

import (
	"fmt"
	"math"
)

// Explanations:
//	D - displacement
//	F - load
//	λ - load proportionality factor (LPF). The dimensionless ``load'' vector
//	F = f(D) - result precision function
//	K = dF/dD = f(F,D) - stiffness matrix, function of displacement and load
//	invert(K) - invertion stiffness matrix
//
// Linear solution displacement by external load:
//	1. Do = 0 - initialize zero displacement
//	   Fo = 0 - initialize zero load
//	   Fe - result load
//	   De - result displacement
//	2. Ko = dF/dD for (Do,Fo)=(0,0)
//	3. Solve : Ko*De = Fe
//	   invert(Ko)*Fe = De
//
// Non-linear solution by "Increment(Euler first order) method":
//	1. Do,   Fo - initialize displacement, load
//	   De,   Fe - result displacement, load
//	   dDe, dFe - step increment displacement, load
//	2. Ko = dF/dD for (Do,Fo)=(0,0)
//	3. dFe = Fe - Fo
//	4. Solve : Ko*dDe = dFe
//	   invert(Ko)*dFe = dDe
//	5. De = Do + dDe
//	   De = Do + invert(Ko)*dFe
//	6. Break criteria by (De, Fe, amount steps)
//	7. Next step:
//	   Do = De
//	   Fo = Fe
//	   Fe - generate new loads
//	   Go to step 1.
//
//
//
// TODO: cylinder arc-lenght method
// TODO: Quasi-Newton method
// TODO: Runge-Cutta 4 order
// TODO: Displacement control
// TODO: Work control
//
// TODO: Automatic load increment
//	L - lambda load value
//	I - amount iterations
//	n - factor, typically n = 1.0/2.0...1.0
//	n subscript - for next iteration
//	o subscript - for last iteration
// dLn = dLo*(In/Io)^n

type (
	Displacements []float64
	Forces        []float64
	solverFunc    func(
		Do Displacements,
		Fo, Fe Forces,
		Kstiff func(F Forces, D Displacements) (matrix, error),
		Update func(F Forces, D Displacements, K matrix),
		Stop func(iter uint, dF, F Forces, dD, D Displacements) (bool, error),
	) (iterations uint, err error)
)

type matrix interface {
	Solve(dF Forces) (Displacements, error)
	Mul(dD Displacements) (Forces, error)
	Determinant() (float64, error)
}

// nr is non-linear solution by "The Newton-Raphson method"(Load control):
//
// Taylor Series:
//
//	f(x) = f(a)+f`(a)/1!*(x-a)+f``(a)/2!*(x-a)^2+f```(a)/3!*(x-a)^3+...
//
//	0. Theory:
//	   base formula:
//	   f(x) = f(a)+f`(a)/1!*(x-a)+f``(a)/2!*(x-a)^2
//	   avoid factorial:
//	   f(x) = f(a)+f`(a)*(x-a)+1/2*f``(a)*(x-a)^2
//	   find f(b) by f(a):
//	   f(b) = f(a)+f`(a)*(b-a)+1/2*f``(a)*(b-a)^2
//	   find f(c) by f(b):
//	   f(c) = f(b)+f`(b)*(c-b)+1/2*f``(b)*(c-b)^2
//	   find f(c) by f(a):
//	   f(c) = f(a)+f`(a)*(b-a)+1/2*f``(a)*(b-a)^2 +
//	              +f`(b)*(c-b)+1/2*f``(b)*(c-b)^2
//	   avoid f``:
//	   f(c) = f(a)+f`(a)*(b-a)+f`(b)*(c-b)
//	   use infinite formula
//	   f(c) = f(a)+f`(a)*(b-a)+f`(b)*(c-b)+...
//	   Problems:
//	   * with negative f`(x)
//	   * f`(x) = 0
//	   * f`(x) = infinite
//	   * Snap-Through under load control
//	   * Snap-Through under displacement control
//	1. Do, Fo - initialize displacement, load
//	   De, Fe - result displacement, load
//	   i = 0
//	   D_(i) = Do
//	   F_(i) = Fe
//	2. Substep for D_(i), F_(i)
//	3. Ki = dFi/dDi for (Di,Fi)
//	4. Solve : Ki*dDi = Fe
//	   invert(Ki)*Fe  = dDi
//	   dD_(i) = invert(Ki)*Fe
//	5. D_(i+1) = D_(i)+dD_(i)
//	   D_(i+1) = D_(i)+invert(Ki)*Fe
//	6. Break criteria by (D_(i+1), amount steps)
//	7. Next step:
//	   i = i+1
//	   Go to step 2.
func nr(
	Do Displacements,
	Fo, Fe Forces,
	Kstiff func(F Forces, D Displacements) (matrix, error),
	Update func(F Forces, D Displacements, K matrix),
	Stop func(iter uint, dF, F Forces, dD, D Displacements) (bool, error),
) (iterations uint, err error) {
	// assemble stiffness matrix with geometric nonlinearity
	K, err := Kstiff(Fo, Do)
	if err != nil {
		return
	}
	size := len(Do)
	for iterations = 1; ; iterations++ {
		// load increment
		dF := make([]float64, size)
		for i := range dF {
			dF[i] = Fe[i] - Fo[i]
		}
		// displacement increment
		dD, err := K.Solve(dF)
		if err != nil {
			return iterations, err
		}
		for i := range Do {
			Do[i] += dD[i]
		}
		// assemble stiffness matrix with geometric nonlinearity
		K, err = Kstiff(Fe, Do)
		if err != nil {
			return iterations, err
		}
		// update load
		dF, err = K.Mul(dD)
		if err != nil {
			return iterations, err
		}
		for i := range Fo {
			Fo[i] += dF[i]
		}
		// update load case data
		Update(Fo, Do, K)
		// stop criteria
		stop, err := Stop(iterations, dF, Fo, dD, Do)
		if err != nil {
			return iterations, err
		}
		if stop {
			break
		}
	}
	return
}

func nrs(
	substeps uint,
	inF solverFunc,
) (outF solverFunc) {
	return func(
		Do Displacements,
		Fo, Fe Forces,
		Kstiff func(F Forces, D Displacements) (matrix, error),
		Update func(F Forces, D Displacements, K matrix),
		Stop func(iter uint, dF, F Forces, dD, D Displacements) (bool, error),
	) (iterations uint, err error) {
		// Example:
		// input:
		//	Fo = 1
		//	Fe = 5
		//	substeps = 3
		// output:
		//	dF = (Fe-Fo)/(substeps+1) = (5-1)/(3+1) = 4/4 = 1
		//	loop:
		//	s = 0.	Fe = Fe - (substeps-s) * dF = 5-(3-0)*1 = 2
		//	s = 1.	Fe = Fe - (substeps-s) * dF = 5-(3-1)*1 = 3
		//	s = 2.	Fe = Fe - (substeps-s) * dF = 5-(3-2)*1 = 4
		//	s = 3.	Fe = Fe - (substeps-s) * dF = 5-(3-3)*1 = 5
		size := len(Do)
		dF := make([]float64, size)
		for i := range dF {
			dF[i] = (Fe[i] - Fo[i]) / float64(substeps+1)
		}
		endF := make([]float64, len(Fe))
		copy(endF, Fe)
		var s uint
		for ; s <= substeps; s++ {
			for i := range Fe {
				Fe[i] = endF[i] - float64(substeps-s)*dF[i]
			}
			iter, err := inF(Do, Fo, Fe, Kstiff, Update, Stop)
			iterations += iter
			if err != nil {
				return iterations, err
			}
		}
		return
	}
}

type arc struct {
	// Hyperellipsoid ratio
	Ksi float64
	// Radius
	Radius float64
}

// Non-linear solution by "The `spherical arc-lenght` method":
//  0. Theory:
//     Based on research
//     Nonlinear Analysis of Structures
//     The Arc Length Method: Formulation, Implementation and Applications
//     Nikolaos Vasios
//     Formula (2.12). (∆u + δu)T*(∆u + δu) + ψ^2*(∆λ + δλ)^2*(𝐪T * 𝐪) = ∆l^2
//     Formula (2.14). δu = δū + δλ*δut
//     δū = -invert[KT](uo+Δu) * (Fint*(uo+Δu)-(λo+Δλ)*𝐪)
//     δut = invert[KT](uo+Δu) * 𝐪
//     Formula (2.15). 𝛼1*δλ^2 + 𝛼2*δλ + 𝛼3 = 0
//     𝛼1 = δutT*δut + ψ^2*(𝐪T * 𝐪)
//     𝛼2 = 2*(∆u+δū)*δut+2*ψ^2*∆λ*(𝐪T * 𝐪)
//     𝛼3 = (∆u + δū)T*(∆u + δū)+ψ^2*∆λ^2*(𝐪T * 𝐪)-∆l^2
//  1. Do, Fo - initialize displacement, load
//     λo     - initialize load proportionality factor
//     i = 0
func (s arc) solver(
	Do Displacements,
	Fo, Fe Forces,
	iKstiff func(F Forces, D Displacements) (matrix, error),
	iUpdate func(F Forces, D Displacements, K matrix),
	iStop func(iter uint, dF, F Forces, dD, D Displacements) (bool, error),
) (iterations uint, err error) {
	// move point (Do, Fo) to (0,0)
	Kstiff := func(F Forces, D Displacements) (matrix, error) {
		return iKstiff(summa(F, Fo), summa(D, Do))
	}
	Update := func(F Forces, D Displacements, K matrix) {
		iUpdate(summa(F, Fo), summa(D, Do), K)
	}
	Stop := func(iter uint, dF, F Forces, dD, D Displacements) (bool, error) {
		return iStop(iter, dF, summa(F, Fo), dD, summa(D, Do))
	}
	for i := range Fe {
		Fe[i] = Fe[i] - Fo[i]
	}
	// input data
	var (
		// Degree of freedom.
		dof = len(Do)
		// Initial displacement.
		// Default value is zero.
		u = make([]float64, dof)
		// End force.
		𝐪 = Fe
		// Load proportionality factor (LPF).
		// The dimensionless ``load'' vector.
		// Default value is zero.
		λ float64
		// Hyperellipsoid ratio
		// if 𝜓 = 1, then Spherical Arc-Length Method
		𝜓 float64 = s.Ksi
		// Radius
		Δl float64 = s.Radius
	)
	// error handling input data
	if 𝜓 <= 0 {
		err = fmt.Errorf("Hyperellipsoid ratio %.5e is not valid", 𝜓)
		return
	}
	if Δl <= 0 {
		err = fmt.Errorf("Radius %.5e is not valid", Δl)
		return
	}
	// prepare output data
	defer func() {
		for i := range Do {
			Do[i] += u[i]
		}
		for i := range Fe {
			Fo[i] = λ * Fe[i]
		}
	}()
	// iteration
	for iterations = 1; ; iterations++ {
		// prepare new data for next iteration
		var (
			Δu = make([]float64, dof)

			δu       []float64
			δu1, δu2 []float64

			det    float64
			Δλ     float64
			fcheck float64 = 1.0 // any default positive value

			δλ       float64
			δλ1, δλ2 float64

			Kt matrix
		)

		begin := func() (err error) {
			Kt, err = Kstiff(scale(λ, 𝐪), summa(u, Δu))
			if err != nil {
				return err
			}
			// For formula (2.14):
			// δū = -invert[KT](uo+Δu) * (Fint*(uo+Δu)-(λo+Δλ)*𝐪)
			// value of Fint is precision value and for we cannot find them
			// only by Jacobi matrix.
			// Fint*(uo+Δu)-(λo+Δλ)*𝐪 is equal R(uo+Δu), but
			// theoretically R(uo+Δu) = 0, then vector is zero always:
			// if Kt == nil {
			// Kt, err = Kstiff(scale(λ, 𝐪), u)
			// if err != nil {
			// 	return err
			// }
			// }
			δū := make([]float64, dof)
			// δū, err = Kt.Solve(summa(scale(λ+Δλ, 𝐪), scale(-1, scale(λ+Δλ, 𝐪))))
			//
			// δū, err = Kt.Solve(scale(Δλ, 𝐪))
			// if err != nil {
			// 	return err
			// }
			// δū = summa(δū, scale(-1, Δu))
			// δū = scale(-1.0, δū)
			// log.Println(δū)
			// -np.dot(dfinv,f) //  TODO  main problem find this
			// TODO : find the solution
			δut, err := Kt.Solve(𝐪)
			if err != nil {
				return err
			}

			// Formula (2.12):
			// (∆u + δu)T*(∆u + δu) + ψ^2*(∆λ + δλ)^2*(𝐪T * 𝐪) = ∆l^2
			//
			// Formula (2.14)
			// δu = δū + δλ*δut
			//
			// Formula (2.15)
			// 𝛼1*δλ^2 + 𝛼2*δλ + 𝛼3 = 0
			//
			// symbolic math:
			// pow(deltau+(δu_+δλ*δut),2)+ψ2*pow(deltaλ+δλ,2)*(q2)-l2
			//
			// deltau*deltau + 2*deltau*δu_ + δu_*δu_ + 2*deltau*δut*δλ + \
			// ::::::::::::::::::::::::::::::::::::::   ---------------   \
			// 2*δu_*δut*δλ + δut*δut*δλ*δλ + deltaλ*deltaλ*q2*ψ2 +       \
			// -------------  =============   :::::::::::::::::::         \
			// 2*deltaλ*q2*δλ*ψ2 + q2*δλ*δλ*ψ2 - l2                       \
			// -----------------   ============::::                       \
			//
			// 𝛼1 = δutT*δut + ψ^2*(𝐪T * 𝐪)
			// 𝛼2 = 2*(∆u+δū)*δut+2*ψ^2*∆λ*(𝐪T * 𝐪)
			// 𝛼3 = (∆u + δū)T*(∆u + δū)+ψ^2*∆λ^2*(𝐪T * 𝐪)-∆l^2
			//
			var (
				// calculate the coefficients of the polynomial
				𝛼1 = dot(δut, δut) +
					math.Pow(𝜓, 2.0)*dot(𝐪, 𝐪)
				𝛼2 = 2.0*dot(summa(Δu, δū), δut) +
					2.0*Δλ*math.Pow(𝜓, 2.0)*dot(𝐪, 𝐪)
				𝛼3 = dot(summa(Δu, δū), summa(Δu, δū)) +
					math.Pow(𝜓, 2.0)*math.Pow(Δλ, 2.0)*dot(𝐪, 𝐪) -
					math.Pow(Δl, 2.0)

				// determinant
				D = 𝛼2*𝛼2 - 4.*𝛼1*𝛼3
			)
			if 0.0 < D {
				// acceptable 2 solutions
				δλ1 = (-𝛼2 - math.Sqrt(D)) / (2.0 * 𝛼1)
				δλ2 = (-𝛼2 + math.Sqrt(D)) / (2.0 * 𝛼1)
			} else {
				δλ1 = -𝛼2 / (2.0 * 𝛼1)
				δλ2 = -𝛼2 / (2.0 * 𝛼1)
				// panic((fmt.Errorf("not implemented: (%e,%e,%e) - %e",
				// 	𝛼1, 𝛼2, 𝛼3, D)))
				// TODO : check coverage for that part of code : D < 0.0
			}
			// After checking - acceptable swap the results, but no need
			// δλ2, δλ1 = δλ1, δλ2

			// Formula (2.14):
			// δu = δū + δλ*δut
			δu1 = summa(δū, scale(δλ1, δut))
			δu2 = summa(δū, scale(δλ2, δut))

			// calculate determinant matrix of stiffiners
			det, err = Kt.Determinant()
			if err != nil {
				return err
			}
			return nil
		}
		if err = begin(); err != nil {
			return iterations, err
		}

		if math.Signbit(det) == math.Signbit(δλ1) {
			δu, δλ = δu1, δλ1
		} else {
			δu, δλ = δu2, δλ2
		}

		finish := func() {
			Δu = summa(Δu, δu)
			Δλ = Δλ + δλ
			fcheck = math.Max(norm(δu), math.Abs(δλ))
			// log.Printf("%.3e %v\n", fcheck, fcheck < 1e-6)
			// log.Println(δu, δλ)
		}
		finish()
		// substeps
		// log.Printf("%.3e %v\n", fcheck, fcheck < 1e-6)
		// log.Println(δu, δλ)
		for 1e-6 < fcheck {
			// log.Println(";;;;;;;")
			if err = begin(); err != nil {
				return iterations, err
			}

			daomag := dot(Δu, Δu)
			if daomag == 0.0 {
				if math.Signbit(Δλ+δλ1) == math.Signbit(det) {
					δu, δλ = δu1, δλ1
				} else {
					δu, δλ = δu2, δλ2
				}
			} else {
				// Formula (2.16)
				DOT1 := dot(summa(Δu, δu1), Δu) + math.Pow(𝜓, 2)*Δλ*(Δλ+δλ1)*dot(𝐪, 𝐪)
				DOT2 := dot(summa(Δu, δu2), Δu) + math.Pow(𝜓, 2)*Δλ*(Δλ+δλ2)*dot(𝐪, 𝐪)
				if DOT1 > DOT2 {
					δu, δλ = δu1, δλ1
				} else {
					δu, δλ = δu2, δλ2
				}
			}
			if δλ1 == δλ2 {
				δu, δλ = δu1, δλ1
			}
			// log.Println("SUBSTEP: ", δu, δλ)
			// log.Println("daomag : ", daomag)
			finish()
		}
		// store values
		u = summa(u, Δu)
		λ += Δλ
		// renaming
		dD, D := Δu, u
		dF, F := scale(Δλ, Fe), scale(λ, Fe)
		// update load case data
		Update(F, D, Kt)
		// stop criteria
		stop, err := Stop(iterations, dF, F, dD, D)
		if err != nil {
			return iterations, err
		}
		if stop {
			break
		}
	}
	return
}

func dot(a, b []float64) float64 {
	var res float64
	for i := range a {
		res += a[i] * b[i]
	}
	return res
}

func scale(f float64, a []float64) []float64 {
	size := len(a)
	res := make([]float64, size)
	for i := 0; i < size; i++ {
		res[i] = f * a[i]
	}
	return res
}

func norm(a []float64) float64 {
	var res float64
	for _, v := range a {
		res += v * v
	}
	return math.Sqrt(res)
}

func summa(a, b []float64) []float64 {
	size := len(a)
	res := make([]float64, size)
	for i := 0; i < size; i++ {
		res[i] = a[i] + b[i]
	}
	return res
}
