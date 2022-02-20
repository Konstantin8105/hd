package hd

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
// Taylor Series:
//	f(x) = f(a)+f`(a)/1!*(x-a)+f``(a)/2!*(x-a)^2+f```(a)/3!*(x-a)^3+...
//
// Non-linear solution by "The Newton-Raphson method"(Load control):
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
//
// Non-linear solution by "The `spherical arc-lenght` method":
// 0. Theory:
//	Based on research
//			Nonlinear Analysis of Structures
//			The Arc Length Method: Formulation, Implementation and Applications
//			Nikolaos Vasios
//	   Formula (2.12). (∆u + δu)T*(∆u + δu) + ψ^2*(∆λ + δλ)^2*(𝐪T * 𝐪) = ∆l^2
//	   Formula (2.14). δu = δū + δλ*δut
//	                   δū = -invert[KT](uo+Δu) * (Fint*(uo+Δu)-(λo+Δλ)*𝐪)
//	                   δut = invert[KT](uo+Δu) * 𝐪
//	   Formula (2.15). 𝛼1*δλ^2 + 𝛼2*δλ + 𝛼3 = 0
//	                   𝛼1 = δutT*δut + ψ^2*(𝐪T * 𝐪)
//	                   𝛼2 = 2*(∆u+δū)*δut+2*ψ^2*∆λ*(𝐪T * 𝐪)
//	                   𝛼3 = (∆u + δū)T*(∆u + δū)+ψ^2*∆λ^2*(𝐪T * 𝐪)-∆l^2
//	1. Do, Fo - initialize displacement, load
//	   λo     - initialize load proportionality factor
//	   i = 0
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
		Kstiff func(F Forces, D Displacements) (interface{}, error),
		Solve func(K interface{}, F Forces) (Displacements, error),
		Update func(F Forces, D Displacements, K interface{}),
		Mul func(K interface{}, dD Displacements) (Forces, error),
		Stop func(iter uint, dF Forces, dD Displacements) (bool, error),
	) (iterations uint, err error)
)

func nr(
	Do Displacements,
	Fo, Fe Forces,
	Kstiff func(F Forces, D Displacements) (interface{}, error),
	Solve func(K interface{}, F Forces) (Displacements, error),
	Update func(F Forces, D Displacements, K interface{}),
	Mul func(K interface{}, dD Displacements) (Forces, error),
	Stop func(iter uint, dF Forces, dD Displacements) (bool, error),
) (iterations uint, err error) {
	// assemble stiffness matrix with geometric nonlinearity
	K, err := Kstiff(Fo, Do)
	if err != nil {
		return
	}
	size := len(Do)
	for iterations = 0; ; iterations++ {
		// load increment
		dF := make([]float64, size)
		for i := range Fo {
			dF[i] = Fe[i] - Fo[i]
		}
		// displacement increment
		dD, err := Solve(K, dF)
		if err != nil {
			return iterations, err
		}
		for i := range Do {
			Do[i] += dD[i]
		}
		// assemble stiffness matrix with geometric nonlinearity
		K, err := Kstiff(Fe, Do)
		if err != nil {
			return iterations, err
		}
		// update load
		dF, err = Mul(K, dD)
		if err != nil {
			return iterations, err
		}
		for i := range Do {
			Fo[i] += dF[i]
		}
		// update load case data
		Update(Fo, Do, K)
		// stop criteria
		stop, err := Stop(iterations, dF, dD)
		if stop {
			break
		}
		if err != nil {
			return iterations, err
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
		Kstiff func(F Forces, D Displacements) (interface{}, error),
		Solve func(K interface{}, F Forces) (Displacements, error),
		Update func(F Forces, D Displacements, K interface{}),
		Mul func(K interface{}, dD Displacements) (Forces, error),
		Stop func(iter uint, dF Forces, dD Displacements) (bool, error),
	) (iterations uint, err error) {
		dF := make([]float64, len(Fo))
		for i := range dF {
			dF[i] = (Fe[i] - Fo[i]) / float64(substeps+1)
		}
		endF := make([]float64, len(Fe))
		copy(endF, Fe)
		var s uint
		for ; s <= substeps; s++ {
			for j := range endF {
				Fe[j] = endF[j] * float64(s+1) / float64(substeps+1)
			}
			iter, err := inF(Do, Fo, Fe, Kstiff, Solve, Update, Mul, Stop)
			iterations += iter
			if err != nil {
				return iterations, err
			}
		}
		return
	}
}
