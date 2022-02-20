package hd

import (
	"fmt"
	"math"
	"os"
)

// Do: 0.0,
// Fo: 0.0,
var nlCases = []struct {
	name string
	F    func(D Displacements) []float64
	K    func(F Forces, D Displacements) [][]float64
	De   Displacements
	Fe   Forces
}{
	{
		name: "parabola /\\",
		F: func(D Displacements) []float64 {
			return []float64{-math.Pow(D[0], 2) + 30*D[0]}
		},
		K: func(F Forces, D Displacements) [][]float64 {
			return [][]float64{{-2*D[0] + 30}}
		},
		De: Displacements([]float64{7.5}),
		Fe: Forces([]float64{168.75}),
	},
	{
		name: "parabola /\\ near top",
		F: func(D Displacements) []float64 {
			return []float64{-math.Pow(D[0], 2) + 30*D[0]}
		},
		K: func(F Forces, D Displacements) [][]float64 {
			return [][]float64{{-2*D[0] + 30}}
		},
		De: Displacements([]float64{10.0}),
		Fe: Forces([]float64{200.0}),
	},
	{
		name: "parabola /\\ after top",
		F: func(D Displacements) []float64 {
			return []float64{-math.Pow(D[0], 2) + 30*D[0]}
		},
		K: func(F Forces, D Displacements) [][]float64 {
			return [][]float64{{-2*D[0] + 30}}
		},
		De: Displacements([]float64{28.0}),
		Fe: Forces([]float64{56.0}),
	},
	{
		name: "curve before top",
		F: func(D Displacements) []float64 {
			return []float64{-0.06*math.Pow(D[0], 3) + 1.2*math.Pow(D[0], 2) + 3*D[0]}
		},
		K: func(F Forces, D Displacements) [][]float64 {
			return [][]float64{{-0.06*3*math.Pow(D[0], 2) + 1.2*2*D[0] + 3}}
		},
		De: Displacements([]float64{10.0}),
		Fe: Forces([]float64{90.0}),
	},
	{
		name: "curve after top",
		F: func(D Displacements) []float64 {
			return []float64{-0.06*math.Pow(D[0], 3) + 1.2*math.Pow(D[0], 2) + 3*D[0]}
		},
		K: func(F Forces, D Displacements) [][]float64 {
			return [][]float64{{-0.06*3*math.Pow(D[0], 2) + 1.2*2*D[0] + 3}}
		},
		De: Displacements([]float64{20.0}),
		Fe: Forces([]float64{60.0}),
	},
}

var nlSolvers = []struct {
	name   string
	solver solverFunc
}{
	{
		name:   "The Newton-Raphson method(Load control)",
		solver: nr,
	},
	{
		name:   "NR method with 2 subteps",
		solver: nrs(2, nr),
	},
	{
		name:   "NR method with 10 subteps",
		solver: nrs(10, nr),
	},
	{
		name:   "NR method with 1000 subteps",
		solver: nrs(1000, nr),
	},
	{
		name:   "NR method with 100000 subteps",
		solver: nrs(100000, nr),
	},
}

func norm(dd []float64) (res float64) {
	for _, d := range dd {
		res = d * d
	}
	return math.Sqrt(res)
}

func ratio(e, o []float64) (res float64) {
	return math.Abs(norm(e)-norm(o)) / math.Max(norm(e), norm(o))
}

type DS struct { // double slice
	value [][]float64
}

func (K DS) Mul(dD Displacements) (Forces, error) {
	ndof := len(K.value)
	res := make([]float64, ndof)
	for i := 0; i < ndof; i++ {
		for k := 0; k < ndof; k++ {
			res[k] += K.value[k][i] * dD[i]
		}
	}
	return res, nil
}

func (K DS) Solve(dF Forces) (Displacements, error) {
	invKd := func(Kd [][]float64) (inv [][]float64) {
		size := len(Kd)
		inv = make([][]float64, size)
		for i := range inv {
			inv[i] = make([]float64, size)
		}
		switch size {
		case 1:
			inv[0][0] = 1.0 / Kd[0][0]
		case 2:
			var (
				a = Kd[0][0]
				b = Kd[0][1]
				c = Kd[1][0]
				d = Kd[1][1]
			)
			det := 1. / (a*d - b*c)
			inv[0][0] = det * d
			inv[0][1] = det * (-b)
			inv[1][0] = det * (-c)
			inv[1][1] = det * a
		default:
			panic("")
		}
		// TODO: divide by zero
		return
	}(K.value)
	k := DS{value: invKd}
	res, err := k.Mul(Displacements(dF))
	if err != nil {
		panic(err)
	}
	return Displacements(res), nil
}

func ExampleNonlinear() {
	for _, c := range nlCases {
		for _, s := range nlSolvers {
			// header
			fmt.Fprintf(os.Stdout, "\n")
			fmt.Fprintf(os.Stdout, "Solver : %s\n", s.name)
			fmt.Fprintf(os.Stdout, "Case   : %s\n", c.name)
			// run
			Do := make([]float64, len(c.Fe))
			Fo := make([]float64, len(c.Fe))
			Kstiff := func(F Forces, D Displacements) (matrix, error) {
				k := DS{value: c.K(F, D)}
				return k, nil
			}
			Update := func(F Forces, D Displacements, K matrix) {
			}
			Stop := func(iter uint, dF, F Forces, dD, D Displacements) (stop bool, err error) {
				const tol = 1e-8
				if 100000 < iter {
					err = fmt.Errorf("Too much iterations")
					return
				}
				if norm(dF) < tol {
					return true, nil
				}
				if norm(dD) < tol {
					return true, nil
				}
				if norm(c.Fe) < norm(F) && norm(c.De) < norm(D) {
					return true, nil
				}
				return
			}
			iterations, err := s.solver(Do, Fo, c.Fe, Kstiff, Update, Stop)
			if err != nil {
				fmt.Fprintf(os.Stdout, "error  : %v\n", err)
				continue
			}
			// compare results
			fmt.Fprintf(os.Stdout, "iters   = %d\n", iterations)
			// big error
			const tol = 1e-2
			if r := ratio(Do, c.De); tol < r {
				fmt.Fprintf(os.Stdout, "error  : tolerance displacements %.3e\n", r)
				continue
			}
			if r := ratio(Fo, c.Fe); tol < r {
				fmt.Fprintf(os.Stdout, "error  : tolerance forces %.3e\n", r)
				continue
			}
			// bottom
			fmt.Fprintf(os.Stdout, "Do      = %.5f\n", Do)
			fmt.Fprintf(os.Stdout, "De      = %.5f\n", c.De)
			fmt.Fprintf(os.Stdout, "norm(D) = %.2e\n", ratio(Do, c.De))
			fmt.Fprintf(os.Stdout, "Fo      = %.5f\n", Fo)
			fmt.Fprintf(os.Stdout, "Fe      = %.5f\n", c.Fe)
			fmt.Fprintf(os.Stdout, "norm(F) = %.2e\n", ratio(Fo, c.Fe))
		}
	}

	// Output:
}
