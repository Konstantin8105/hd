package hd

import (
	"fmt"
	"math"
	"os"
)

// func TestNonlinearCases(t *testing.T) {
// 	for i := range nlCases {
// 		t.Run(nlCases[i].name, func(t *testing.T) {
// 			tol := 1e-8
// 			dif := make([]float64, len(nlCases[i].De))
// 			for h := range dif {
// 				dif[h] = nlCases[i].Fe[h] - nlCases[i].F(nlCases[i].De)[h]
// 			}
// 			if tol < norm(dif) {
// 				t.Errorf("Big diffirence: %.e", dif)
// 			}
// 		})
// 	}
// }

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
		name: "curve _/ before top",
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
		name: "curve _/ after top",
		F: func(D Displacements) []float64 {
			return []float64{-0.06*math.Pow(D[0], 3) + 1.2*math.Pow(D[0], 2) + 3*D[0]}
		},
		K: func(F Forces, D Displacements) [][]float64 {
			return [][]float64{{-0.06*3*math.Pow(D[0], 2) + 1.2*2*D[0] + 3}}
		},
		De: Displacements([]float64{20.0}),
		Fe: Forces([]float64{60.0}),
	},
	{
		name: "2 unknown curve",
		F: func(D Displacements) []float64 {
			return []float64{
				10*D[0] + 0.4*math.Pow(D[1], 3) - 5*math.Pow(D[1], 2),
				0.4*math.Pow(D[0], 3) - 3*math.Pow(D[0], 2) + 10*D[1],
			}
		},
		K: func(F Forces, D Displacements) [][]float64 {
			return [][]float64{
				{10, 1.2*math.Pow(D[1], 2) - 10*D[1]},
				{1.2*math.Pow(D[0], 2) - 6*D[0], 10},
			}
		},
		De: Displacements([]float64{2.340759528952, 1.593101851107}),
		Fe: Forces([]float64{40, 15}),
	},
}

var nlSolvers = []struct {
	name   string
	solver solverFunc
}{
	// NR
	{
		name:   "The Newton-Raphson method(Load control)",
		solver: nr,
	},
	// NR substep
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
	// Arc method
	{
		name:   "Arc method",
		solver: arc{Ksi: 1.0, Radius: 0.01}.solver,
	},
}

// func printData(data []row, filename string,
// 	q []float64,
// 	Kt func(a []float64) (df [][]float64),
// 	F func(x []float64, lambda float64) []float64) {
// 	// gnuplot graph
// 	// plot "data.txt" using 2:1 title "rotation", \
// 	//      "data.txt" using 3:1 title "vertical disp"
// 	var buf bytes.Buffer
// 	yintegral := make([]float64, len(q))
// 	for index, r := range data {
// 		fmt.Fprintf(&buf, "%.12f", r.lambda)
// 		for i := range r.u {
// 			fmt.Fprintf(&buf, " %.12f", r.u[i])
// 		}
// 		if F != nil {
// 			ys := F(r.u, 0) // r.lambda)
// 			if 0 < index {
// 				yintegral = summa(
// 					yintegral,
// 					dotm(Kt(r.u), summa(data[index].u, scale(-1, data[index-1].u))))
// 			}
// 			for i := range ys {
// 				fmt.Fprintf(&buf, " %.12f %.12f", yintegral[i], ys[i])
// 			}
// 		}
// 		fmt.Fprintf(&buf, "\n")
// 	}
// 	if err := os.WriteFile(filename, buf.Bytes(), 0644); err != nil {
// 		panic(err)
// 	}
// }

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

func (K DS) Determinant() (float64, error) {
	size := len(K.value)
	switch size {
	case 1:
		return K.value[0][0], nil
	case 2:
		var (
			a = K.value[0][0]
			b = K.value[0][1]
			c = K.value[1][0]
			d = K.value[1][1]
		)
		return (a*d - b*c), nil
	}
	panic("")
}

// TODO gonum Dense
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

func (_ DS) Middle(K1, K2 matrix) (matrix, error) {
	K1d := K1.(DS)
	K2d := K2.(DS)
	size := len(K1d.value)
	var mid DS
	mid.value = make([][]float64, size)
	for i := 0; i < size; i++ {
		mid.value[i] = make([]float64, size)
	}
	for i := 0; i < size; i++ {
		for j := 0; j < size; j++ {
			mid.value[i][j] = 0.5*K1d.value[i][j] + 0.5*K2d.value[i][j]
		}
	}
	return mid, nil
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
				defer func() {
					if err != nil {
						err = fmt.Errorf("%v. Last coordinate: D=%.3e,F=%.3e", err, D, F)
					}
				}()
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
				for i := range D {
					if 1e6 < math.Abs(D[i]) {
						err = fmt.Errorf("Too big displacement %.e", D[i])
						return
					}
				}
				for i := range F {
					if 1e6 < math.Abs(F[i]) {
						err = fmt.Errorf("Too big force %.e", F[i])
						return
					}
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
	// Solver : The Newton-Raphson method(Load control)
	// Case   : parabola /\
	// iters   = 8
	// error  : tolerance displacements 3.769e-01
	//
	// Solver : NR method with 2 subteps
	// Case   : parabola /\
	// iters   = 16
	// error  : tolerance displacements 1.435e-01
	//
	// Solver : NR method with 10 subteps
	// Case   : parabola /\
	// iters   = 47
	// error  : tolerance displacements 4.439e-02
	//
	// Solver : NR method with 1000 subteps
	// Case   : parabola /\
	// iters   = 3003
	// Do      = [7.50389]
	// De      = [7.50000]
	// norm(D) = 5.19e-04
	// Fo      = [168.75000]
	// Fe      = [168.75000]
	// norm(F) = 0.00e+00
	//
	// Solver : NR method with 100000 subteps
	// Case   : parabola /\
	// iters   = 200002
	// Do      = [7.50004]
	// De      = [7.50000]
	// norm(D) = 5.20e-06
	// Fo      = [168.75000]
	// Fe      = [168.75000]
	// norm(F) = 0.00e+00
	//
	// Solver : Arc method
	// Case   : parabola /\
	// iters   = 16893
	// Do      = [7.50022]
	// De      = [7.50000]
	// norm(D) = 2.98e-05
	// Fo      = [168.75681]
	// Fe      = [168.75000]
	// norm(F) = 4.04e-05
	//
	// Solver : The Newton-Raphson method(Load control)
	// Case   : parabola /\ near top
	// iters   = 5
	// error  : tolerance displacements 7.028e-01
	//
	// Solver : NR method with 2 subteps
	// Case   : parabola /\ near top
	// iters   = 17
	// error  : tolerance displacements 8.748e-01
	//
	// Solver : NR method with 10 subteps
	// Case   : parabola /\ near top
	// iters   = 50
	// error  : tolerance displacements 9.846e-02
	//
	// Solver : NR method with 1000 subteps
	// Case   : parabola /\ near top
	// iters   = 3003
	// Do      = [10.01098]
	// De      = [10.00000]
	// norm(D) = 1.10e-03
	// Fo      = [200.00000]
	// Fe      = [200.00000]
	// norm(F) = 0.00e+00
	//
	// Solver : NR method with 100000 subteps
	// Case   : parabola /\ near top
	// iters   = 200002
	// Do      = [10.00011]
	// De      = [10.00000]
	// norm(D) = 1.10e-05
	// Fo      = [200.00000]
	// Fe      = [200.00000]
	// norm(F) = 0.00e+00
	//
	// Solver : Arc method
	// Case   : parabola /\ near top
	// iters   = 20028
	// Do      = [10.00002]
	// De      = [10.00000]
	// norm(D) = 1.52e-06
	// Fo      = [200.00563]
	// Fe      = [200.00000]
	// norm(F) = 2.82e-05
	//
	// Solver : The Newton-Raphson method(Load control)
	// Case   : parabola /\ after top
	// iters   = 5
	// error  : tolerance displacements 9.237e-01
	//
	// Solver : NR method with 2 subteps
	// Case   : parabola /\ after top
	// iters   = 12
	// error  : tolerance displacements 9.268e-01
	//
	// Solver : NR method with 10 subteps
	// Case   : parabola /\ after top
	// iters   = 44
	// error  : tolerance displacements 9.281e-01
	//
	// Solver : NR method with 1000 subteps
	// Case   : parabola /\ after top
	// iters   = 3003
	// error  : tolerance displacements 9.286e-01
	//
	// Solver : NR method with 100000 subteps
	// Case   : parabola /\ after top
	// iters   = 200002
	// error  : tolerance displacements 9.286e-01
	//
	// Solver : Arc method
	// Case   : parabola /\ after top
	// iters   = 39627
	// Do      = [28.00030]
	// De      = [28.00000]
	// norm(D) = 1.07e-05
	// Fo      = [56.03242]
	// Fe      = [56.00000]
	// norm(F) = 5.79e-04
	//
	// Solver : The Newton-Raphson method(Load control)
	// Case   : curve _/ before top
	// iters   = 1
	// error  : tolerance displacements 6.667e-01
	//
	// Solver : NR method with 2 subteps
	// Case   : curve _/ before top
	// iters   = 13
	// error  : tolerance displacements 9.928e-02
	//
	// Solver : NR method with 10 subteps
	// Case   : curve _/ before top
	// error  : Too big force 5e+09. Last coordinate: D=[-2.951e+03],F=[4.654e+09]
	//
	// Solver : NR method with 1000 subteps
	// Case   : curve _/ before top
	// iters   = 3003
	// error  : tolerance displacements 7.216e-01
	//
	// Solver : NR method with 100000 subteps
	// Case   : curve _/ before top
	// iters   = 200002
	// error  : tolerance displacements 7.215e-01
	//
	// Solver : Arc method
	// Case   : curve _/ before top
	// iters   = 9061
	// error  : tolerance forces 8.182e-01
	//
	// Solver : The Newton-Raphson method(Load control)
	// Case   : curve _/ after top
	// iters   = 3
	// error  : tolerance displacements 7.148e-01
	//
	// Solver : NR method with 2 subteps
	// Case   : curve _/ after top
	// error  : Too big displacement -2e+06. Last coordinate: D=[-1.543e+06],F=[6.612e+17]
	//
	// Solver : NR method with 10 subteps
	// Case   : curve _/ after top
	// iters   = 51
	// error  : tolerance displacements 7.581e-01
	//
	// Solver : NR method with 1000 subteps
	// Case   : curve _/ after top
	// iters   = 3009
	// error  : tolerance displacements 7.383e-01
	//
	// Solver : NR method with 100000 subteps
	// Case   : curve _/ after top
	// iters   = 201168
	// error  : tolerance displacements 7.382e-01
	//
	// Solver : Arc method
	// Case   : curve _/ after top
	// iters   = 16774
	// error  : tolerance forces 3.336e-01
	//
	// Solver : The Newton-Raphson method(Load control)
	// Case   : 2 unknown curve
	// iters   = 2
	// error  : tolerance displacements 8.287e-01
	//
	// Solver : NR method with 2 subteps
	// Case   : 2 unknown curve
	// error  : Too big force -3e+17. Last coordinate: D=[-8.477e+05 -6.093e+05],F=[-2.715e+17 -7.310e+17]
	//
	// Solver : NR method with 10 subteps
	// Case   : 2 unknown curve
	// error  : Too big force -1e+07. Last coordinate: D=[-3.055e+02 -1.972e+02],F=[-9.700e+06 -3.512e+07]
	//
	// Solver : NR method with 1000 subteps
	// Case   : 2 unknown curve
	// error  : Too big force -5e+06. Last coordinate: D=[-2.093e+02 -1.553e+02],F=[-4.793e+06 -1.140e+07]
	//
	// Solver : NR method with 100000 subteps
	// Case   : 2 unknown curve
	// iters   = 204415
	// error  : tolerance displacements 1.910e-02
	//
	// Solver : Arc method
	// Case   : 2 unknown curve
	// iters   = 1374
	// Do      = [2.34856 1.59900]
	// De      = [2.34076 1.59310]
	// norm(D) = 3.44e-03
	// Fo      = [12.34909 4.63091]
	// Fe      = [12.34160 4.62810]
	// norm(F) = 6.07e-04
}
