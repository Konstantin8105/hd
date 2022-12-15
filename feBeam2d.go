package hd

import (
	"fmt"
	"math"

	"github.com/Konstantin8105/pow"
	"gonum.org/v1/gonum/mat"
)

// Ke - elastic matrix stiffiner
// Kg - geometric stiffness matrix
// second-order elastic analysis: [Ke+Kg]*D=P
// Km - plastic reduction matrix
// first-order inelestic analysis: [Ke+Km]*D=P
// second-order inelestic analysis: [Ke+Kg+Km]*D=P

// distance between 2 point of beam
func (m Model) distance(st, en int) (res float64) {
	defer func() {
		if r := recover(); r != nil {
			res = math.NaN()
		}
	}()
	return math.Hypot(m.Points[st][0]-m.Points[en][0],
		m.Points[st][1]-m.Points[en][1])
}

// matrix of stiffiner for beam 2d
func (m Model) getStiffBeam2d(pos int) *mat.Dense {
	data := make([]float64, 36)
	kr := mat.NewDense(6, 6, data)

	length := m.distance(m.Beams[pos].N[0], m.Beams[pos].N[1])

	EFL := m.Beams[pos].E * m.Beams[pos].A / length
	kr.Set(0, 0, +EFL)
	kr.Set(0, 3, -EFL)
	kr.Set(3, 0, -EFL)
	kr.Set(3, 3, +EFL)

	EJL3 := 12.0 * m.Beams[pos].E * m.Beams[pos].J / pow.E3(length)
	kr.Set(1, 1, +EJL3)
	kr.Set(4, 4, +EJL3)
	kr.Set(1, 4, -EJL3)
	kr.Set(4, 1, -EJL3)

	EJL2 := 6.0 * m.Beams[pos].E * m.Beams[pos].J / pow.E2(length)
	kr.Set(1, 2, +EJL2)
	kr.Set(2, 1, +EJL2)
	kr.Set(1, 5, +EJL2)
	kr.Set(5, 1, +EJL2)
	kr.Set(2, 4, -EJL2)
	kr.Set(4, 2, -EJL2)
	kr.Set(4, 5, -EJL2)
	kr.Set(5, 4, -EJL2)

	EJL := 2.0 * m.Beams[pos].E * m.Beams[pos].J / length
	kr.Set(2, 5, EJL)
	kr.Set(5, 2, EJL)

	EJL = 2 * EJL
	kr.Set(2, 2, EJL)
	kr.Set(5, 5, EJL)

	for fr := 0; fr < 6; fr++ {
		if !m.Pins[pos][fr] {
			// DoF is rigid
			continue
		}
		// DoF is free
		// Algorithm:
		// from
		// [ KAA , KAB ]
		// [ KBA , KBB ]
		// to
		// [ 0   , 0   ]
		// [ 0   , KS  ]
		//
		// where KAA is stiffiner coefficient with free DoF
		//
		// KS = [ KBB ] - [ KBA ] * [ KAA ]^(-1) * [ KAB ]
		dataKBB := make([]float64, 5*5)
		KBB := mat.NewDense(5, 5, dataKBB)

		dataKBB2 := make([]float64, 5*5)
		KBB2 := mat.NewDense(5, 5, dataKBB2)

		dataKS := make([]float64, 5*5)
		KS := mat.NewDense(5, 5, dataKS)

		dataKBA := make([]float64, 5)
		KBA := mat.NewDense(5, 1, dataKBA)

		dataKAB := make([]float64, 5)
		KAB := mat.NewDense(1, 5, dataKAB)

		// initialization data
		KAA := kr.At(fr, fr)
		if KAA == 0.0 {
			uncorrect := false
			for i := 0; i < 6; i++ {
				uncorrect = uncorrect || kr.At(i, fr) != 0.0 || kr.At(fr, i) != 0.0
			}
			if uncorrect {
				panic("KAA")
			}
			continue
		}
		for i := 0; i < 6; i++ {
			ia := i
			if i == fr {
				continue
			} else if i > fr {
				ia--
			}
			KBA.Set(ia, 0, kr.At(i, fr))
			KAB.Set(0, ia, kr.At(fr, i))
			for j := 0; j < 6; j++ {
				ja := j
				if j == fr {
					continue
				} else if j > fr {
					ja--
				}
				KBB.Set(ia, ja, kr.At(i, j))
			}
		}
		// calculate KS
		KBB2.Mul(KBA, KAB)
		KBB2.Scale(1/KAA, KBB2)
		KS.Sub(KBB, KBB2)

		for i := 0; i < 5; i++ {
			ia := i
			if i == fr {
				continue
			} else if i > fr {
				ia++
			}
			for j := 0; j < 5; j++ {
				ja := j
				if j == fr {
					continue
				} else if j > fr {
					ja++
				}
				kr.Set(ia, ja, KS.At(i, j))
			}
		}
		for i := 0; i < 6; i++ {
			kr.Set(i, fr, 0.0)
			kr.Set(fr, i, 0.0)
		}
	}

	return kr
}

func (m Model) getCoordTransStiffBeam2d(pos int) *mat.Dense {
	data := make([]float64, 36)
	tr := mat.NewDense(6, 6, data)

	length := m.distance(m.Beams[pos].N[0], m.Beams[pos].N[1])

	start := m.Beams[pos].N[0]
	end := m.Beams[pos].N[1]

	lambdaXX := (m.Points[end][0] - m.Points[start][0]) / length
	lambdaXY := (m.Points[end][1] - m.Points[start][1]) / length
	lambdaYX := -lambdaXY
	lambdaYY := +lambdaXX

	tr.Set(0, 0, lambdaXX)
	tr.Set(3, 3, lambdaXX)

	tr.Set(0, 1, lambdaXY)
	tr.Set(3, 4, lambdaXY)

	tr.Set(1, 0, lambdaYX)
	tr.Set(4, 3, lambdaYX)

	tr.Set(1, 1, lambdaYY)
	tr.Set(4, 4, lambdaYY)

	tr.Set(2, 2, 1.0)
	tr.Set(5, 5, 1.0)

	return tr
}

// matrix of geometric stiffner for beam 2d
// TODO: it is look like P-DELTA analysis
func (m Model) getGeometricBeam2d(pos int, lc *LoadCase) *mat.Dense {
	data := make([]float64, 36)
	kr := mat.NewDense(6, 6, data)

	// compress axial force
	var axialForce float64
	if lc.BeamForces[pos][0] > 0 {
		// compress at the begin point of beam
		axialForce = lc.BeamForces[pos][0]
	}
	if -lc.BeamForces[pos][3] > axialForce {
		// compress at the end point of beam
		axialForce = -lc.BeamForces[pos][3]
	}

	defer func() {
		//TODO : propably : minus axial force ??
		// if axialForce < 0 {
		// 	axialForce = 0
		// }

		kr.Scale(axialForce, kr)
	}()

	length := m.distance(m.Beams[pos].N[0], m.Beams[pos].N[1])

	var found bool

	// number conversion:
	// 0 1 2 3 4 5  for full size finite element
	//   0 1   2 3  for   bending finite element
	for _, v := range []struct {
		pins [4]bool
		f    func()
	}{
		{
			pins: [4]bool{false, false, false, false},
			f: func() {
				// -0.03333 * (l * N) : [ 1, 3] :
				// -0.03333 * (l * N) : [ 3, 1] :
				//
				//  0.09990 * N       : [ 0, 1] :
				//  0.09990 * N       : [ 1, 0] :
				//  0.09996 * N       : [ 0, 3] :
				//  0.09996 * N       : [ 3, 0] :
				// -0.09990 * N       : [ 1, 2] :
				// -0.09990 * N       : [ 2, 1] :
				// -0.09996 * N       : [ 2, 3] :
				// -0.09996 * N       : [ 3, 2] :
				//
				//  0.13333 * (l * N) : [ 1, 1] :
				//  0.13333 * (l * N) : [ 3, 3] :
				//
				//  1.19988 * N / l   : [ 0, 0] :
				//  1.19988 * N / l   : [ 2, 2] :
				// -1.19988 * N / l   : [ 0, 2] :
				// -1.19988 * N / l   : [ 2, 0] :

				G := length / 30.
				kr.Set(5, 2, -G)
				kr.Set(2, 5, -G)

				G = 0.1
				kr.Set(1, 2, +G)
				kr.Set(2, 1, +G)
				kr.Set(1, 5, +G)
				kr.Set(5, 1, +G)
				kr.Set(2, 4, -G)
				kr.Set(4, 2, -G)
				kr.Set(5, 4, -G)
				kr.Set(4, 5, -G)

				G = 2.0 / 15.0 * length
				kr.Set(2, 2, +G)
				kr.Set(5, 5, +G)

				G = 6.0 / 5.0 * length
				kr.Set(1, 4, -G)
				kr.Set(4, 1, -G)
				kr.Set(1, 1, +G)
				kr.Set(4, 4, +G)

				//  TODO: is it oK?
				G = 1.0 / length
				kr.Set(0, 0, +G)
				kr.Set(3, 3, +G)
				kr.Set(3, 0, -G)
				kr.Set(0, 3, -G)
			},
		},
		{
			pins: [4]bool{true, false, false, false},
			f: func() {
				// -0.04165 * (l * N) : [ 1, 3] :
				// -0.04165 * (l * N) : [ 3, 1] :
				//
				// 0.12500 * (l * N)  : [ 3, 3] :
				// 0.12501 * (l * N)  : [ 1, 1] :
				G := 1.0 / 24.0 * length
				kr.Set(2, 5, -G)
				kr.Set(5, 2, -G)

				G = 1.0 / 8.0 * length
				kr.Set(2, 2, +G)
				kr.Set(5, 5, +G)
			},
		},
		{
			pins: [4]bool{false, true, false, false},
			f: func() {
				// -0.12493 * N       : [ 0, 3] :  -0.125
				// -0.12493 * N       : [ 3, 0] :  -0.125
				//  0.12493 * N       : [ 2, 3] :   0.125
				//  0.12493 * N       : [ 3, 2] :   0.125
				//
				// -1.12503 * (N / l) : [ 0, 2] :  -1.125 / l
				// -1.12503 * (N / l) : [ 2, 0] :  -1.125 / l
				//  1.12503 * (N / l) : [ 0, 0] :   1.125 / l
				//  1.12503 * (N / l) : [ 2, 2] :   1.125 / l
				//
				// 0.12500 * (l * N)  : [ 3, 3] :   0.125 * l
				//
				// 0.00000            : [ 0, 1] :
				// 0.00000            : [ 1, 0] :
				// 0.00000            : [ 1, 1] :
				// 0.00000            : [ 1, 2] :
				// 0.00000            : [ 1, 3] :
				// 0.00000            : [ 2, 1] :
				// 0.00000            : [ 3, 1] :

				// 0.12493 * N        : [ 0, 3] :
				// 0.12493 * N        : [ 3, 0] :
				// -0.12493 * N       : [ 2, 3] :
				// -0.12493 * N       : [ 3, 2] :
				//
				// 0.12500 * (l * N)  : [ 3, 3] :
				//
				// -1.12503 * (N / l) : [ 0, 2] :
				// -1.12503 * (N / l) : [ 2, 0] :
				// 1.12503 * (N / l)  : [ 0, 0] :
				// 1.12503 * (N / l)  : [ 2, 2] :

				// pins at the begin of beam
				G := 1.0 / 8.0
				kr.Set(1, 5, +G)
				kr.Set(5, 1, +G)
				kr.Set(4, 5, -G)
				kr.Set(5, 4, -G)

				G = 9.0 / 8.0 * 1.0 / length
				kr.Set(1, 4, -G)
				kr.Set(4, 1, -G)
				kr.Set(1, 1, +G)
				kr.Set(4, 4, +G)

				G = 1.0 / 8.0 * length
				kr.Set(5, 5, +G)
				// kr.Set(5, 5, length/5.0)
				// TODO: why it is different матрица потенциала
			},
		},
		{
			pins: [4]bool{false, false, true, false},
			f: func() {
				// -0.04165 * (l * N) : [ 1, 3] :
				// -0.04165 * (l * N) : [ 3, 1] :
				//
				// 0.12500 * (l * N)  : [ 3, 3] :
				// 0.12501 * (l * N)  : [ 1, 1] :
				G := 1.0 / 24.0 * length
				kr.Set(2, 5, -G)
				kr.Set(5, 2, -G)

				G = 1.0 / 8.0 * length
				kr.Set(5, 5, +G)
				kr.Set(2, 2, +G)
			},
		},
		{
			pins: [4]bool{false, false, false, true},
			f: func() {
				// 0.12489 * N        : [ 0, 1] :
				// 0.12489 * N        : [ 1, 0] :
				// -0.12489 * N       : [ 1, 2] :
				// -0.12489 * N       : [ 2, 1] :
				//
				// -1.12494 * (N / l) : [ 0, 2] :
				// -1.12494 * (N / l) : [ 2, 0] :
				// 1.12494 * (N / l)  : [ 0, 0] :
				// 1.12494 * (N / l)  : [ 2, 2] :
				//
				// 0.12500 * (l * N)  : [ 1, 1] :

				// pins at the end of beam
				G := 0.125
				kr.Set(1, 2, +G)
				kr.Set(2, 1, +G)
				kr.Set(2, 4, -G)
				kr.Set(4, 2, -G)
				// 	kr.Set(1, 2, +1./5.)
				// 	kr.Set(2, 1, +1./5.)
				// 	kr.Set(4, 2, -1./5.)
				// 	kr.Set(2, 4, -1./5.)

				G = 9.0 / 8.0 * 1.0 / length
				kr.Set(1, 4, -G)
				kr.Set(4, 1, -G)
				kr.Set(1, 1, +G)
				kr.Set(4, 4, +G)
				// 	G1 := 6. / (5. * length)
				// 	kr.Set(1, 1, +G1)
				// 	kr.Set(4, 4, +G1)
				// 	kr.Set(1, 4, -G1)
				// 	kr.Set(4, 1, -G1)

				G = 1.0 / 8.0 * length
				kr.Set(2, 2, +G)
				// 	kr.Set(2, 2, length/5.0)
			},
		},
		{
			pins: [4]bool{false, true, false, true},
			f: func() {
				// -1.00016 * (N / l) : [ 0, 2] :
				// -1.00016 * (N / l) : [ 2, 0] :
				// 1.00016 * (N / l)  : [ 0, 0] :
				// 1.00016 * (N / l)  : [ 2, 2] :

				// pins on the start and end of beam
				G := 1.0 / length
				kr.Set(0, 0, +G)
				kr.Set(1, 1, +G)
				kr.Set(3, 3, +G)
				kr.Set(4, 4, +G)
				kr.Set(0, 3, -G)
				kr.Set(3, 0, -G)
				kr.Set(1, 4, -G)
				kr.Set(4, 1, -G)
				// kr.Set(1, 4, -G)
				// kr.Set(4, 1, -G)
				// kr.Set(1, 1, +G)
				// kr.Set(4, 4, +G)
			},
		},
		{
			pins: [4]bool{true, true, false, false},
			f: func() {
				G := 10.0 / 9.0 * length
				kr.Set(5, 5, +G)
			},
		},
		{
			pins: [4]bool{false, false, true, true},
			f: func() {
				G := 10.0 / 9.0 * length
				kr.Set(2, 2, +G)
			},
		},
		// {
		// 	pins: [4]bool{true, true, true, true},
		// 	f: func() {
		// 		// truss
		// 		panic("not acceptable case")
		// 	},
		// },
	} {
		if (v.pins[0] == m.Pins[pos][1] &&
			v.pins[1] == m.Pins[pos][2] &&
			v.pins[2] == m.Pins[pos][4] &&
			v.pins[3] == m.Pins[pos][5]) ||
			// truss
			(m.Pins[pos][1] && m.Pins[pos][2] && m.Pins[pos][4] && m.Pins[pos][5] &&
				!v.pins[0] && v.pins[1] && !v.pins[2] && v.pins[3]) {
			v.f()
			found = true
			break
		}
	}
	if !found {
		panic(fmt.Errorf("not found geometric matrix case: %v", m.Pins[pos]))
	}

	// check symmetrical kr[6,6] -- data[36]
	for r := 0; r < 6; r++ {
		for c := 0; c < 6; c++ {
			v1, v2 := data[r*6+c], data[c*6+r]
			if v1 != v2 {
				panic(fmt.Errorf("[%d,%d]-{%f,%f}", r, c, v1, v2))
			}
		}
	}

	return kr
}
