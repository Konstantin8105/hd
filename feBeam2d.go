package hd

import (
	"math"

	"github.com/Konstantin8105/pow"
	"gonum.org/v1/gonum/mat"
)

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

// matrix of stiffiner for beam 2d.
// local linear-elastic stiffiness matrix
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
		kr.Scale(axialForce, kr)
	}()

	length := m.distance(m.Beams[pos].N[0], m.Beams[pos].N[1])

	// no pins
	if m.Pins[pos][2] == false && m.Pins[pos][5] == false {
		G1 := 6. / (5. * length)
		kr.Set(1, 1, +G1)
		kr.Set(4, 4, +G1)
		kr.Set(1, 4, -G1)
		kr.Set(4, 1, -G1)

		kr.Set(1, 2, +0.1)
		kr.Set(2, 1, +0.1)
		kr.Set(1, 5, +0.1)
		kr.Set(5, 1, +0.1)

		kr.Set(2, 4, -0.1)
		kr.Set(4, 2, -0.1)
		kr.Set(5, 4, -0.1)
		kr.Set(4, 5, -0.1)

		G2 := 2. * length / 15.
		kr.Set(2, 2, +G2)
		kr.Set(5, 5, +G2)

		G3 := -length / 30.
		kr.Set(5, 2, +G3)
		kr.Set(2, 5, +G3)

		return kr
	}

	// pins at the begin of beam
	if m.Pins[pos][2] && m.Pins[pos][5] == false {
		G1 := 6. / (5. * length)
		kr.Set(1, 1, +G1)
		kr.Set(4, 4, +G1)
		kr.Set(1, 4, -G1)
		kr.Set(4, 1, -G1)

		kr.Set(1, 5, +1./5.)
		kr.Set(5, 1, +1./5.)

		kr.Set(4, 5, -1./5.)
		kr.Set(5, 4, -1./5.)

		kr.Set(5, 5, length/5.0)

		return kr
	}

	// pins at the end of beam
	if m.Pins[pos][2] == false && m.Pins[pos][5] {
		G1 := 6. / (5. * length)
		kr.Set(1, 1, +G1)
		kr.Set(4, 4, +G1)
		kr.Set(1, 4, -G1)
		kr.Set(4, 1, -G1)

		kr.Set(1, 2, +1./5.)
		kr.Set(2, 1, +1./5.)

		kr.Set(4, 2, -1./5.)
		kr.Set(2, 4, -1./5.)

		kr.Set(2, 2, length/5.0)

		return kr
	}

	// pins on the start and end of beam
	if m.Pins[pos][2] && m.Pins[pos][5] {
		kr.Set(2, 2, +1.0/length)
		kr.Set(5, 5, +1.0/length)
		kr.Set(2, 5, -1.0/length)
		kr.Set(5, 2, -1.0/length)

		return kr
	}

	// TODO: add implementation
	panic("add implementation")
}

// TODO: The tangent stiffness matrix contains information regarding both the linear elastic and geometric stiffness of the structure. 
