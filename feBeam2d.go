package hd

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

// distance between 2 point of beam
func (m Model) distance(st, en int) (res float64) {
	defer func() {
		if r := recover(); r != nil {
			res = math.NaN()
		}
	}()
	var sum float64
	for i := 0; i < len(m.Points[st]); i++ {
		sum += math.Pow(m.Points[st][i]-m.Points[en][i], 2.0)
	}
	return math.Sqrt(sum)
}

func view(k *mat.Dense) {
	r, c := k.Dims()
	for j := 0; j < r; j++ {
		for i := 0; i < c; i++ {
			fmt.Printf(" %10.3e", k.At(j, i))
		}
		fmt.Printf("\n")
	}
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

	EJL3 := 12.0 * m.Beams[pos].E * m.Beams[pos].J / math.Pow(length, 3.0)
	kr.Set(1, 1, +EJL3)
	kr.Set(4, 4, +EJL3)
	kr.Set(1, 4, -EJL3)
	kr.Set(4, 1, -EJL3)

	EJL2 := 6.0 * m.Beams[pos].E * m.Beams[pos].J / (length * length)
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
