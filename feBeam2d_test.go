package hd

import (
	"fmt"
	"math"
	"testing"

	"github.com/Konstantin8105/errors"
)

func TestGeometricBeam2d(t *testing.T) {
	tcs := []struct {
		m      Model
		matrix [6][6]float64
	}{
		{
			m: Model{
				Points: [][2]float64{
					{0.0, 0.0},
					{1.0, 0.0},
				},
				Beams: []BeamProp{
					{N: [2]int{0, 1}, A: 40e-4, J: 1, E: 2.0e11},
				},
				Pins: [][6]bool{
					{false, false, false, false, false, false},
				},
				LoadCases: []LoadCase{
					{BeamForces: [][6]float64{{0, 0, 0, -1, 0, 0}}},
				},
			},
			matrix: [6][6]float64{
				{+0.0},
				{+0.0, +1.2},
				{+0.0, +0.1, +2. / 15.},
				{+0.0},
				{+0.0, -1.2, -0.1, +0.0, +1.2},
				{+0.0, +0.1, -1. / 30., 0.0, -0.1, +2. / 15.},
			},
		},
		{
			m: Model{
				Points: [][2]float64{
					{0.0, 0.0},
					{2.0, 0.0},
				},
				Beams: []BeamProp{
					{N: [2]int{0, 1}, A: 40e-4, J: 1, E: 2.0e11},
				},
				Pins: [][6]bool{
					{false, false, false, false, false, false},
				},
				LoadCases: []LoadCase{
					{BeamForces: [][6]float64{{1, 0, 0, 0, 0, 0}}},
				},
			},
			matrix: [6][6]float64{
				{+0.0},
				{+0.0, +1.2 / 2.},
				{+0.0, +0.1, +2. / 15. * 2.},
				{+0.0},
				{+0.0, -1.2 / 2., -0.1, +0.0, +1.2 / 2.},
				{+0.0, +0.1, -1. / 30. * 2., 0.0, -0.1, +2. / 15. * 2.},
			},
		},
		{
			m: Model{
				Points: [][2]float64{
					{0.0, 0.0},
					{2.0, 0.0},
				},
				Beams: []BeamProp{
					{N: [2]int{0, 1}, A: 40e-4, J: 1, E: 2.0e11},
				},
				Pins: [][6]bool{
					{false, false, true, false, false, false},
				},
				LoadCases: []LoadCase{
					{BeamForces: [][6]float64{{1, 0, 0, -1, 0, 0}}},
				},
			},
			matrix: [6][6]float64{
				{+0.0},
				{+0.0, +1.2 / 2.},
				{+0.0},
				{+0.0},
				{+0.0, -1.2 / 2., +0.0, +0.0, +1.2 / 2.},
				{+0.0, +0.2, +0.0, +0.0, -0.2, +0.2 * 2},
			},
		},
		{
			m: Model{
				Points: [][2]float64{
					{0.0, 0.0},
					{2.0, 0.0},
				},
				Beams: []BeamProp{
					{N: [2]int{0, 1}, A: 40e-4, J: 1, E: 2.0e11},
				},
				Pins: [][6]bool{
					{false, false, false, false, false, true},
				},
				LoadCases: []LoadCase{
					{BeamForces: [][6]float64{{1, 0, 0, -1, 0, 0}}},
				},
			},
			matrix: [6][6]float64{
				{+0.0},
				{+0.0, +1.2 / 2.},
				{+0.0, +0.2, 1. / 5. * 2.},
				{+0.0},
				{+0.0, -1.2 / 2., -0.2, 0, +1.2 / 2.},
				{+0.0},
			},
		},
	}

	for i := range tcs {
		for r := 0; r < 6; r++ {
			for c := 0; c < 6; c++ {
				if r > c {
					continue
				}
				tcs[i].matrix[r][c] = tcs[i].matrix[c][r]
			}
		}
	}

	for i := range tcs {
		t.Run(fmt.Sprintf("case %d", i), func(t *testing.T) {
			et := errors.New("comparing")

			actual := tcs[i].m.getGeometricBeam2d(0, &(tcs[i].m.LoadCases[0]))
			r, c := actual.Dims()
			if r != 6 || c != 6 {
				t.Fatalf("size of matrix is not ok. case %d", i)
			}

			// compare matrix
			expect := tcs[i].matrix

			eps := 1e-16
			for r := 0; r < 6; r++ {
				for c := 0; c < 6; c++ {
					if expect[r][c] == 0.0 && math.Abs(actual.At(r, c)) > eps {
						_ = et.Add(fmt.Errorf("case [%d,%d]: Zero checking. %+.5e",
							r, c, actual.At(r, c)))
						continue
					}
					if math.Abs(actual.At(r, c)-expect[r][c])/math.Abs(expect[r][c]) > eps {
						_ = et.Add(fmt.Errorf("case [%d,%d]: Value checking. %+.5e != %+.5e",
							r, c, actual.At(r, c), expect[r][c]))
						continue
					}
				}
			}

			if et.IsError() {
				t.Fatal(et)
			}
		})
	}
}
