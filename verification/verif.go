package verif

import (
	"github.com/Konstantin8105/hd"
)

func G() (model hd.Model, lc hd.LoadCase, name string, isOk func(lc *hd.LoadCase) (tol []float64)) {
	name = `
Book:
William McGuire, Richard H.Gallagher, Ronald D.Ziemian
Matrix Structural Analysis

EXAMPLE 9.1
Page 247
`
	model = hd.Model{
		Points: [][2]float64{
			{0.0, 4.0}, // a
			{4.0, 4.0}, // b
			{4.0, 0.0}, // c
		},
		Beams: []hd.BeamProp{
			{N: [2]int{0, 1}, A: 2e-6, J: 1, E: 200000e6},    // ab
			{N: [2]int{1, 2}, A: 5000e-6, J: 1, E: 200000e6}, // bc
		},
		Pins: [][6]bool{
			{false, true, true, false, true, true},
			{false, true, true, false, true, true},
		},
		Supports: [][3]bool{
			{true, true, false},
			{false, false, false},
			{true, true, false},
		},
	}
	P := 500e3
	a := 0.05
	lc = hd.LoadCase{
		LoadNodes: []hd.LoadNode{
			{N: 1, Forces: [3]float64{a * P, -P, 0}},
		},
	}
	lc.NonlinearNR.MaxIterations = 29000
	lc.NonlinearNR.Substep = 50

	lc.NonlinearNK.MaxIterations = 29000
	lc.NonlinearNK.Substep = 50.0

	// TODO: add arc-method

	isOk = func(lc *hd.LoadCase) (tol []float64) {
		var (
			expectD = 1692.0 / 1000.0
			expectP = 340.0 * 1000.0
		)
		for i := range lc.NonlinearNK.Results {
			if i == 0 {
				continue
			}
			lastD := lc.NonlinearNK.Results[i-1].PointDisplacementGlobal[1][0]
			presD := lc.NonlinearNK.Results[i].PointDisplacementGlobal[1][0]
			if lastD <= expectD && expectD <= presD {
				lastP := lc.NonlinearNK.Results[i-1].Reactions[2][1]
				presP := lc.NonlinearNK.Results[i].Reactions[2][1]
				p := lastP + (presP-lastP)*(expectD-lastD)/(presD-lastD)
				tol = append(tol, (p-expectP)/expectP*100.0)
			}
		}

		for i := range lc.NonlinearNR.Results {
			if i == 0 {
				continue
			}
			lastD := lc.NonlinearNR.Results[i-1].PointDisplacementGlobal[1][0]
			presD := lc.NonlinearNR.Results[i].PointDisplacementGlobal[1][0]
			if lastD <= expectD && expectD <= presD {
				lastP := lc.NonlinearNR.Results[i-1].Reactions[2][1]
				presP := lc.NonlinearNR.Results[i].Reactions[2][1]
				p := lastP + (presP-lastP)*(expectD-lastD)/(presD-lastD)
				tol = append(tol, (p-expectP)/expectP*100.0)
			}
		}
		return
	}
	return
}
