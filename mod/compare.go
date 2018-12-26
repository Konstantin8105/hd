package mod

import (
	"fmt"
	"math"

	"github.com/Konstantin8105/errors"
	"github.com/Konstantin8105/hd"
)

// Compare reactions and displacement, return nil is load cases are same or
// return error for other cases.
func Compare(exp, act []hd.LoadCase, eps float64) error {
	et := errors.New("comparing load cases")

	if len(exp) != len(act) {
		_ = et.Add(fmt.Errorf("length of load cases is not same : %d != %d", len(exp), len(act)))
	}
	if eps <= 0.0 {
		_ = et.Add(fmt.Errorf("epsilon is not valid: %14e", eps))
	}

	if et.IsError() {
		return et
	}

	// compare reactions
	for e := range exp {
		for r := range exp[e].Reactions {
			expReact := exp[e].Reactions[r]
			actReact := act[e].Reactions[r]
			for j := 0; j < 3; j++ {
				if expReact[j] == 0.0 {
					if math.Abs(actReact[j]) > eps {
						_ = et.Add(fmt.Errorf("expReact is zero, actReact != 0"))
					}
					continue
				}
				diff := math.Abs((expReact[j] - actReact[j]) / expReact[j])
				if diff > eps {
					_ = et.Add(fmt.Errorf("reactions    `%d` is not same: %15e != %15e", j, expReact[j], actReact[j]))
				}
			}
		}
	}

	// compare displacement
	for e := range exp {
		for d := range exp[e].PointDisplacementGlobal {
			expDisp := exp[e].PointDisplacementGlobal[d]
			actDisp := act[e].PointDisplacementGlobal[d]
			for j := 0; j < 3; j++ {
				if expDisp[j] == 0.0 {
					if math.Abs(actDisp[j]) > eps {
						_ = et.Add(fmt.Errorf("expDisp is zero, actDisp != 0\n"+
							"actDisp = %10.8e", actDisp[j]))
					}
					continue
				}
				diff := math.Abs((expDisp[j] - actDisp[j]) / expDisp[j])
				if diff > eps {
					_ = et.Add(fmt.Errorf("displacement `%d` is not same: %15e != %15e", j, expDisp[j], actDisp[j]))
				}
			}
		}
	}

	if et.IsError() {
		return et
	}

	return nil
}
