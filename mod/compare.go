package mod

import (
	"fmt"
	"math"

	"github.com/Konstantin8105/hd"
)

func Compare(exp, act []hd.LoadCase, eps float64) error {

	// compare reactions
	for e := range exp {
		for r := range exp[e].Reactions {
			expReact := exp[e].Reactions[r]
			actReact := act[e].Reactions[r]
			for j := 0; j < 3; j++ {
				if expReact[j] == 0.0 {
					if math.Abs(actReact[j]) > eps {
						return fmt.Errorf("expReact is zero, actReact != 0")
					}
					continue
				}
				diff := math.Abs((expReact[j] - actReact[j]) / expReact[j])
				if diff > eps {
					return fmt.Errorf("Diff[%d] is not ok : %15.5e", j, diff)
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
						return fmt.Errorf("expDisp is zero, actDisp != 0\n"+
							"actDisp = %10.8e", actDisp[j])
					}
					continue
				}
				diff := math.Abs((expDisp[j] - actDisp[j]) / expDisp[j])
				if diff > eps {
					return fmt.Errorf("Diff[%d] is not ok : %15.5e", j, diff)
				}
			}
		}
	}

	return nil
}
