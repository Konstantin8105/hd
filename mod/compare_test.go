package mod_test

import (
	"fmt"
	"testing"

	"github.com/Konstantin8105/hd"
	"github.com/Konstantin8105/hd/mod"
)

func TestCompare(t *testing.T) {
	tcs := []struct {
		exp, act []hd.LoadCase
		eps      float64
		isError  bool
	}{
		{
			exp:     nil,
			act:     nil,
			eps:     1.0,
			isError: false,
		},
		{
			exp:     nil,
			act:     nil,
			eps:     -1.0,
			isError: true,
		},
		{
			exp:     nil,
			act:     nil,
			eps:     0.0,
			isError: true,
		},
		{
			exp: []hd.LoadCase{
				hd.LoadCase{
					PointDisplacementGlobal: [][3]float64{
						[3]float64{0, 0, 0},
					},
					Reactions: [][3]float64{
						[3]float64{0, 0, 0},
					},
				},
			},
			act: []hd.LoadCase{
				hd.LoadCase{
					PointDisplacementGlobal: [][3]float64{
						[3]float64{1, 0, 0},
					},
					Reactions: [][3]float64{
						[3]float64{1, 0, 0},
					},
				},
			},
			eps:     0.001,
			isError: true,
		},
		{
			exp: []hd.LoadCase{
				hd.LoadCase{
					PointDisplacementGlobal: [][3]float64{
						[3]float64{1, 0, 0},
					},
					Reactions: [][3]float64{
						[3]float64{1, 0, 0},
					},
				},
			},
			act: []hd.LoadCase{
				hd.LoadCase{
					PointDisplacementGlobal: [][3]float64{
						[3]float64{2, 0, 0},
					},
					Reactions: [][3]float64{
						[3]float64{2, 0, 0},
					},
				},
			},
			eps:     0.001,
			isError: true,
		},
		{
			exp: []hd.LoadCase{
				hd.LoadCase{},
				hd.LoadCase{},
			},
			act: []hd.LoadCase{
				hd.LoadCase{},
			},
			eps:     0.001,
			isError: true,
		},
	}

	for i, tc := range tcs {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			err := mod.Compare(tc.exp, tc.act, tc.eps)
			if (err != nil) != tc.isError {
				t.Log(err)
				t.Fatalf("is not same")
			}
		})
	}
}
