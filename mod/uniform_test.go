package mod

import (
	"bytes"
	"fmt"
	"math"
	"testing"

	"github.com/Konstantin8105/hd"
)

func TestLoadUniform(t *testing.T) {
	t.Run("error checking", func(t *testing.T) {
		// error checking
		errs := []error{
			func() error {
				var m *hd.Model
				_, err := LoadUniform(m, 0, false, [2]float64{0, 0})
				return err
			}(),
			func() error {
				m := hd.Model{}
				_, err := LoadUniform(&m, 0, false, [2]float64{0, 0})
				return err
			}(),
			func() error {
				m := hd.Model{}
				_, err := LoadUniform(&m, -1, false, [2]float64{0, 0})
				return err
			}(),
			func() error {
				m := hd.Model{}
				_, err := LoadUniform(&m, 100, false, [2]float64{0, 0})
				return err
			}(),
			func() error {
				m := hd.Model{
					Beams: []hd.BeamProp{
						{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
					},
				}
				_, err := LoadUniform(&m, 1, false, [2]float64{0, 0})
				return err
			}(),
			func() error {
				m := hd.Model{
					Beams: []hd.BeamProp{
						{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
					},
				}
				_, err := LoadUniform(&m, 0, false, [2]float64{math.Inf(0), math.NaN()})
				return err
			}(),
		}

		for i := range errs {
			if errs[i] == nil {
				t.Errorf("case %d is not correct", i)
			}
			t.Log(errs[i])
		}
	})

	// calculation checking
	t.Run("check on beam with part load", func(t *testing.T) {
		for _, proj := range []bool{false, true} {
			for _, mirror := range []bool{false, true} {
				for _, d := range []float64{200, -0.5, 0.5, -200, 0} {
					size := 3
					t.Run(fmt.Sprintf("Mirror:%v/Proj:%v/%3.1f/Size:%d", mirror, proj, d, size), func(t *testing.T) {
						m := hd.Model{
							Points: [][2]float64{
								{0.0, 0.0},
								{1.0, d},
								{2.0, 0.0},
							},
							Beams: []hd.BeamProp{
								{N: [2]int{0, 1}, A: 12e-4, J: 120e-6, E: 2.0e11},
								{N: [2]int{1, 2}, A: 12e-4, J: 120e-6, E: 2.0e11},
							},
							Supports: [][3]bool{
								{true, true, true},
								{false, false, false},
								{false, true, false},
							},
							LoadCases: []hd.LoadCase{
								{LoadNodes: []hd.LoadNode{{N: 1, Forces: [3]float64{100.0, 0.0, 0.0}}}},
							},
						}

						if mirror {
							m.Points[0][0] = 2.0
							m.Points[2][0] = 0.0
						}

						if err := SplitBeam(&m, 0, size); err != nil {
							t.Fatal(err)
						}
						if err := SplitBeam(&m, 1, size); err != nil {
							t.Fatal(err)
						}

						// reset and allocate memory
						m.LoadCases = make([]hd.LoadCase, 1)
						// loads
						ux := -10.0
						uy := -25.0
						// uniform load
						for i := range m.Beams {
							if m.Points[m.Beams[i].N[0]][0] > 1.0 || m.Points[m.Beams[i].N[1]][0] > 1.0 {
								// not add on right part of beam
								continue
							}
							un, err := LoadUniform(&m, i, proj, [2]float64{ux, uy})
							if err != nil {
								t.Fatal(err)
							}
							m.LoadCases[0].LoadNodes = append(m.LoadCases[0].LoadNodes, un...)
						}
						// node load
						mn := m
						mn.LoadCases = make([]hd.LoadCase, 1)
						dx := math.Abs(m.Points[0][0] - m.Points[1][0])
						dy := math.Abs(m.Points[0][1] - m.Points[1][1])
						if !proj {
							// not projection
							dx = math.Sqrt(dx*dx + dy*dy)
							dy = dx
						}
						for i := range mn.Points {
							if m.Points[i][0] > 1.0 {
								// not add on right part of beam
								continue
							}
							if m.Points[i][0] == 0.0 || m.Points[i][0] == 1.0 {
								mn.LoadCases[0].LoadNodes = append(mn.LoadCases[0].LoadNodes, hd.LoadNode{
									N: i,
									Forces: [3]float64{
										ux * dy / float64(size+1) / 2.0, // X
										uy * dx / float64(size+1) / 2.0, // Y
										0.0,                             // M
									},
								})
								continue
							}
							mn.LoadCases[0].LoadNodes = append(mn.LoadCases[0].LoadNodes, hd.LoadNode{
								N: i,
								Forces: [3]float64{
									ux * dy / float64(size+1), // X
									uy * dx / float64(size+1), // Y
									0.0,                       // M
								},
							})
						}

						// calculation
						var buf bytes.Buffer
						if err := m.Run(&buf); err != nil {
							t.Fatalf("m model error : %v", err)
						}
						buf.Reset()
						if err := mn.Run(&buf); err != nil {
							t.Fatalf("mn model error : %v", err)
						}
						buf.Reset()

						// comparing displacement
						eps := 0.05 // 5%
						var actual, sum float64
						for i := 0; i < 3; i++ {
							actual += math.Pow(
								m.LoadCases[0].PointDisplacementGlobal[1][i]-
									mn.LoadCases[0].PointDisplacementGlobal[1][i],
								2)
							sum += math.Pow(m.LoadCases[0].PointDisplacementGlobal[1][i], 2)
						}
						actual = math.Sqrt(actual)
						sum = math.Sqrt(sum)
						if actual/sum >= eps {
							t.Log(actual / sum)
							t.Log(m.LoadCases[0].PointDisplacementGlobal[1])
							t.Log(mn.LoadCases[0].PointDisplacementGlobal[1])
							t.Errorf("displacement precision of calculation is not ok: %10.5f >= 0.05", actual/sum)
						}
						t.Logf("displacament precition : %14f <= %14f", actual/sum, eps)

						// comparing reactions
						actual = 0.0
						sum = 0.0
						for i := 0; i < 3; i++ {
							actual += math.Pow(
								m.LoadCases[0].Reactions[0][i]-
									mn.LoadCases[0].Reactions[0][i],
								2)
							sum += math.Pow(mn.LoadCases[0].Reactions[0][i], 2)
						}
						actual = math.Sqrt(actual)
						sum = math.Sqrt(sum)
						if actual/sum >= eps {
							t.Log(actual / sum)
							t.Log(m.LoadCases[0].Reactions[0])
							t.Log(mn.LoadCases[0].Reactions[0])
							t.Errorf("reaction     precision of calculation is not ok: %10.5f >= 0.05", actual/sum)
						}
						t.Logf("reaction     precition : %14f <= %14f", actual/sum, eps)
					})
				}
			}
		}
	})
}
