package hd

import (
	"encoding/json"
	"testing"
)

func TestJsonModel(t *testing.T) {
	m := Model{
		Points: [][2]float64{
			[2]float64{0.0, 0.0},
			[2]float64{1.0, 0.0},
		},
		Beams: []BeamProp{
			{
				N: [2]int{0, 1},
				A: 12e-4,
				J: 120e-6,
				E: 2.0e11,
			},
		},
		Supports: [][3]bool{
			[3]bool{true, true, true},
			[3]bool{false, false, false},
		},
		LoadCases: []LoadCase{
			LoadCase{
				LoadNodes: []LoadNode{
					{
						N:      1,
						Forces: [3]float64{0, 2.3, 0},
					},
					{
						N:      1,
						Forces: [3]float64{10, 0, 0},
					},
				},
			},
		},
	}
	b, err := json.Marshal(m)
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("%s", b)

	var mu Model
	err = json.Unmarshal(b, &mu)
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("%#v", mu)

	err = mu.Run(nil)
	if err != nil {
		t.Fatal(err)
	}
}
