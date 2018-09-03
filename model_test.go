package hd

import (
	"bytes"
	"encoding/json"
	"testing"
)

func baseModel() Model {
	return Model{
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
}

func TestJsonModel(t *testing.T) {
	m := baseModel()
	b, err := json.Marshal(m)
	if err != nil {
		t.Fatal(err)
	}

	var mu Model
	err = json.Unmarshal(b, &mu)
	if err != nil {
		t.Fatal(err)
	}

	bu, err := json.Marshal(m)
	if err != nil {
		t.Fatal(err)
	}

	if !bytes.Equal(b, bu) {
		t.Fatalf("Is not same:\n%s\n%s", b, bu)
	}
}

func TestRun(t *testing.T) {
	m := baseModel()
	var b bytes.Buffer
	if err := m.Run(&b); err != nil {
		t.Fatalf("Error : %v", err)
	}
}
