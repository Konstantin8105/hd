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
		Beams: [][2]int{
			[2]int{0, 1},
		},
		Supports: [][3]bool{
			[3]bool{true, true, true},
			[3]bool{false, false, false},
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
	t.Logf("%+v", mu)
}
