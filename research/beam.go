// +build ignore
package main

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"os"

	"github.com/Konstantin8105/hd"
)

func baseBeam() hd.Model {
	return hd.Model{
		Points: [][2]float64{
			{0.0, 0.0},
			{1.0, 0.0},
		},
		Beams: []hd.BeamProp{
			{
				N: [2]int{0, 1},
				A: 12e-4,
				J: 120e-6,
				E: 2.0e11,
			},
		},
		Supports: [][3]bool{
			{true, true, true},
			{false, false, false},
		},
		LoadCases: []hd.LoadCase{
			{
				LoadNodes: []hd.LoadNode{
					{N: 1, Forces: [3]float64{0, 2.3, 0}},
				},
			},
		},
	}
}

type result struct {
	AmountIntermediantPoints int
	Displacement             string
}

// Main points of software:
// * example of using Golang package hd
// * create graph precision by amount of separation points
func main() {

	resultFile := "beam.json"

	//
	// Only for begin create resultFile
	//
	if _, err := os.Stat(resultFile); os.IsNotExist(err) {
		var rs []result
		bo, err := json.Marshal(rs)
		if err != nil {
			panic(err)
		}
		err = ioutil.WriteFile(resultFile, bo, 0644)
		if err != nil {
			panic(err)
		}
	}

next:
	// get file with results in json format
	b, err := ioutil.ReadFile(resultFile)
	if err != nil {
		panic(err)
	}

	// read all data
	var rs []result
	err = json.Unmarshal(b, &rs)
	if err != nil {
		panic(err)
	}

	// show all results in stdout
	for _, r := range rs {
		fmt.Fprintf(os.Stdout, "%10d %30s\n", r.AmountIntermediantPoints, r.Displacement)
	}

	// calculate next amount separation points
	a := 1
	for _, r := range rs {
		if r.AmountIntermediantPoints >= a {
			a = r.AmountIntermediantPoints * 2
		}
	}
	fmt.Fprintf(os.Stdout, "Calculate with %d intermediant points\n", a)
	m := baseBeam()
	err = m.SplitBeam(0, a)
	if err != nil {
		panic(err)
	}

	// run model calculation
	err = m.Run(nil)
	if err != nil {
		panic(err)
	}

	// write result
	rs = append(rs, result{
		AmountIntermediantPoints: a,
		Displacement: fmt.Sprintf("%20.19e",
			m.LoadCases[0].PointDisplacementGlobal[1][1]),
	})

	bo, err := json.Marshal(rs)
	if err != nil {
		panic(err)
	}
	err = ioutil.WriteFile(resultFile, bo, 0644)
	if err != nil {
		panic(err)
	}

	// go to next iteration
	goto next
}
