// +build ignore
package main

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"os"

	"github.com/Konstantin8105/hd"
)

// frame design:
//	0--0--0--0--0--0--0--0--0
//	|\ |\ |\ |\ |\ |\ |\ |\ |
//	| \| \| \| \| \| \| \| \|
//	0--0--0--0--0--0--0--0--0
func baseFrame(size int) hd.Model {
	// even size points
	size = size + size%2

	// generate frame
	var m hd.Model

	rows := size / 2

	for i := 0; i < rows; i++ {
		m.Points = append(m.Points, [2]float64{float64(i) * 10.0, 0.0})
		m.Points = append(m.Points, [2]float64{float64(i) * 10.0, 10.0})
		if i > 0 {
			m.Beams = append(m.Beams, hd.BeamProp{
				N: [2]int{2 * i, 2 * (i - 1)},
				A: 12e-4, J: 120e-6, E: 2.0e11,
			})
			m.Beams = append(m.Beams, hd.BeamProp{
				N: [2]int{2*i + 1, 2*(i-1) + 1},
				A: 12e-4, J: 120e-6, E: 2.0e11,
			})
			m.Beams = append(m.Beams, hd.BeamProp{
				N: [2]int{2 * i, 2*i + 1},
				A: 12e-4, J: 120e-6, E: 2.0e11,
			})
			m.Beams = append(m.Beams, hd.BeamProp{
				N: [2]int{2*(i-1) + 1, 2 * i},
				A: 12e-4, J: 120e-6, E: 2.0e11,
			})
		}
	}

	for i := range m.Points {
		if i < 2 {
			m.Supports = append(m.Supports, [3]bool{true, true, true})
			continue
		}
		m.Supports = append(m.Supports, [3]bool{false, false, false})
	}

	m.LoadCases = []hd.LoadCase{
		{
			LoadNodes: []hd.LoadNode{
				{N: size - 1, Forces: [3]float64{100.0, 250.0, 0}},
			},
		},
	}

	return m
}

type result struct {
	AmountIntermediantPoints int
	Displacement             string
}

// Main points of software:
// * example of using Golang package hd
// * create graph precision by amount of separation points
func main() {

	resultFile := "frame.json"

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
	a := 10
	for _, r := range rs {
		if r.AmountIntermediantPoints >= a {
			a = r.AmountIntermediantPoints * 2
		}
	}
	fmt.Fprintf(os.Stdout, "Calculate with %d intermediant points\n", a)
	m := baseFrame(a)

	// run model calculation
	err = m.Run(nil)
	if err != nil {
		panic(err)
	}

	// write result
	pdg := m.LoadCases[0].PointDisplacementGlobal
	rs = append(rs, result{
		AmountIntermediantPoints: a,
		Displacement:             fmt.Sprintf("%20.19e", pdg[len(pdg)-1][1]),
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
