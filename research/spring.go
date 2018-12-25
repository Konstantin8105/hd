// +build ignore
package main

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"math"
	"os"

	"github.com/Konstantin8105/hd"
)

func baseSpring(size int) hd.Model {
	var m hd.Model

	for i := 0; i < size; i++ {
		x := float64(i)
		y := math.Sin(float64(i)) * 100.0
		m.Points = append(m.Points, [2]float64{x, y})
	}

	m.Points = append(m.Points, [2]float64{0, 100})

	for i := 0; i < size; i++ {
		if i == 0 {
			continue
		}
		s := i - 1
		f := i
		m.Beams = append(m.Beams, hd.BeamProp{
			N: [2]int{s, f},
			A: 12e-4,
			J: 120e-6,
			E: 2.0e11,
		})
		if i < 2 {
			continue
		}
		m.Beams = append(m.Beams, hd.BeamProp{
			N: [2]int{size, f},
			A: 12e-4,
			J: 120e-6,
			E: 2.0e11,
		})
	}

	for i := range m.Points {
		if i == 0 || i == len(m.Points)-1 {
			m.Supports = append(m.Supports, [3]bool{true, true, true})
			continue
		}
		m.Supports = append(m.Supports, [3]bool{false, false, false})
	}

	m.LoadCases = []hd.LoadCase{
		{
			LoadNodes: []hd.LoadNode{
				{N: size - 1, Forces: [3]float64{0, 2.3, 0}},
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

	resultFile := "spring.json"

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
	m := baseSpring(a)

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
