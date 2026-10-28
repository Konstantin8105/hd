package main

import (
	"fmt"
	"os"

	"github.com/Konstantin8105/sm"
)

// run :
//	go run ./theory/axisymBeam2d

func cal(name, str string) string {
	fmt.Fprintf(os.Stdout, "\nName    : %s\n", name)
	fmt.Fprintf(os.Stdout, "\nFormula : %s\n", str)
	val, err := sm.Sexpr(
		// os.Stdout,
		nil,
		str)
	if err != nil {
		panic(err)
	}
	fmt.Fprintf(os.Stdout, "Value   : %s\n", val)
	return val
}

func main() {

	sm.MaxIteration = -1
	sm.FloatFormat = 5

	// -----------------------
	// CONICAL FINITE ELEMENT STRESS ANALYSIS OF AN AXISYMMETRIC
	// TANK ON AN ELASTIC FOUNDATION
	// by
	// HOWARD RICHARD HORN, JR., B.S. in C.E., M.S. In C.E.

	// Matrix of elastic properties
	D := "matrix(" +
		"E*t/(1-v*v)*1, E*t/(1-v*v)*v, 0, 0," +
		"E*t/(1-v*v)*v, E*t/(1-v*v)*1, 0, 0," +
		"0, 0, E*t*t*t/(12*pow(1-v,2))*1, E*t*t*t/(12*pow(1-v,2))*v," +
		"0, 0, E*t*t*t/(12*pow(1-v,2))*v, E*t*t*t/(12*pow(1-v,2))*1," +
		"4, 4)"
	fmt.Fprintf(os.Stdout, "D = %s\n", D)

	// Matrix of shape functions
	a := "(-e+3*e*e-2*e*e*e)"
	b := "(6-12*e)"
	c := "(6*e-6*e*e)"
	B := "matrix(" +
		"-cos(fi)/l, -sin(fi)/l, 0, cos(fi)/l, -sin(fi)/l, 0," +
		"sin(fi)*cos(fi)/r*" + a + "," +
		"1/r*(1-e-" + a + "*pow(cos(fi),2))," +
		"l*cos(fi)/r*(e-2*e*e+e*e*e)," +
		"sin(fi)*cos(fi)/r*(-" + a + ")," +
		"1/r*(e+" + a + "*pow(cos(fi),2))," +
		"l*cos(fi)/r*(-e*e+e*e*e)," +
		"sin(fi)/pow(l,2)*(+(" + b + "))," +
		"cos(fi)/pow(l,2)*(-(" + b + "))," +
		"1/l*(-4+6*e)," +
		"sin(fi)/pow(l,2)*(-(" + b + "))," +
		"cos(fi)/pow(l,2)*(+(" + b + "))," +
		"1/l*(-2+6*e)," +
		"sin(fi)*sin(fi)/(r*l)*(+(" + c + "))," +
		"sin(fi)*cos(fi)/(r*l)*(-(" + c + "))," +
		"sin(fi)/r*(1-4*e+3*e*e*e)," +
		"sin(fi)*sin(fi)/(r*l)*(-(" + c + "))," +
		"sin(fi)*cos(fi)/(r*l)*(+(" + c + "))," +
		"sin(fi)/r*(-2*e+3*e*e*e)," +
		"4, 6)"
	fmt.Fprintf(os.Stdout, "B = %s\n", B)

	K := cal("K", "transpose("+B+")*("+D+")*("+B+")")
	_ = K
}
