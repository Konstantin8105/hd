// +build ignore

package main

import (
	"fmt"
	"go/ast"
	"os"

	"github.com/Konstantin8105/sm"
)

// run :
//	go run axisym.go
//	go run -mod=mod axisym.go

func cal(name, str string, view ...bool) string {
	hide := false
	if 0 < len(view) {
		hide = !view[0]
	}
	if !hide {
		fmt.Fprintf(os.Stdout, "\nName    : %s\n", name)
		fmt.Fprintf(os.Stdout, "\nFormula : %s\n", str)
	}
	var val string
	var err error
	val, err = sm.Sexpr(
		// os.Stdout,
		nil,
		str)
	if err != nil {
		panic(err)
	}
	if !hide {
		fmt.Fprintf(os.Stdout, "Value   : %s\n", val)
		if m, ok := sm.ParseMatrix(val); ok {
			fmt.Fprintf(os.Stdout, "%s\n", m)
		}
	}
	return val
}

func main() {

	sm.MaxIteration = -1
	sm.FloatFormat = 5


	constants := ";constant(a1,a2,a3,a4,a5,a6,l);"
	u := "a1+a2*s"
	w := "a3+a4*s+a5*s*s+a6*s*s*s"
	b := cal("b", "d("+w+",s)"+constants+";variable(s);")

	fmt.Println("---------------------------------")
	fmt.Println("s = 0")
	cal("u", "inject("+u+",s,0)")
	cal("w", "inject("+w+",s,0)")
	cal("b", "inject("+b+",s,0)")
	fmt.Println("---------------------------------")
	fmt.Println("s = l")
	cal("u", "inject("+u+",s,l)")
	cal("w", "inject("+w+",s,l)")
	cal("b", "inject("+b+",s,l)")
	fmt.Println("---------------------------------")

	// -----------------------
	// CONICAL FINITE ELEMENT STRESS ANALYSIS OF AN AXISYMMETRIC
	// TANK ON AN ELASTIC FOUNDATION
	// by
	// HOWARD RICHARD HORN, JR., B.S. in C.E., M.S. In C.E.

	D := sm.CreateMatrix(4, 4)
	{
		D.Args[D.Position(0, 0)] = ast.NewIdent("E*t/(1-v*v)*1")
		D.Args[D.Position(1, 1)] = ast.NewIdent("E*t/(1-v*v)*1")
		D.Args[D.Position(1, 0)] = ast.NewIdent("E*t/(1-v*v)*v")
		D.Args[D.Position(0, 1)] = ast.NewIdent("E*t/(1-v*v)*v")

		D.Args[D.Position(2, 2)] = ast.NewIdent("E*t*t*t/(12*pow(1-v,2))*1")
		D.Args[D.Position(3, 3)] = ast.NewIdent("E*t*t*t/(12*pow(1-v,2))*1")
		D.Args[D.Position(3, 2)] = ast.NewIdent("E*t*t*t/(12*pow(1-v,2))*v")
		D.Args[D.Position(2, 3)] = ast.NewIdent("E*t*t*t/(12*pow(1-v,2))*v")
	}
	fmt.Fprintf(os.Stdout, "D = %s\n", D)

	B := sm.CreateMatrix(4, 6)
	{
		a := "(-e+3*e*e-2*e*e*e)"
		b := "(6-12*e)"
		c := "(6*e-6*e*e)"
		B.Args[B.Position(0,0)] = ast.NewIdent("-cos(fi)/l")
		B.Args[B.Position(0,1)] = ast.NewIdent("-sin(fi)/l")
		B.Args[B.Position(0,3)] = ast.NewIdent("cos(fi)/l")
		B.Args[B.Position(0,4)] = ast.NewIdent("-sin(fi)/l")

		B.Args[B.Position(1,0)] = ast.NewIdent("sin(fi)*cos(fi)/r*"+a)
		B.Args[B.Position(1,1)] = ast.NewIdent("1/r*(1-e-"+a+"*pow(cos(fi),2))")
		B.Args[B.Position(1,2)] = ast.NewIdent("l*cos(fi)/r*(e-2*e*e+e*e*e)")
		B.Args[B.Position(1,3)] = ast.NewIdent("sin(fi)*cos(fi)/r*(-"+a+")")
		B.Args[B.Position(1,4)] = ast.NewIdent("1/r*(e+"+a+"*pow(cos(fi),2))")
		B.Args[B.Position(1,5)] = ast.NewIdent("l*cos(fi)/r*(-e*e+e*e*e)")

		B.Args[B.Position(2,0)] = ast.NewIdent("sin(fi)/pow(l,2)*(+("+b+"))")
		B.Args[B.Position(2,1)] = ast.NewIdent("cos(fi)/pow(l,2)*(-("+b+"))")
		B.Args[B.Position(2,2)] = ast.NewIdent("1/l*(-4+6*e)")
		B.Args[B.Position(2,3)] = ast.NewIdent("sin(fi)/pow(l,2)*(-("+b+"))")
		B.Args[B.Position(2,4)] = ast.NewIdent("cos(fi)/pow(l,2)*(+("+b+"))")
		B.Args[B.Position(2,5)] = ast.NewIdent("1/l*(-2+6*e)")

		B.Args[B.Position(3,0)] = ast.NewIdent("sin(fi)*sin(fi)/(r*l)*(+("+c+"))")
		B.Args[B.Position(3,1)] = ast.NewIdent("sin(fi)*cos(fi)/(r*l)*(-("+c+"))")
		B.Args[B.Position(3,2)] = ast.NewIdent("sin(fi)/r*(1-4*e+3*e*e*e)")
		B.Args[B.Position(3,3)] = ast.NewIdent("sin(fi)*sin(fi)/(r*l)*(-("+c+"))")
		B.Args[B.Position(3,4)] = ast.NewIdent("sin(fi)*cos(fi)/(r*l)*(+("+c+"))")
		B.Args[B.Position(3,5)] = ast.NewIdent("sin(fi)/r*(-2*e+3*e*e*e)")
	}
	fmt.Fprintf(os.Stdout, "B = %s\n", B.String())

	K := cal("K", "transpose("+sm.AstToStr(B.Ast())+")*("+sm.AstToStr(D.Ast())+")*("+sm.AstToStr(B.Ast())+")")
	_ = K
}
