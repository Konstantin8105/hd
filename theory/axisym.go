// +build ignore

package main

import (
	"fmt"
	"os"

	"github.com/Konstantin8105/sm"
)

// run : 
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
	constants := ";constant(a1,a2,a3,a4,a5,a6,l);"
	u := "a1+a2*s"
	w := "a3+a4*s+a5*s*s+a6*s*s*s"
	b := cal("b","d("+w+",s)"+constants+";variable(s);")

	fmt.Println("---------------------------------")
	fmt.Println("s = 0")
	cal("u","inject("+u+",s,0)")
	cal("w","inject("+w+",s,0)")
	cal("b","inject("+b+",s,0)")
	fmt.Println("---------------------------------")
	fmt.Println("s = l")
	cal("u","inject("+u+",s,l)")
	cal("w","inject("+w+",s,l)")
	cal("b","inject("+b+",s,l)")
	fmt.Println("---------------------------------")
}
