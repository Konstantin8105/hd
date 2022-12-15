// +build ignore

package main

import (
	"fmt"
	"os"
	"strings"

	"github.com/Konstantin8105/sm"
)

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

func solve(code string) {
	lines := strings.Split(code, "\n")
	vars := map[string]string{}
	show := true
	for _, line := range lines {
		line = strings.TrimSpace(line)
		if line == "" {
			continue
		}
		epos := strings.Index(line, "=")
		name := strings.TrimSpace(line[:epos])
		eq := strings.TrimSpace(line[epos+1:])
		for k, v := range vars {
			eq = strings.ReplaceAll(eq, k, "("+v+")")
		}
		v := cal(name, eq, show)
		vars[name] = v
	}
}

func main() {
	sm.MaxIteration = -1
	sm.FloatFormat = 5

	// Research
	// Higher-order striffness matrices in nonlinear finite element analysis of plate truss Structures

	solve(`
vNu   = matrix(1-x/L,0,x/L,0,1,4);
vdNu  = d(vNu,x);variable(x);
vNv   = matrix(0,1-x/L,0,x/L,1,4);
vdNv  = d(vNv,x);variable(x);
Dv1    = Dv2-Δv
Du1    = Du2-Δu
vD    = matrix(Du1,Dv1,Du2,Dv2,4,1);
vKo   = integral(E*A*transpose(vdNu)*vdNu,x,0,L);variable(x);constant(L, E, A)
vKp   = P*integral(transpose(vdNu)*vdNu+transpose(vdNv)*vdNv,x,0,L);variable(x);constant(L, E, A, P)
shortKp = vKp*(L/P)
vK1   = integral(E*A/2*transpose(vdNu)*( transpose(vD)*transpose(vdNu)*vdNu  +   transpose(vD)*transpose(vdNv)*vdNv    ),x,0,L);variable(x); constant(L,E,A,P)
shortK1 = vK1*(2*L*L/(E*A))
vK2   = integral(E*A*( transpose(vdNu)*vdNu+transpose(vdNv)*vdNv   )*vD*vdNu ,x ,0,L   );variable(x); constant(L,E,A,P)
shortK2 = vK2*(2*L*L/(E*A))
vK3   = integral(E*A/2*(transpose(vdNu)*vdNu+transpose(vdNv)*vdNv) * vD*transpose(vD) * (transpose(vdNu)*vdNu+transpose(vdNv)*vdNv),x,0,L);variable(x); constant(L,E,A,P)
shortK3 = vK3*(2*L*L*L/(E*A))
`)
}
