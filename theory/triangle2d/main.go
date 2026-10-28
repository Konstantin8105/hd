package main

import (
	"fmt"
	"os"
	"strings"

	"github.com/Konstantin8105/sm"
)

// run :
//	go run ./theory/triangle2d

func cal(name, str string, show bool) string {
	if show {
		fmt.Fprintf(os.Stdout, "\nName    : %s\n", name)
		fmt.Fprintf(os.Stdout, "\nFormula : %s\n", str)
	}
	val, err := sm.Sexpr(
		// os.Stdout,
		nil,
		str)
	if err != nil {
		panic(err)
	}
	if show {
		fmt.Fprintf(os.Stdout, "Value   : %s\n", val)
	}
	return val
}

// vector returns the `sm` expression of a 6x1 unit vector with one
// at the position `f`.
func vector(f int) string {
	var b strings.Builder
	b.WriteString("matrix(")
	for i := 0; i < 6; i++ {
		if 0 < i {
			b.WriteString(", ")
		}
		if i == f {
			b.WriteString("1")
		} else {
			b.WriteString("0")
		}
	}
	b.WriteString(", 6, 1)")
	return b.String()
}

// value returns the single element of a 1x1 matrix expression.
// `sm` prints such a matrix as: matrix(<element>, 1.00000, 1.00000)
func value(m string) string {
	m = strings.TrimSpace(m)
	const prefix = "matrix("
	if !strings.HasPrefix(m, prefix) || !strings.HasSuffix(m, ")") {
		panic("not a matrix: " + m)
	}
	body := m[len(prefix) : len(m)-1]
	i := strings.LastIndex(body, ",")
	j := strings.LastIndex(body[:i], ",")
	return strings.TrimSpace(body[:j])
}

func main() {

	sm.MaxIteration = -1
	sm.FloatFormat = 4

	show := false // true

	// -----------------------------------------------------------------
	// Constant Strain Triangle (CST) finite element in a 2D plane.
	//
	// Nodes : 1 = (x1, y1), 2 = (x2, y2), 3 = (x3, y3)
	// DoF   : (u1, v1, u2, v2, u3, v3)
	//
	// The displacement field is linear inside the triangle
	//
	//	u(x,y) = a + b*x + c*y
	//	v(x,y) = d + e*x + f*y
	//
	// so the strains (and the B matrix) are constant over the element.
	//
	// Geometry symbols used below:
	//
	//	b1 = y2-y3, b2 = y3-y1, b3 = y1-y2
	//	c1 = x3-x2, c2 = x1-x3, c3 = x2-x1
	//	d  = 2A = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)  (twice the area)
	// -----------------------------------------------------------------

	// strain-displacement matrix without the 1/(2A) factor:
	//	B = (1/d) * B0
	B0 := "matrix(" +
		"b1, 0, b2, 0, b3, 0," +
		"0, c1, 0, c2, 0, c3," +
		"c1, b1, c2, b2, c3, b3," +
		"3, 6)"

	// geometric definitions (informative)
	_ = cal("d = 2A", "(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)", show)
	_ = cal("b1", "y2-y3", show)
	_ = cal("c1", "x3-x2", show)

	// dimensionless matrix of elastic properties (plane stress)
	Dnum := "matrix(1,v,0, v,1,0, 0,0,(1-v)/2, 3,3)"

	// matrix of elastic properties (plane stress)
	_ = cal("D (plane stress)", "E/(1-v*v)*("+Dnum+")", show)

	// matrix of elastic properties (plane strain)
	_ = cal("D (plane strain)",
		"E/((1+v)*(1-2*v))*matrix(1-v,v,0, v,1-v,0, 0,0,(1-2*v)/2, 3,3)", show)

	// strain-displacement matrix
	_ = cal("B", "("+B0+")/d", show)

	// The stiffness matrix is
	//
	//	K = t*A*Bᵀ*D*B = s * Knum
	//
	// where
	//
	//	s    = E*t / (2*d*(1-v*v))
	//	Knum = B0ᵀ*Dnum*B0
	//
	// The common scalar `s` is factored out, so that the (symbolically
	// heavy) pin condensation runs on the simpler matrix `Knum`.
	s := "E*t/(2*d*(1-v*v))"

	Knum := cal("Knum", "transpose("+B0+")*("+Dnum+")*("+B0+")", show)

	K := cal("K", "("+s+")*("+Knum+")", show)
	fmt.Fprintf(os.Stdout, "K = %s\n", K)

	// -----------------------------------------------------------------
	// Pins (release of a DoF) by static condensation.
	//
	// For a released DoF `f`:
	//
	//	Knew = K - K[:,f] * K[f,:] / K[f,f]
	//
	// which zeroes the row and the column `f`. The same transform is
	// applied to `Knum` (the scalar `s` cancels), then `s` is restored.
	// -----------------------------------------------------------------
	for _, pins := range [][]int{
		{0},
		{1},
		{2},
		{3},
		{4},
		{5},
	} {
		fmt.Fprintln(os.Stdout, "|---------------------------------------------|")
		fmt.Fprintln(os.Stdout, "Pins:", pins)

		KN := Knum
		for _, f := range pins {
			free := vector(f)
			fmt.Fprintf(os.Stdout, "free:\n%s\n", free)

			col := cal("column", "("+KN+") * "+free, show)
			row := cal("row", "transpose("+free+") * ("+KN+")", show)
			kff := value(cal("kff", "transpose("+free+") * ("+KN+") * "+free, show))

			KN = cal("Knum with pin", "("+KN+") - ("+col+") * ("+row+") / ("+kff+")", show)

			Kpin := cal("K with pin", "("+s+")*("+KN+")", show)

			fmt.Fprintf(os.Stdout, "K with pin %d\n%s\n", f, Kpin)
		}
	}
}
