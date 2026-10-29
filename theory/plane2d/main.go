package main

import (
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/Konstantin8105/sm"
)

// run :
//	go run ./theory/plane2d

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

// vector returns the `sm` expression of a 8x1 unit vector with one
// at the position `f`.
func vector(f int) string {
	var b strings.Builder
	b.WriteString("matrix(")
	for i := 0; i < 8; i++ {
		if 0 < i {
			b.WriteString(", ")
		}
		if i == f {
			b.WriteString("1")
		} else {
			b.WriteString("0")
		}
	}
	b.WriteString(", 8, 1)")
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

// splitTopLevel splits by commas that are not inside brackets.
func splitTopLevel(s string) []string {
	var (
		parts []string
		depth int
		start int
	)
	for i, r := range s {
		switch r {
		case '(', '[', '{':
			depth++
		case ')', ']', '}':
			depth--
		case ',':
			if depth == 0 {
				parts = append(parts, s[start:i])
				start = i + 1
			}
		}
	}
	return append(parts, s[start:])
}

// parseMatrix parses the `sm` matrix expression:
//
//	matrix(<elem>, <elem>, ..., <rows>, <cols>)
func parseMatrix(str string) (rows, cols int, elems []string, ok bool) {
	str = strings.TrimSpace(str)
	const prefix = "matrix("
	if !strings.HasPrefix(str, prefix) || !strings.HasSuffix(str, ")") {
		return 0, 0, nil, false
	}
	parts := splitTopLevel(str[len(prefix) : len(str)-1])
	if len(parts) < 3 {
		return 0, 0, nil, false
	}
	r, err := strconv.ParseFloat(strings.TrimSpace(parts[len(parts)-2]), 64)
	if err != nil {
		return 0, 0, nil, false
	}
	c, err := strconv.ParseFloat(strings.TrimSpace(parts[len(parts)-1]), 64)
	if err != nil {
		return 0, 0, nil, false
	}
	rows, cols = int(r), int(c)
	if rows <= 0 || cols <= 0 || len(parts)-2 != rows*cols {
		return 0, 0, nil, false
	}
	for i := range parts[:len(parts)-2] {
		parts[i] = strings.TrimSpace(parts[i])
	}
	return rows, cols, parts[:len(parts)-2], true
}

func main() {

	sm.MaxIteration = -1
	sm.FloatFormat = 12

	show := false // true

	// -----------------------------------------------------------------
	// 4-node bilinear quadrilateral (Q4) finite element in a 2D plane.
	//
	// Nodes : 1, 2, 3, 4 (counter clockwise)
	// DoF   : (u1, v1, u2, v2, u3, v3, u4, v4)
	//
	// Natural coordinates: (xi, eta) in [-1, 1] x [-1, 1]
	// Bilinear shape functions:
	//
	//	N1 = (1-xi)*(1-eta)/4     N2 = (1+xi)*(1-eta)/4
	//	N3 = (1+xi)*(1+eta)/4     N4 = (1-xi)*(1+eta)/4
	//
	// The element is a parallelogram (the isoparametric map is affine):
	//
	//	node 1 = (0, 0)
	//	node 2 = (u1, u2)
	//	node 3 = (u1+v1, u2+v2)
	//	node 4 = (v1, v2)
	//
	// so the Jacobian J is constant and the strain field is linear.
	// -----------------------------------------------------------------

	// bilinear shape functions
	N1 := "(1-xi)*(1-eta)/4"
	N2 := "(1+xi)*(1-eta)/4"
	N3 := "(1+xi)*(1+eta)/4"
	N4 := "(1-xi)*(1+eta)/4"
	fmt.Fprintf(os.Stdout, "N = [%s, %s, %s, %s]\n", N1, N2, N3, N4)

	// isoparametric map
	x := N2 + "*u1+" + N3 + "*(u1+v1)+" + N4 + "*v1"
	y := N2 + "*u2+" + N3 + "*(u2+v2)+" + N4 + "*v2"

	cd := ";constant(eta);constant(u1);constant(u2);constant(v1);constant(v2)"
	ce := ";constant(xi);constant(u1);constant(u2);constant(v1);constant(v2)"

	dxdxi := cal("dx/dxi", "d("+x+",xi);variable(xi)"+cd, show)
	dydxi := cal("dy/dxi", "d("+y+",xi);variable(xi)"+cd, show)
	dxdeta := cal("dx/deta", "d("+x+",eta);variable(eta)"+ce, show)
	dydeta := cal("dy/deta", "d("+y+",eta);variable(eta)"+ce, show)

	// Jacobian of the map and its determinant
	J := "matrix(" + dxdxi + ", " + dydxi + ", " + dxdeta + ", " + dydeta + ", 2, 2)"
	fmt.Fprintf(os.Stdout, "J = %s\n", J)
	detJ := cal("detJ", "det("+J+")", show)
	fmt.Fprintf(os.Stdout, "detJ = %s\n", detJ)

	// derivatives of the shape functions in (xi, eta)
	p := [4]string{"-(1-eta)/4", "(1-eta)/4", "(1+eta)/4", "-(1+eta)/4"} // dNi/dxi
	q := [4]string{"-(1-xi)/4", "-(1+xi)/4", "(1+xi)/4", "(1-xi)/4"}     // dNi/deta

	// numerators of the shape function derivatives in (x, y):
	//
	//	[dNi/dx]         [ dydeta  -dydxi ] [dNi/dxi ]
	//	[dNi/dy] = 1/J * [ -dxdeta  dxdxi ] [dNi/deta]
	//
	// here only the numerators are kept (without 1/J).
	gnumx := [4]string{}
	gnumy := [4]string{}
	for i := 0; i < 4; i++ {
		gnumx[i] = cal(fmt.Sprintf("dN%d/dx", i+1),
			"("+dydeta+")*("+p[i]+")-("+dydxi+")*("+q[i]+")", show)
		gnumy[i] = cal(fmt.Sprintf("dN%d/dy", i+1),
			"-("+dxdeta+")*("+p[i]+")+("+dxdxi+")*("+q[i]+")", show)
	}

	// strain-displacement matrix without the 1/J factor:
	//	B = (1/detJ) * Bnum
	Bnum := "matrix(" +
		gnumx[0] + ", 0, " + gnumx[1] + ", 0, " + gnumx[2] + ", 0, " + gnumx[3] + ", 0," +
		"0, " + gnumy[0] + ", 0, " + gnumy[1] + ", 0, " + gnumy[2] + ", 0, " + gnumy[3] + "," +
		gnumy[0] + ", " + gnumx[0] + ", " + gnumy[1] + ", " + gnumx[1] + ", " +
		gnumy[2] + ", " + gnumx[2] + ", " + gnumy[3] + ", " + gnumx[3] + "," +
		"3, 8)"
	fmt.Fprintf(os.Stdout, "Bnum = %s\n", Bnum)

	// dimensionless matrix of elastic properties (plane stress)
	Dnum := "matrix(1,v,0, v,1,0, 0,0,(1-v)/2, 3,3)"

	// matrix of elastic properties (plane stress)
	_ = cal("D (plane stress)", "E/(1-v*v)*("+Dnum+")", show)

	// matrix of elastic properties (plane strain)
	_ = cal("D (plane strain)",
		"E/((1+v)*(1-2*v))*matrix(1-v,v,0, v,1-v,0, 0,0,(1-2*v)/2, 3,3)", show)

	// The stiffness matrix is
	//
	//	K = t * ∫∫ Bᵀ*D*B * detJ dxi deta = s * M
	//
	// where
	//
	//	s = E*t / (1-v*v) / detJ
	//	M = ∫∫ Bnumᵀ * Dnum * Bnum dxi deta
	//
	// The double integral is the exact 2x2 Gauss integration, which is
	// exact because the integrand is quadratic in xi and eta.
	Q := cal("Q", "transpose("+Bnum+")*("+Dnum+")*("+Bnum+")", show)
	M2 := cal("M2",
		"integral(integral("+Q+", xi, -1, 1), eta, -1, 1);variable(xi);variable(eta)", show)

	s := "E*t/((1-v*v)*(" + detJ + "))"

	K := cal("K", "("+s+")*("+M2+")", show)
	fmt.Fprintf(os.Stdout, "K = %s\n", K)

	// -----------------------------------------------------------------
	// Pins (release of a DoF) by static condensation.
	//
	// For a released DoF `f`:
	//
	//	Knew = K - K[:,f] * K[f,:] / K[f,f]
	//
	// which zeroes the row and the column `f`. The common scalar `s`
	// cancels, so the condensation is done on `M2`.
	//
	// The condensation is symbolically heavy for the full `M2`, so it
	// is done once on a small matrix of symbols m_i_j, which are then
	// replaced back by the elements of `M2`.
	// -----------------------------------------------------------------
	rows, cols, elems, ok := parseMatrix(M2)
	if !ok {
		panic("M2 is not a matrix")
	}
	symbol := func(i, j int) string { return fmt.Sprintf("m_%d_%d", i, j) }
	{
		var list []string
		for i := 0; i < rows; i++ {
			for j := 0; j < cols; j++ {
				list = append(list, symbol(i, j))
			}
		}
		smMatrix := "matrix(" + strings.Join(list, ", ") +
			fmt.Sprintf(", %d, %d)", rows, cols)

		for _, pins := range [][]int{
			{0},
			{1},
			{2},
			{3},
			{4},
			{5},
			{6},
			{7},
		} {
			fmt.Fprintln(os.Stdout, "|---------------------------------------------|")
			fmt.Fprintln(os.Stdout, "Pins:", pins)

			MS := smMatrix
			for _, f := range pins {
				free := vector(f)
				fmt.Fprintf(os.Stdout, "free:\n%s\n", free)

				col := cal("column", "("+MS+") * "+free, show)
				row := cal("row", "transpose("+free+") * ("+MS+")", show)
				kff := value(cal("kff", "transpose("+free+") * ("+MS+") * "+free, show))

				MS = cal("M2 sym with pin", "("+MS+") - ("+col+") * ("+row+") / ("+kff+")", show)
			}

			// replace symbols m_i_j by the elements of M2
			out := MS
			for i := 0; i < rows; i++ {
				for j := 0; j < cols; j++ {
					out = strings.ReplaceAll(out, symbol(i, j), elems[i*cols+j])
				}
			}

			fmt.Fprintf(os.Stdout, "K with pin %v\n(%s) * (%s)\n", pins, s, out)
		}
	}
}
