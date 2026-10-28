package main

import (
	"fmt"
	"os"
	"strings"

	"github.com/Konstantin8105/sm"
)

// run :
//	go run ./theory/beam2d

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

// vector returns the `sm` expression of a 4x1 unit vector with one
// at the position `f`.
func vector(f int) string {
	var b strings.Builder
	b.WriteString("matrix(")
	for i := 0; i < 4; i++ {
		if 0 < i {
			b.WriteString(", ")
		}
		if i == f {
			b.WriteString("1")
		} else {
			b.WriteString("0")
		}
	}
	b.WriteString(", 4, 1)")
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

	show := true // false

	var (
		// матрица коэффициентов
		L = cal("[L]", `
		matrix(
			1,      0,      0,     0,
			0,      1,      0,     0,
			1,      l,    l*l, l*l*l,
			0,      1,    2*l, 3*l*l,
			4, 4)
		`, show)

		// обратная матрица коэффициентов
		invL = cal("[L]^-1", "inverse("+L+")", show)

		// функция формы конечного элемента
		Ψ      = cal("Ψbend", "matrix(1,x,x*x,x*x*x,1,4)*"+invL, show)
		dFdx   = cal("dFdx", "d("+Ψ+", x); variable(x)", show)
		d2Fdx2 = cal("d2Fdx2", "d("+dFdx+", x); variable(x); constant(l);", show)

		// матрица жесткости
		Ko = cal("Ko", "EJ*integral( transpose("+d2Fdx2+") * "+d2Fdx2+",x,0,l);variable(x); constant(l);", show)

		// геометрическая матрица жесткости
		Kg = cal("Kg", "N*integral( transpose("+dFdx+") * "+dFdx+",x,0,l);variable(x); constant(l);", show)
	)

	for matrixIndex, s := range []struct {
		name string
		K    string
	}{
		{"Matrix stiffnes", Ko},
		{"Geometric matrix stiffness", Kg},
	} {
		fmt.Fprintln(os.Stdout, "|***************************************************|")
		fmt.Fprintf(os.Stdout, "%s without free:\n%s\n", s.name, s.K)

		// степень свободы
		for _, freeIndex := range [][]int{
			{0},
			{1},
			{2},
			{3},
			{3, 1},
			{0, 1},
			{2, 3},
		} {
			fmt.Fprintln(os.Stdout, "|---------------------------------------------|")
			fmt.Fprintln(os.Stdout, "Name :", s.name, " with free ", freeIndex)
			K := s.K

			for _, f := range freeIndex {
				// create vector with free dof
				free := vector(f)
				fmt.Fprintf(os.Stdout, "free:\n%s\n", free)

				// column and row of the freedom `f`
				col := cal("column", "("+K+") * "+free, show)
				row := cal("row", "transpose("+free+") * ("+K+")", show)
				kff := value(cal("kff", "transpose("+free+") * ("+K+") * "+free, show))

				// release of the freedom `f`:
				//   Knew = K - K[:,f] * K[f,:] / K[f,f]
				K = cal("KwithFree", "("+K+") - ("+col+") * ("+row+") / ("+kff+")", show)

				// for the stiffness matrix the released diagonal term
				// is printed as 1
				out := K
				if matrixIndex == 0 {
					out = cal("KwithFree diagonal", "("+K+") + ("+free+") * transpose("+free+")", show)
				}
				fmt.Fprintf(os.Stdout, "K with free\n%s\n", out)
			}
		}
		fmt.Fprintln(os.Stdout, "||||||==========================================||||||")
	}
}
