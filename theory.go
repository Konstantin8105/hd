// +build ignore

package main

import (
	"fmt"
	"os"

	"github.com/Konstantin8105/sm"
)

func main() {
	sm.MaxIteration = -1
	sm.FloatFormat = 5

	bendBeam()
}

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
	val, err = sm.Sexpr(nil, str)
	if err != nil {
		panic(err)
	}
	if !hide {
		fmt.Fprintf(os.Stdout, "Value     : %s\n", val)
	}
	return val
}

func createMatrix(r, c int) *sm.Matrix {
	return sm.CreateMatrix(r, c)
}

func bendBeam() {
	// w    = cal("displ", "matrix(u1, v1, O1, u2, v2, O2, 6,1)")
	// load = cal("load", "matrix( N1, V1, M1, N2, V2, M2, 6,1)")
	// a    = cal("a", "matrix(a1,a2,a3,a4,a5,a6,6,1)")
	// v    = cal("перемещения v", "transpose("+a+") * matrix(1,x,x*x,x*x*x, 0,0 , 6,1)")
	// u    = cal("перемещения u", "transpose("+a+") * matrix(0,0,0,0,1,x, 6,1)")

	var (
		// матрица коэффициентов
		L = cal("[L]", `
		matrix(
			1,      0,      0,     0,
			0,     -1,      0,     0,
			1,      l,    l*l, l*l*l,
			0,     -1,   -2*l,-3*l*l,
			4, 4)
		`)

		// обратная матрица коэффициентов
		invL = cal("[L]^-1", "inverse("+L+")")

		// функция формы конечного элемента
		Ψ      = cal("Ψbend", "matrix(1,x,x*x,x*x*x,1,4)*"+invL)
		dFdx   = cal("dFdx", "d("+Ψ+", x); variable(x)")
		d2Fdx2 = cal("d2Fdx2", "d("+dFdx+", x); variable(x); constant(l);")

		// матрица жесткости
		Ko = cal("Ko", "EJ*integral( transpose("+d2Fdx2+") * "+d2Fdx2+",x,0,l);variable(x); constant(l);")

		// геометрическая матрица жесткости
		Kg = cal("Kg", "N*integral( transpose("+dFdx+") * "+dFdx+",x,0,l);variable(x); constant(l);")
	)

	show := false

	for matrixIndex, s := range []struct {
		name string
		K    string
	}{
		{"Matrix stiffnes", Ko},
		{"Geometric matrix stiffness", Kg},
	} {
		fmt.Fprintln(os.Stdout, "|***************************************************|")
		mK, ok := sm.ParseMatrix(s.K)
		if !ok {
			panic("not valid matrix")
		}
		fmt.Fprintf(os.Stdout, "Matrix K without free:\n%s\n", mK)

		// степень свободы
		for _, freeIndex := range [][]int{
			{0},
			{1},
			{2},
			{3},
			{3, 1},
		} {
			fmt.Fprintln(os.Stdout, "|---------------------------------------------|")
			fmt.Fprintln(os.Stdout, "Name :", s.name, " with free ", freeIndex)
			K := s.K

			for _, f := range freeIndex {
				// wbend  = cal("displ", "matrix(v1, O1, v2, O2, 4,1)")
				// MVbend = cal("MV load", Kbend+" * "+wbend)

				// create vector with free dof
				free := createMatrix(4, 1)
				free.Args[free.Position(f, 0)] = sm.CreateFloat(1.0)
				// fmt.Fprintf(os.Stdout,"free:\n%s\n", free)

				gfactor := cal("gfactor", "( "+K+" * "+sm.AstToStr(free.Ast())+")", show)

				mG, ok := sm.ParseMatrix(gfactor)
				if !ok {
					panic("not valid matrix")
				}
				if show {
					fmt.Fprintf(os.Stdout, "Matrix gfactor:\n%s\n", mG)
				}
				minGfactor := sm.AstToStr(mG.Args[mG.Position(f, 0)])
				if show {
					fmt.Fprintln(os.Stdout, "minGfactor = ", minGfactor)
				}
				gfactor = cal("gfactor minimal", gfactor+"/(-1*("+minGfactor+"))", show)

				mG, ok = sm.ParseMatrix(gfactor)
				if !ok {
					panic("not valid matrix")
				}
				if show {
					fmt.Fprintf(os.Stdout, "Matrix gfactor:\n%s\n", mG)
				}

				G := sm.CreateMatrix(4, 4)
				for i := 0; i < 4; i++ {
					G.Args[G.Position(i, i)] = sm.CreateFloat(1.0)
				}
				for i := 0; i < 4; i++ {
					G.Args[G.Position(f, i)] = mG.Args[mG.Position(i, 0)]
				}
				G.Args[G.Position(f, f)] = sm.CreateFloat(0.0)
				if show {
					fmt.Fprintf(os.Stdout, "Matrix G:\n%s\n", G)
				}

				KwithFree := cal("KwithFree", "transpose("+sm.AstToStr(G.Ast())+") * "+K+
					" * "+sm.AstToStr(G.Ast()), show)

				K = KwithFree

				mKF, ok := sm.ParseMatrix(KwithFree)
				if !ok {
					panic("not valid matrix")
				}
				if matrixIndex == 0 {
					mKF.Args[mKF.Position(f, f)] = sm.CreateFloat(1.0)
				}
				fmt.Fprintf(os.Stdout, "K with free\n%s\n", mKF)

			}
		}
		fmt.Fprintln(os.Stdout, "||||||==========================================||||||")
	}
}
