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

	// bendBeam()
	// trussHighOrder()
	// beam()

	Care()
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
	val, err = sm.Sexpr(
		os.Stdout,
		// nil,
		str)
	if err != nil {
		panic(err)
	}
	if !hide {
		fmt.Fprintf(os.Stdout, "Value     : %s\n", val)
		if m, ok := sm.ParseMatrix(val); ok {
			fmt.Fprintf(os.Stdout, "%s\n", m)
		}
	}
	return val
}

func createMatrix(r, c int) *sm.Matrix {
	return sm.CreateMatrix(r, c)
}

func bendBeam() {
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
		mK, ok := sm.ParseMatrix(s.K)
		if !ok {
			panic("not valid matrix")
		}
		fmt.Fprintf(os.Stdout, "%s without free:\n%s\n", s.name, mK)

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
				// wbend  = cal("displ", "matrix(v1, O1, v2, O2, 4,1)")
				// MVbend = cal("MV load", Kbend+" * "+wbend)

				// create vector with free dof
				free := createMatrix(4, 1)
				free.Args[free.Position(f, 0)] = sm.CreateFloat(1.0)
				// fmt.Fprintf(os.Stdout,"free:\n%s\n", free)

				gfactor := cal("gfactor", "( "+K+" * "+sm.AstToStr(free.Ast())+"); variable(x);constant(l)", show)

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
				gfactor = cal("gfactor minimal", gfactor+"/(-1*("+minGfactor+")) ; variable(x);constant(l)", show)

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
					" * "+sm.AstToStr(G.Ast())+"; variable(x);constant(l)", show)

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

func trussHighOrder() {
	// ELASTIC NONLINEAR ANALYSIS OF PLANE TRUSS BRIDGES
	// Morteza A.M. Torkamani , and Jyh-Hung Shieh

	show := true // false

	var (
		Nu  = cal("Nu", "matrix(1-x/L,0,x/L,0,1,4)", show)
		Nv  = cal("Nv", "matrix(0,1-x/L,0,x/L,1,4)", show)
		d   = cal("d", "matrix(u1,v1,u2,v2,4,1)", show)
		ux  = cal("ux", Nu+" * "+d, show)
		vx  = cal("vx", Nv+" * "+d, show)
		f   = cal("f", "matrix(Fx1,Fy1,Fx2,Fy2,4,1)", show)
		dNu = cal("dNu", "d("+Nu+",x);variable(x);", show)
		dNv = cal("dNv", "d("+Nv+",x);variable(x);", show)
		Ko  = cal("Ko", "integral(E*A*transpose("+dNu+")*"+dNu+",x,0,L);variable(x)", show)
		Kp  = cal("Kp", "P*integral(transpose("+dNu+")*"+dNu+
			"+ transpose("+dNv+")*"+dNv+",x,0,L);variable(x)", show)

		// Мне не нравиться что сумма матриц К1 и К2 не симметричная
		K1 = cal("K1", "integral(E*A/2*"+
			"transpose("+dNu+")*"+
			"("+
			"transpose("+d+")*"+"transpose("+dNu+")*"+dNu+"+"+
			"transpose("+d+")*"+"transpose("+dNv+")*"+dNv+
			")"+
			",x,0,L);variable(x);", show)
		K2 = cal("K2", "integral(E*A*("+
			"transpose("+dNu+")*"+dNu+"+transpose("+dNv+")*"+dNv+")*"+
			d+"*"+
			dNu+
			",x,0,L);variable(x)", show)

		K2v2 = cal("K2v2", "integral(E*A*"+
			dNu+"*"+"("+
			"transpose("+dNu+")*"+dNu+"*"+d+"+"+
			"transpose("+dNv+")*"+dNv+"*"+d+")"+
			",x,0,L);variable(x)", show)

		// for book:
		//
		// Theory and Analysis of Nonlinear Framed Structures
		// Yeong-Bin Yang, Shyh-Rong Kuo
		//
		// have another values.
		// TODO: - need approve -
		// нравиться что симметричная
		K3 = cal("K3", "integral(E*A/2*"+
			"(transpose("+dNu+")*"+dNu+"+transpose("+dNv+")*"+dNv+")*"+
			d+"*"+
			"transpose("+d+")*"+
			"(transpose("+dNu+")*"+dNu+"+transpose("+dNv+")*"+dNv+")*"+
			"1.0"+
			",x,0,L);variable(x)", show)

		_ = cal("Ko simplification", Ko+"*L/(E*A)", show)
		_ = cal("Kp simplification", Kp+"*L/P", show)
		_ = cal("K1 simplification", K1+"*2*L*L/E*1/A", show)
		_ = cal("K2 simplification", K2+"*2*L*L/E*1/A", show)
		_ = cal("K2v2 simplification", K2v2+"*2*L*L/E*1/A", show)
		_ = cal("K3 simplification", K3+"*2*L*L*L/E*1/A", show)
	)
	_ = ux
	_ = vx
	_ = f
	_ = Ko
	_ = Kp
	_ = K1
	_ = K2
	_ = K3
}

func beam() {
	// Theory and Analysis of Nonlinear Framed Structures
	// Yeong-Bin Yang, Shyh-Rong Kuo

	a1 := cal("a1", "(pow(-u1+u2,2)+pow(-v1+v2,2))*du")
	_ = a1

	a2 := cal("a2", "1.0*d("+a1+",u1); variable(u1)")
	_ = a2

	//  a3 := cal("a3", "1.0*d("+a1+",v1); variable(v1)")
	//  _ = a3

	//  a4 := cal("a4", "1.0*d("+a1+",u2); variable(u2)")
	//  _ = a4

	//  a5 := cal("a5", "1.0*d("+a1+",v2); variable(v2)")
	//  _ = a5

	s1 := cal("s1", "du*u1+dv*v1-du*u2-dv*v2")
	_ = s1

	//  s2 := cal("s2", "-(u2-u1)*u1-(v2-v1)*v1+(u2-u1)*u2+(v2-v1)*v2")
	//  _ = s2

	// show := true// false
	//
	// var (
	// 	// матрица коэффициентов
	// 	L = cal("[L]", `
	// 	matrix(
	// 		1,      0,      0,     0,
	// 		0,      1,      0,     0,
	// 		1,      l,    l*l, l*l*l,
	// 		0,      1,    2*l, 3*l*l,
	// 		4, 4)
	// 	`, show)
	//
	// 	// обратная матрица коэффициентов
	// 	invL = cal("[L]^-1", "inverse("+L+")", show)
	//
	// 	// функция формы конечного элемента
	// 	Ψ = cal("Ψbend", "matrix(1,x,x*x,x*x*x,1,4)*"+invL, show)
	// )
	//
	// _ = Ψ
}

func Care() {
	show := true // false
	u := cal("u", "a1 + a2*x", show)
	w := cal("w", "a3 + a4*x + a5*x*x + a6*x*x*x", show)
	constants := ";constant(a1,a2,a3,a4,a5,a6,L);variable(x);"
	// du := cal("du", "d("+u+",x)"+constants, show)
	dw := cal("dw", "d("+w+",x)"+constants, show)
	// ddw := cal("ddw", "d("+dw+",x)"+constants, show)
	// e := cal("e", du+"+0.5*"+dw+constants, show)
	// ro := cal("ro", "-("+ddw+")/pow("+"("+dw+")*("+dw+")"+",3/2)", show)

	fmt.Println("The generalised displacements")
	q1 := cal("q1 =  u(0)", "inject("+u+",x,0)"+constants, show)
	q2 := cal("q2 =  w(0)", "inject("+w+",x,0)"+constants, show)
	q3 := cal("q3 = dw(0)", "inject("+dw+",x,0)"+constants, show)
	q4 := cal("q4 =  u(L)", "inject("+u+",x,L)"+constants, show)
	q5 := cal("q5 =  w(L)", "inject("+w+",x,L)"+constants, show)
	q6 := cal("q5 = dw(L)", "inject("+dw+",x,L)"+constants, show)
	q := []string{q1, q2, q3, q4, q5, q6}

	fmt.Println("The unknown coefficients")
	a := "matrix(a1,a2,a3,a4,a5,a6,6,1)"
	am, _ := sm.ParseMatrix(a)

	fmt.Println("The connectivity matrix")
	C := "matrix("
	for row := 0; row < 6; row++ {
		fmt.Println("========= row ", row, " : ", q[row])
		for col := 0; col < 6; col++ {
			acol := sm.AstToStr(am.Args[col])
			cij := cal(fmt.Sprintf("[%d,%d]", row, col),
				"integral(d("+q[row]+","+acol+"), "+acol+",0,"+acol+")/"+acol+";"+
					"variable("+acol+")",
				show)
			C += cij + ","
		}
	}
	C += "6,6)"
	fmt.Println("C = ", C)
	invC := cal("inverse C", "inverse("+C+")", show)

	ai := cal("invC*q", invC+"*matrix(q1,q2,q3,q4,q5,q6,6,1)", show)

	A1 := "matrix(1,x,0,0,0,0,6,1)"
	u = cal("u new", "transpose("+A1+")*"+ai, show)
	A2 := "matrix(0,0,1,x,x*x,x*x*x,6,1)"
	w = cal("w new", "transpose("+A2+")*"+ai, show)
	du := cal("du new", "d("+u+",x)"+constants, show)
	_ = du
	dw = cal("dw new", "d("+w+",x)"+constants, show)
	ddw := cal("ddw new", "d("+dw+",x)"+constants, show)
	_ = ddw

	{
		V := cal("Axial strain classic",
			"EA/2*integral(pow("+du+",2),x,0,L)"+constants+
				";constant(q1,q2,q3,q4,q5,q6);",
			show)
		K1 := cal("K1",
			"d("+V+",q1);"+
				"constant(EA,L,x); variable(q1)",
			show)
		K11 := cal("K11",
			"integral(d("+K1+",q1),q1,0,q1)/q1;"+
				"constant(EA,L,x,q4); variable(q1)",
			show)
		K14 := cal("K14",
			"integral(d("+K1+",q4),q4,0,q4)/q4;"+
				"constant(EA,L,x,q1); variable(q4)",
			show)
		_ = V
		_ = K1
		_ = K11
		_ = K14
	}
	{
		V := cal("Axial strain with added part",
			"EA/2*integral(pow("+du+"+0.5*pow("+dw+",2),2),x,0,L)"+constants+
				";constant(q1,q2,q3,q4,q5,q6);",
			show)
		// 		K1 := cal("K1",
		// 			"d("+V+",q1);"+
		// 				"constant(EA,L,x); variable(q1)",
		// 			show)
		// 		K11 := cal("K11",
		// 			"integral(d("+K1+",q1),q1,0,q1)/q1;"+
		// 				"constant(EA,L,x,q4); variable(q1)",
		// 			show)
		// 		K14 := cal("K14",
		// 			"integral(d("+K1+",q4),q4,0,q4)/q4;"+
		// 				"constant(EA,L,x,q1); variable(q4)",
		// 			show)
		_ = V
		// 		_ = K1
		// 		_ = K11
		// 		_ = K14
	}

	//	V := cal("V1", "EA/2*integral(pow(0.5*pow("+dw+",2),2),x,0,L)"+constants+
	//		";constant(q1,q2,q3,q4,q5,q6);",
	//		show)
	/// V := cal("V1", "EA/2*integral(pow("+du+"+0.5*pow("+dw+",2),2),x,0,L)"+constants, show)
	// V := cal("V1", "EA/2*integral(pow("+ddw+",2),x,0,L)"+constants, show)

	// cal("u = A1^T*a", "transpose("+A1+")*"+a, show)
	// cal("w = A2^T*a", "transpose("+A2+")*"+a, show)
	//
	// da1 := cal("dA1/dx*invC", "transpose(d("+A1+",x))*"+invC+constants, show)
	// cal("dA2/dx*invC = 0", "inject("+da1+",x,0)"+constants, show)
	// cal("dA2/dx*invC = L", "inject("+da1+",x,L)"+constants, show)
	// da2 := cal("dA2/dx*invC", "transpose(d("+A2+",x))*"+invC+constants, show)
	// cal("dA2/dx*invC = 0", "inject("+da2+",x,0)"+constants, show)
	// cal("dA2/dx*invC = L", "inject("+da2+",x,L)"+constants, show)
	// dda2 := cal("ddA2/dx*invC", "transpose(d(d("+A2+",x),x))*"+invC+constants, show)
	// cal("ddA2/dx*invC = 0", "inject("+dda2+",x,0)"+constants, show)
	// cal("ddA2/dx*invC = L", "inject("+dda2+",x,L)"+constants, show)

	// w2 := cal("w2", "transpose("+dda2+")*("+dda2+")", show)
	// V := cal("V short", "0.5*EA*integral("+
	// 	// "transpose("+da1+")*"+da1+"+"+
	// 	//"0.5*transpose("+da1+")*("+w2+")+"+
	// 	"0.25*transpose("+w2+")*("+w2+")+"+
	// 	"0"+
	// 	",x,0,L)"+constants, show)
	// _ = w2
	// // V := cal("V short", "EA/2*integral(pow("+da1+"+0.5*pow("+da2+",2),2),x,0,L)+"+
	// // "EJ/2*integral(pow("+dda2+",2),x,0,L)"+constants, show)
	// _ = V

	// {
	// 	U := cal("U", "EJ/2*integral(("+ddw+")*("+ddw+"),x,0,L)"+constants, show)
	// 	_ = U
	// }
}
