package main

import (
	"fmt"
	"math"
)

func main() {
	// octave main_arclength.m

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// clear all;
	// clc;
	// more off;
	// format long;
	//
	// %fname = "input_Truss_2D_3members_model1.txt";
	// %fname = "input_Truss_3D_2members.txt";
	// fname = "input_Truss_3D_12members.txt";
	// %fname = "input_LeeFrame-nelem20.txt";
	// %fname = "input_arch-215deg.txt";
	// %fname = "input_Arch_semicircle-nelem50-sym.txt";
	// %fname = "input_Arch_semicircle-nelem50-unsym.txt";
	// %fname = "input-beamEndMoment-nelem10.txt";
	//
	// [ndim, ndof, nnode, nelem, coords, elemConn, elemData, LM, neq, assy4r, \
	// dof_force, Fext, maxloadSteps, loadincr, outputlist] = processfile(fname)
	ndim, ndof, nnode, nelem, coords, elemConn, elemData, LM, neq, assy4r,
		dof_force, Fext, maxloadSteps, loadincr, outputlist := Processfile()
	//
	// disp = zeros(neq,1);
	disp := make([]float64, neq)
	//
	dispPrev := make([]float64, neq)  // disp;
	dispPrev2 := make([]float64, neq) // disp;
	dispPrev3 := make([]float64, neq) // disp;
	dispPrev4 := make([]float64, neq) // disp;
	// 	//
	// 	var Kglobal [][]float64 = make([][]float64, neq)
	// 	for i := 0; i < neq; i++ {
	// 		Kglobal[i] = make([]float64, neq)
	// 	}
	// 	// Kglobal = zeros(neq,neq);
	// 	Rglobal := make([]float64, neq) // = zeros(neq,1);
	//
	// bf=[0.0 0.0];
	var bf [2]float64
	//
	var (
		Ds     = loadincr
		DsPrev = Ds
		DsMax  = Ds
		DsMin  = Ds

		loadfactor      = loadincr
		loadfactorPrev2 = 0.0
		loadfactorPrev  = 0.0

		converged     = false
		convergedPrev = false

		loadStepConverged = 0
	)
	// output = [disp(outputlist)];
	// llist = [0.0];
	//
	// dispFull = [disp];
	//
	for loadStep := 1; loadStep <= maxloadSteps; loadStep++ {
		fmt.Printf("load step = %d\n", loadStep)

		if loadStep > 1 {
			//  Ds
			//  DsPrev
			DsFactor1 := Ds / DsPrev
			for i := range disp {
				disp[i] = (1.0+DsFactor1)*dispPrev[i] - DsFactor1*dispPrev2[i]
			}
			loadfactor = (1.0+DsFactor1)*loadfactorPrev - DsFactor1*loadfactorPrev2
		}
		//
		Du := make([]float64, neq)
		for i := range Du {
			Du[i] = disp[i] - dispPrev[i]
		}
		Dl := loadfactor - loadfactorPrev

		convergedPrev = converged
		converged = false
		//
		//     for iter = 1:10
		for iter := 1; iter <= 10; iter++ {
			//
			var Kglobal [][]float64 = make([][]float64, neq)
			for i := 0; i < neq; i++ {
				Kglobal[i] = make([]float64, neq)
			}
			// Kglobal = zeros(neq,neq);
			Rglobal := make([]float64, neq) // = zeros(neq,1);
			// Kglobal(1:end,1:end) = 0.0;
			// Rglobal(1:end) = 0.0;
			//
			// if(ndim == 2)
			//   if(ndof == 2) % Truss element
			//     for e = 1:nelem
			//         [Klocal, Flocal] = Truss_2D_model1(elemData, elemConn, e, coords, disp, bf);
			//         Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
			//         Rglobal = Assembly_Vector(Rglobal,Flocal,LM,e);
			//     end
			//   else % Beam element
			for e := 0; e < nelem; e++ {
				// for e = 1:nelem
				Klocal, Flocal := GeomExactBeam_2D(elemData, elemConn, e, coords, disp, bf)
				Kglobal = Assembly_Matrix(Kglobal, Klocal, LM, e)
				Rglobal = Assembly_Vector(Rglobal, Flocal, LM, e)
				// end
			}
			//   end
			// else
			//   if(ndof == 3) % Truss element
			//     for e = 1:nelem
			//         [Klocal, Flocal] = Truss_3D_model2(elemData, elemConn, e, coords, disp, bf);
			//         Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
			//         Rglobal = Assembly_Vector(Rglobal,Flocal,LM,e);
			//     end
			//   end
			// end
			//
			for i := range Rglobal {
				Rglobal[i] += loadfactor * Fext[i]
			}
			//
			// %        [converged, du, dl] = solve_arclength(loadStep, neq, iter, Kglobal, Rglobal, dof_force, Fext, assy4r, Du, Dl, Ds);
			converged, du, dl, _ := solve_arclength_split(
				loadStep, neq, iter, Kglobal, Rglobal,
				dof_force, Fext, assy4r, Du, Dl, Ds)
			//
			if converged {
				break
			}
			//
			disp[assy4r] = disp[assy4r] + du
			loadfactor = loadfactor + dl
			//
			Du[assy4r] = Du[assy4r] + du
			Dl = Dl + dl
			//     end
		}
		if converged {
			// %      disp
			if loadStep == 1 {
				{
					// Ds = math.Sqrt(Du'*Du + loadfactor*loadfactor*Fext'*Fext);
					var d float64
					for i := range Du {
						d += Du[i] * Du[i]
					}
					var f float64
					for i := range Fext {
						f += Fext[i] * Fext[i]
					}
					Ds = math.Sqrt(d*d + loadfactor*loadfactor*f*f)
				}
				DsMax = Ds
				DsMin = Ds / 1024.0
			} // end

			loadfactorPrev2 = loadfactorPrev
			loadfactorPrev = loadfactor
			dispPrev2 = dispPrev
			dispPrev = disp

			DsPrev = Ds
			if convergedPrev {
				Ds = math.Min(math.Max(2.0*Ds, DsMin), DsMax)
			} // end

			// dispFull = [dispFull; disp];
			// output = [output disp(outputlist)];
			// llist = [llist; loadfactor];

			// %      plot(abs(output(1,:)), llist,'bx-');
			// %      hold on
			//       plot(abs(output(2,:)), llist,'bs-');
			//       hold on
			//
			// %      linestr = strsplit(fname, ".");
			// %      figname = strcat(linestr(1){:}, "_", num2str(loadStepConverged), ".pdf")
			// %      plot_semicircular_arch(coords, disp, figname);
			// %      plot(abs(output(1,:)), llist, 'bs-'); hold on;
			// %      plot(abs(output(end,:)), llist,'bx-');
			// %      plot(abs(dy)/R, (E*I/R/R)*llist.^(-1),'ko-')
			// %      plot(coords(:,1)+disp(1:3:neq-1), coords(:,2)+disp(2:3:neq-1), 'ko-')
			// %      axis([-150 150 -150 150])
			// %      hold on
			//
			// %      for e=1:nelem
			// %        n1 = elemConn(e,3);
			// %        n2 = elemConn(e,4);
			// %        xx = [coords(n1,1)+disp(ndof*(n1-1)+1) coords(n2,1)+disp(ndof*(n2-1)+1)];
			// %        yy = [coords(n1,2)+disp(ndof*(n1-1)+2) coords(n2,2)+disp(ndof*(n2-1)+2)];
			// %        plot(xx, yy, 'ko-')
			// %        hold on
			// %      endfor
			// %      axis([-0.6 0.6 -3 3])
			// %      hold off
			//
			loadStepConverged = loadStepConverged + 1
		} else {
			if convergedPrev {
				Ds = math.Max(Ds*0.5, DsMin)
			} else {
				Ds = math.Max(Ds*0.25, DsMin)
			} // end
		}
		//     end
	}
	// %    waitforbuttonpress
	// end
	//
	// %plot(abs(output(1,:)), llist,'bx-');
	// %hold on
	// %plot(abs(output(2,:)), llist,'ko-');
	//
	// %plot(t, dy,'k-')
	// %figure(1)
	// %plot(coords(:,1)+disp(1:3:neq-1), coords(:,2)+disp(2:3:neq-1), 'ko-')
	// %figure(2)
	// %plot(dx, llist, 'bs-'); hold on; plot(abs(dy),llist,'ko-')
	// % axis([0,10,-10,10])
	//
	// fileID = fopen('solution.dat','w');
	// for ii=1:size(dispFull,1)
	//     fprintf(fileID,'%12.8f \n', dispFull(ii));
	// end
	// fclose(fileID)
	//
	// fileID = fopen('path.dat','w');
	// for ii=1:size(llist,1)
	//      fprintf(fileID,'%12.8f \t %12.8f \t %12.8f \n', llist(ii), output(1,ii), output(2,ii));
	// %    fprintf(fileID,'%12.8f \t %12.8f \t %12.8f \t %12.8f \n', llist(ii), output(1,ii), output(2,ii), output(3,ii));
	// end
	// fclose(fileID)
}

func assembly_Matrix() {
	// function global_matrix=Assembly_Matrix(global_matrix,local_matrix,LM,e)
	// %%% Assemble global_matrix
	// nen=size(LM,2);
	// for aa=1:nen
	//     mm=LM(e,aa);
	//     if mm~=0
	//         for bb=1:nen
	//             nn=LM(e,bb);
	//             if nn~=0
	//                 global_matrix(mm,nn)=global_matrix(mm,nn)+local_matrix(aa,bb);
	//             end
	//         end
	//     end
	// end
}

func assembly_Vector() {
	// function global_vector=Assembly_Vector(global_vector,local_vector,LM,e)
	// %%% Assemble F_global
	// nen=size(LM,2);
	// for aa=1:nen
	//     mm=LM(e,aa);
	//     if mm~=0
	//         global_vector(mm,1)=global_vector(mm,1)+local_vector(aa,1);
	//     end
	// end
}

func GeomExactBeam_2D(elemData []float64, elemConn [][4]int, e int, coords [][2]float64, soln []float64, bf [2]float64) (Klocal [6][6]float64, Flocal [6]float64) {
	// function [Klocal, Flocal]=GeomExactBeam_2D(elemData, elemConn, e, coords, soln, bf)
	// %%% Shape Function Routine for a 1D Lagrange polynomials
	var (
		af = 1.0

		p      = 1
		nlocal = 2
		ndof   = 3
		nsize  = 6

		rho   = elemData[2]
		A     = elemData[3]
		I     = elemData[4]
		E     = elemData[5]
		nu    = elemData[6]
		kappa = elemData[7]

		G  = E / 2.0 / (1.0 + nu)
		EA = E * A
		EI = E * I
		GA = G * A * kappa
	)
	// 	 Klocal=zeros(nsize,nsize); // % Local stiffness matrix
	// 	 Flocal=zeros(nsize,1);   //% Local load vector

	var x0, y0 [2]float64
	x0[1-1] = coords[elemConn[e][3-1]][1-1]
	y0[1-1] = coords[elemConn[e][3-1]][2-1]
	x0[2-1] = coords[elemConn[e][4-1]][1-1]
	y0[2-1] = coords[elemConn[e][4-1]][2-1]

	var (
		dx = x0[2-1] - x0[1-1]
		dy = y0[2-1] - y0[1-1]
		h0 = math.Sqrt(dx*dx + dy*dy)

		cth0 = dx / h0
		sth0 = dy / h0
	)
	var RotMat [6][6]float64 // =zeros(6,6);

	RotMat[1-1][1-1] = cth0
	RotMat[1-1][2-1] = -sth0
	RotMat[2-1][1-1] = sth0
	RotMat[2-1][2-1] = cth0
	RotMat[3-1][3-1] = 1.0
	RotMat[4-1][4-1] = cth0
	RotMat[4-1][5-1] = -sth0
	RotMat[5-1][4-1] = sth0
	RotMat[5-1][5-1] = cth0
	RotMat[6-1][6-1] = 1.0

	var (
		uxn [2]float64 // [0.0,0.0];
		uzn [2]float64 // [0.0,0.0];
		btn [2]float64 // [0.0,0.0];

		res [3]float64    // zeros(3,1);
		B   [3][6]float64 // =zeros(3,6);
		D   [3][3]float64 // =zeros(3,3);
	)
	uxn[1-1] = soln[ndof*(elemConn[e][3-1]-1)+1]
	uzn[1-1] = soln[ndof*(elemConn[e][3-1]-1)+2]
	btn[1-1] = soln[ndof*(elemConn[e][3-1]-1)+3]

	uxn[2-1] = soln[ndof*(elemConn[e][4-1]-1)+1]
	uzn[2-1] = soln[ndof*(elemConn[e][4-1]-1)+2]
	btn[2-1] = soln[ndof*(elemConn[e][4-1]-1)+3]

	{
		var dummy [6]float64
		vector := [6]float64{uxn[0], uzn[0], btn[0], uxn[1], uzn[1], btn[1]}
		for c := 0; c < 6; c++ {
			for r := 0; r < 6; r++ {
				dummy[r] = RotMat[r][c] * vector[c]
			}
		}
		// dummy = RotMat'*[uxn(1); uzn(1); btn(1); uxn(2); uzn(2); btn(2)];
		uxn[1-1] = dummy[1-1]
		uzn[1-1] = dummy[2-1]
		btn[1-1] = dummy[3-1]
		uxn[2-1] = dummy[4-1]
		uzn[2-1] = dummy[5-1]
		btn[2-1] = dummy[6-1]
	}

	nGP := 1
	gpvec, gwvec := Get_Gauss_points(nGP)

	for gp := 0; gp < nGP; gp++ { // 1:nGP
		// [N,dN_dx,d2N_dx2,J,xcoord]=shape_functions_Lagrange_1D(
		//		elemConn(e,3:end), coords, p, gpvec(gp));
		N, dN_dx, d2N_dx2, J, xcoord := Shape_functions_Lagrange_1D(
			elemConn[e][3:], coords, p, gpvec[gp])

		var (
			ux  = 0.0
			uz  = 0.0
			bt  = 0.0
			dux = 0.0
			duz = 0.0
			dbt = 0.0
		)

		for ii := 0; ii < nlocal; ii++ { //ii=1:nlocal
			ux = ux + uxn[ii]*N[ii]
			uz = uz + uzn[ii]*N[ii]
			bt = bt + btn[ii]*N[ii]
			dux = dux + uxn[ii]*dN_dx[ii]
			duz = duz + uzn[ii]*dN_dx[ii]
			dbt = dbt + btn[ii]*dN_dx[ii]
		} // end

		var (
			sbt = math.Sin(bt)
			cbt = math.Cos(bt)

			// %compute average normal strain, shear strain and curvature

			fact = (1.0+dux)*cbt - duz*sbt

			E = dux + 0.5*(dux*dux+duz*duz)
			G = (1.0+dux)*sbt + duz*cbt
			K = dbt * fact

			//  % compute material response (elastic)

			NF = EA * E // % normal force
			SF = GA * G // % shear force
			BM = EI * K // % bending moment

			// % multiply with volume element

			dvol  = J * gwvec[gp]
			fact1 = dvol * af
		)

		NF = NF * dvol
		SF = SF * dvol
		BM = BM * dvol
		var (
			EAdv = EA * fact1
			GAdv = GA * fact1
			EIdv = EI * fact1
		)

		B[1-1][1-1] = (1.0 + dux) * dN_dx[1-1]
		B[1-1][2-1] = duz * dN_dx[1-1]
		B[1-1][3-1] = 0.0

		B[2-1][1-1] = sbt * dN_dx[1-1]
		B[2-1][2-1] = cbt * dN_dx[1-1]
		B[2-1][3-1] = fact * N[1-1]

		B[3-1][1-1] = dbt * cbt * dN_dx[1-1]
		B[3-1][2-1] = -dbt * sbt * dN_dx[1-1]
		B[3-1][3-1] = fact*dN_dx[1-1] - G*dbt*N[1-1]

		B[1-1][4-1] = (1.0 + dux) * dN_dx[2-1]
		B[1-1][5-1] = duz * dN_dx[2-1]
		B[1-1][6-1] = 0.0

		B[2-1][4-1] = sbt * dN_dx[2-1]
		B[2-1][5-1] = cbt * dN_dx[2-1]
		B[2-1][6-1] = fact * N[2-1]

		B[3-1][4-1] = dbt * cbt * dN_dx[2-1]
		B[3-1][5-1] = -dbt * sbt * dN_dx[2-1]
		B[3-1][6-1] = fact*dN_dx[2-1] - G*dbt*N[2-1]

		D[1-1][1-1] = EAdv
		D[2-1][2-1] = GAdv
		D[3-1][3-1] = EIdv

		res[1-1] = NF
		res[2-1] = SF
		res[3-1] = BM

		// TODO: Klocal = Klocal + ( B'*D*B );
		// TODO: Flocal = Flocal - ( B'*res );

		fact1 = (+SF*cbt - BM*dbt*sbt) * af
		fact2 := (-SF*sbt - BM*dbt*cbt) * af

		for ii := 0; ii < nlocal; ii++ { // ii=1:nlocal
			TI := 3*(ii-1) + 1
			TIp1 := TI + 1
			TIp2 := TI + 2

			for jj := 0; jj < nlocal; jj++ { // jj=1:nlocal
				TJ := 3*(jj-1) + 1
				TJp1 := TJ + 1
				TJp2 := TJ + 2

				Klocal[TI][TJ] = Klocal[TI][TJ] + dN_dx[ii]*NF*dN_dx[jj]*af
				Klocal[TIp1][TJp1] = Klocal[TIp1][TJp1] + dN_dx[ii]*NF*dN_dx[jj]*af

				fact3 := +dN_dx[ii] * BM * cbt * dN_dx[jj] * af
				fact4 := -dN_dx[ii] * BM * sbt * dN_dx[jj] * af

				Klocal[TI][TJp2] = Klocal[TI][TJp2] + (fact3 + dN_dx[ii]*fact1*N[jj])
				Klocal[TIp1][TJp2] = Klocal[TIp1][TJp2] + (fact4 + dN_dx[ii]*fact2*N[jj])
				Klocal[TIp2][TJ] = Klocal[TIp2][TJ] + (fact3 + N[ii]*fact1*dN_dx[jj])
				Klocal[TIp2][TJp1] = Klocal[TIp2][TJp1] + (fact4 + N[ii]*fact2*dN_dx[jj])
				Klocal[TIp2][TJp2] = Klocal[TIp2][TJp2] + (N[ii]*(-SF*G-BM*dbt*fact)*N[jj]-dN_dx[ii]*BM*G*N[jj]-N[ii]*BM*G*dN_dx[jj])*af
			} // end
		} // end
	} //  end

	h := coords[elemConn[e][4-1]] - coords[elemConn[e][3-1]]

	Flocal[1-1] = Flocal[1-1] + 0.5*h*bf[1-1]
	Flocal[4-1] = Flocal[4-1] + 0.5*h*bf[1-1]

	Flocal[2-1] = Flocal[2-1] + 0.5*h*bf[2-1]
	Flocal[5-1] = Flocal[5-1] + 0.5*h*bf[2-1]

	Flocal = RotMat * Flocal
	// TODO : Klocal = RotMat*Klocal*RotMat';

	return
}

func Get_Gauss_points(nGP int) (gp, gw []float64) {
	// function [gp, gw] = get_Gauss_points(nGP)
	switch nGP {
	case 1: // if(nGP == 1)
		gp = []float64{0.0}
		gw = []float64{2.0}
	case 2: // elseif(nGP == 2)
		//     %%% 2 Point quadrature rule
		gp = []float64{-0.577350269, 0.577350269}
		gw = []float64{1.0, 1.0}
	case 3: // elseif(nGP == 3)
		//     %%% 3 Point quadrature rule
		gp = []float64{-math.Sqrt(3 / 5), 0, math.Sqrt(3 / 5)}
		gw = []float64{5 / 9, 8 / 9, 5 / 9}
	case 4: // elseif(nGP == 4)
		//     %%% 4 Point quadrature rule
		gp = []float64{-0.861136311594953, -0.339981043584856, 0.339981043584856, 0.861136311594953}
		gw = []float64{0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454}
	case 5: // else
		//     %%% 5 Point quadrature rule
		gp = []float64{-0.9061798459, -0.5384693101, 0, 0.5384693101, 0.9061798459}
		gw = []float64{0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851}
		// end
		// % create and assemble element matrices
	}
	return
}

func get_global_matrix_vector() {
	// function [Kglobal, Fglobal] = get_global_matrix_vector(elemData, elemConn, coords, LM, td, disp, veloCur, acceCur, bf)
	//
	// neq = max(size(disp))
	// nelem = size(elemConn)(1)
	//
	// Kglobal = zeros(neq,neq);
	// Fglobal = zeros(neq,1);
	//
	// for e = 1:nelem
	//     [Klocal, Flocal] = GeomExactBeam_2D(elemData, elemConn, e, coords, td, disp, veloCur, acceCur, bf);
	//     Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
	//     Fglobal = Assembly_Vector(Fglobal,Flocal,LM,e);
	// end
}

func Processfile() (
	ndim, ndof, nnode, nelem int, coords [][2]float64, elemConn [][4]int,
	elemData []float64,
	LM, neq int,
	assy4r, dof_force, Fext []float64, maxloadSteps int, loadincr float64,
	outputlist []int,
) {
	// function [ndim, ndof, nnode, nelem, coords, elemConn, elemData, LM, neq,\
	// assy4r, dof_force, Fext, maxloadSteps, loadincr, outputlist] =          \
	// processfile(fname)
	//
	// clc
	// % coords: global coordinates of the nodes, x, y, and z
	// % elemConn: element connectivities
	// % nelem: total number of elements
	// % nnode: total number of nodes
	// % nperelem: number of nodes per element for all elements
	// % ndof: number of degrees of per node
	//
	// fid=fopen(fname,'r');
	//
	// % ndim
	//
	// line=fgets(fid);
	// linestr = strsplit(line, ",");
	// ndim = int32(str2num(linestr{1,2}))
	ndim = 2
	//
	// % ndof
	//
	// line=fgets(fid);
	// linestr = strsplit(line, ",");
	// ndof = int32(str2num(linestr{1,2}))
	ndof = 3
	//
	// % nodes
	//
	// line=fgets(fid);
	// linestr = strsplit(line, ",");
	// nnode = int32(str2num(linestr{1,2}))
	nnode = 21
	//
	// nperelem := 2
	// nsize := nperelem * ndof
	neq = nnode * ndof
	// %if(arclen)
	// %  neq = neq+1;
	// %end
	// %neq
	//
	// coords = zeros(nnode,ndim);
	// for i=1:nnode
	//     line = fgets(fid);
	//     linestr = strsplit(line, ",");
	//
	//     coords(i,1) = double(str2num(linestr{1,2}));
	//     coords(i,2) = double(str2num(linestr{1,3}));
	//     if(ndim == 3)
	//       coords(i,3) = double(str2num(linestr{1,4}));
	//     end
	// end
	coords = [][2]float64{
		{0, 0},
		{0, 12},
		{0, 24},
		{0, 36},
		{0, 48},
		{0, 60},
		{0, 72},
		{0, 84},
		{0, 96},
		{0, 108},
		{0, 120},
		{12, 120},
		{24, 120},
		{36, 120},
		{48, 120},
		{60, 120},
		{72, 120},
		{84, 120},
		{96, 120},
		{108, 120},
		{120, 120},
	}
	//
	//
	// % element data
	//
	// line=fgets(fid);
	// linestr = strsplit(line, ",");
	// nelemData = int32(str2num(linestr{1,2}))
	// nelemData := 1
	// elemData = zeros(nelemData, 10);
	// for i=1:nelemData
	//     line = fgets(fid);
	//     linestr = strsplit(line, ",");
	//     for j=1:10
	//       elemData(i,j) = double(str2num(linestr{1,j+1}));
	//     end
	// end
	elemData = []float64{1, 1, 0.0, 6.0, 2.0, 720.0, 0.3, 1.0, 0.0, 0.0, 0.0}
	//
	// % elements
	//
	// line=fgets(fid);
	// linestr = strsplit(line, ",");
	// nelem = int32(str2num(linestr{1,2}))
	nelem = 20
	//
	// elemConn = zeros(nelem, 4, "int32");
	// for i=1:nelem
	//     line = fgets(fid);
	//     linestr = strsplit(line, ",");
	//     elemConn(i,1) = int32(str2num(linestr{1,2}));
	//     elemConn(i,2) = int32(str2num(linestr{1,3}));
	//     elemConn(i,3) = int32(str2num(linestr{1,4}));
	//     elemConn(i,4) = int32(str2num(linestr{1,5}));
	// end
	elemConn = [][4]int{
		{1, 1, 1, 2},
		{1, 1, 2, 3},
		{1, 1, 3, 4},
		{1, 1, 4, 5},
		{1, 1, 5, 6},
		{1, 1, 6, 7},
		{1, 1, 7, 8},
		{1, 1, 8, 9},
		{1, 1, 9, 10},
		{1, 1, 10, 11},
		{1, 1, 11, 12},
		{1, 1, 12, 13},
		{1, 1, 13, 14},
		{1, 1, 14, 15},
		{1, 1, 15, 16},
		{1, 1, 16, 17},
		{1, 1, 17, 18},
		{1, 1, 18, 19},
		{1, 1, 19, 20},
		{1, 1, 20, 21},
	}
	//
	// % Dirichlet boundary conditions
	//
	// line=fgets(fid);
	// linestr = strsplit(line, ",");
	// nDBC    = int32(str2num(linestr{1,2}))

	nDBC := 4

	//
	// %dbclist = zeros(nDBC, 3);
	// dbcnodes = zeros(nDBC, 1, "int32");

	dbcnodes := make([]int, nDBC)

	// for i=1:nDBC
	//     line = fgets(fid);
	//     linestr = strsplit(line, ",");
	//     n1 = int32(str2num(linestr{1,1}));
	//     n2 = int32(str2num(linestr{1,2}));
	// %    dbclist(i,3) = double(str2num(linestr{1,3}));
	//     dbcnodes(i) = (n1-1)*ndof+n2;
	// end

	dbcnodes[0] = (1-1)*ndof + 1
	dbcnodes[1] = (1-1)*ndof + 2
	dbcnodes[2] = (21-1)*ndof + 1
	dbcnodes[3] = (21-1)*ndof + 2

	//
	// assy4r = setdiff([1:neq], dbcnodes)';
	//
	// % Force boundary conditions
	//
	// line=fgets(fid);
	// linestr = strsplit(line, ",");
	// nFBC    = int32(str2num(linestr{1,2}))

	// nFBC := 1

	// fbclist = zeros(nFBC, 3);
	// dof_force = zeros(nFBC, 1, "int32");
	// Fext = zeros(neq, 1);
	Fext = make([]float64, neq)
	// for i=1:nFBC
	//     line = fgets(fid);
	//     linestr = strsplit(line, ",");
	//     n1 = int32(str2num(linestr{1,1}));
	//     n2 = int32(str2num(linestr{1,2}));
	//     ind = (n1-1)*ndof + n2;
	//     dof_force(i) = ind;
	//     Fext(ind) = double(str2num(linestr{1,3}));
	// end
	Fext[(13-1)*ndof+2] = -1.0
	//
	// % for output
	//
	// line=fgets(fid);
	// linestr = strsplit(line, ",");
	// nOutput = int32(str2num(linestr{1,2}))
	// nOutput := 2
	// outputlist = zeros(nOutput, 1, "int32");
	outputlist = []int{
		(13-1)*ndof + 1,
		(13-1)*ndof + 2,
	}
	// for i=1:nOutput
	//     line = fgets(fid);
	//     linestr = strsplit(line, ",");
	//     n1 = int32(str2num(linestr{1,1}));
	//     n2 = int32(str2num(linestr{1,2}));
	//     outputlist(i,1) = (n1-1)*ndof+n2;
	// end
	//
	//
	// % Arclength parameters
	//
	// line=fgets(fid);
	// linestr = strsplit(line, ",");
	// arclen  = (int32(linestr{1,2}) == 1);
	// arclen := true
	//
	// line=fgets(fid)
	// linestr = strsplit(line, ",");
	// maxloadSteps = int32(str2num(linestr{1,1}));
	maxloadSteps = 50
	//
	// line=fgets(fid)
	// linestr = strsplit(line, ",");
	// loadincr = double(str2num(linestr{1,1}));
	loadincr = 0.5
	//
	// fclose(fid);
	//
	// % data structures
	//
	// LM  = zeros(nelem, nsize);
	//
	// for e=1:nelem
	//     count = 1;
	//     for jj=1:nperelem
	//         ind = ndof*(elemConn(e,jj+2)-1)+1;
	//         for kk=1:ndof
	//             LM(e,count) = ind;
	//             ind = ind + 1;
	//             count = count + 1;
	//         end
	//     end
	//     count = count - 1;
	// end
	return
}

func Shape_functions_Lagrange_1D(NodeNums []int, coords [][]float64, p int, xi float64) (
	N, dN_dx, d2N_dx2 []float64, J, x float64,
) {
	// function [N,dN_dx,d2N_dx2,J,x]=shape_functions_Lagrange_1D(NodeNums, coords, p, xi)
	// %%% Shape Function Routine for a 1D Lagrange polynomials
	// nen = p+1;
	nen := p + 1
	// N =zeros(nen,1);
	N = make([]float64, nen)
	// dN_dxi   = zeros(nen,1);
	dN_dxi := make([]float64, nen)
	// d2N_dxi2 = zeros(nen,1);
	d2N_dxi2 := make([]float64, nen)

	switch p {
	case 1: // if(p == 1)
		N[0] = 0.5 * (1.0 - xi) // (1)
		N[1] = 0.5 * (1.0 + xi) // (2)

		dN_dxi[0] = -0.5 // (1)
		dN_dxi[1] = 0.5  // (2)

		d2N_dxi2[0] = 0.0 // (1)
		d2N_dxi2[1] = 0.0 // (2)

	case 2: // elseif(p == 2)
		val := xi * xi
		N[0] = -0.5 * (xi - val) // (1)
		N[1] = 1.0 - val         // (2)
		N[2] = 0.5 * (xi + val)  // (3)

		val = 2.0 * xi
		dN_dxi[0] = -0.5 * (1.0 - val) // (1)
		dN_dxi[1] = -val               // (2)
		dN_dxi[2] = 0.5 * (1.0 + val)  // (3)

		d2N_dxi2[0] = 1.0  // (1)
		d2N_dxi2[1] = -2.0 // (2)
		d2N_dxi2[2] = 1.0  // (3)

	case 3: // elseif(p == 3)
		fact1 := 9. / 16.
		fact2 := 27. / 16.
		val1 := xi * xi

		N[0] = -fact1 * (1 - xi) * (1/9 - val1) // (1)
		N[1] = fact2 * (1 - val1) * (1/3 - xi)  // (2)
		N[2] = fact2 * (1 - val1) * (1/3 + xi)  // (3)
		N[3] = -fact1 * (1 + xi) * (1/9 - val1) // (4)

		val2 := 3.0 * val1
		dN_dxi[0] = -fact1 * (-1/9 - 2.0*xi + val2) // (1)
		dN_dxi[1] = fact2 * (-1 - 2.0/3*xi + val2)  // (2)
		dN_dxi[2] = fact2 * (1 - 2.0/3*xi - val2)   // (3)
		dN_dxi[3] = -fact1 * (1/9 - 2.0*xi - val2)  // (4)

		val2 = 6.0 * xi
		d2N_dxi2[0] = -fact1 * (-2 + val2)  // (1)
		d2N_dxi2[1] = fact2 * (-2/3 + val2) // (2)
		d2N_dxi2[2] = fact2 * (-2/3 - val2) // (3)
		d2N_dxi2[3] = -fact1 * (-2 - val2)  // (4)

	case 4: // elseif(p == 4)
		var (
			fact1 = 2. / 3.
			fact2 = 8. / 3.
			val1  = xi * xi
			val2  = val1 * xi
			val3  = val2 * xi
		)
		N[0] = fact1 * (0.25*xi - 0.25*val1 - val2 + val3)  // (1)
		N[1] = -fact2 * (0.5*xi - val1 - 0.5*val2 + val3)   // (2)
		N[2] = 4.0 * (0.25 - 1.25*val1 - 0 + val3)          // (3)
		N[3] = fact2 * (0.5*xi + val1 - 0.5*val2 - val3)    // (4)
		N[4] = -fact1 * (0.25*xi + 0.25*val1 - val2 - val3) // (5)

		val4 := 4.0 * val2
		dN_dxi[0] = fact1 * (0.25 - 0.5*xi - 3.0*val1 + val4)  // (1)
		dN_dxi[1] = -fact2 * (0.5 - 2.0*xi - 1.5*val1 + val4)  // (2)
		dN_dxi[2] = 4.0 * (0 - 2.5*xi - 0.0 + val4)            // (3)
		dN_dxi[3] = fact2 * (0.5 + 2.0*xi - 1.5*val1 - val4)   // (4)
		dN_dxi[4] = -fact1 * (0.25 + 0.5*xi - 3.0*val1 - val4) // (5)

		val4 = 12.0 * val1
		d2N_dxi2[0] = fact1 * (-0.5 - 6.0*xi + val4)  // (1)
		d2N_dxi2[1] = -fact2 * (-2.0 - 3.0*xi + val4) // (2)
		d2N_dxi2[2] = 4.0 * (-2.5 - 0.0 + val4)       // (3)
		d2N_dxi2[3] = fact2 * (2.0 - 3.0*xi - val4)   // (4)
		d2N_dxi2[4] = -fact1 * (0.5 - 6.0*xi - val4)  // (5)
	default:
		// else
		//     sprintf('
		panic("no basis functions defined for this degree")
		// end
	}

	Jx := 0.0
	Jy := 0.0
	x = 0.0
	for kk := 0; kk < nen; kk++ { //1:nen
		Jx = Jx + dN_dxi[kk]*coords[NodeNums[kk]][0] //1]
		Jy = Jy + dN_dxi[kk]*coords[NodeNums[kk]][1] // 2]
		x = x + N[kk]*coords[NodeNums[kk]][0]        // TODO not clear?
	} // end
	J = math.Sqrt(Jx*Jx + Jy*Jy)

	for i := range dN_dxi {
		dN_dx[i] = (1.0 / J) * dN_dxi[i]
	}
	for i := range d2N_dxi2 {
		d2N_dx2[i] = math.Pow(1.0/J, 2) * d2N_dxi2[i]
	}
	return
}

// func solve_arclength() {
// function  [converged, du, dl] = solve_arclength(timeStep, neq, iter, \
//                 Kglobal, Fglobal, dof_force, Fext, assy4r, Du, Dl, ds)
//
//     if(timeStep > 1)
//         dummy = Du'*Du + Dl*Dl - ds*ds;
//         Fglobal(neq) = Fglobal(neq) - dummy;
//
//         Kglobal(:, neq)   = Kglobal(:, neq) - Fext;
//
//         Kglobal(neq, 1:neq-1) = Kglobal(neq, 1:neq-1) + 2.0*Du';
//
//         Kglobal(neq, neq)  = Kglobal(neq, neq)  + 2.0*Dl;
//     else
//         Fglobal(neq) = 0.0 ;
//         Kglobal(neq, neq) = 1.0;
//     endif
//
//     %%% Applying Boundary Conditions
//
//     K1 = Kglobal(assy4r,assy4r);
//     F1 = Fglobal(assy4r);
//
//     rNorm = norm(F1,2);
//
//     printf(' rNorm : %5d ...  %12.6E \n', iter, rNorm);
//     du = F1*0.0;
//     dl = 0.0;
//     converged = false;
//
//     if(rNorm < 1.0e-8)
//        converged = true;
//        return;
//     end
//
//     %% solve the matrix system
//     dutemp = K1\F1;
//
//     du = dutemp(1:end-1);
//     dl = dutemp(end);
// end
// }

func solve_arclength_split(timeStep, neq, iter int,
	Kglobal [][]float64, Rglobal [][]float64,
	dof_force, Fext, assy4r, Du, Dl, ds, du1 float64) (converged bool, du []float64, dl, du1 float64) {
	// function[converged, du, dl, du1] = solve_arclength_split(timeStep, neq, iter, Kglobal, Rglobal, dof_force, Fext, assy4r, Du, Dl, ds, du1)
	//
	//     psi = 1.0;
	//
	//     FextReduced = Fext(assy4r);
	//     FtF = Fext'*Fext;
	//     if(timeStep > 1)
	//         A = Du'*Du + psi*Dl*Dl*FtF - ds*ds;
	//         a = 2.0*Du(assy4r)';
	//         b = 2.0*psi*Dl*FtF;
	//     else
	//         A = 0.0;
	//         a = 0.0*Du(assy4r)';
	//         b = 1.0;
	//     end
	//
	//     %%% Applying Boundary Conditions
	//
	//     R = Rglobal(assy4r);
	//
	//     rNorm = norm(R,2);
	//     rNorm = math.Sqrt(rNorm*rNorm + A*A);
	//
	//     fprintf(' rNorm : %5d ...  %12.6E \n', iter, rNorm);
	//     du = R*0.0;
	//     dl = 0.0;
	//     converged = false;
	//
	//     if(rNorm < 1.0e-6)
	//        converged = true;
	//        return;
	//     end
	//
	//     K1 = Kglobal(assy4r,assy4r);
	//     [L, U, P] = lu(sparse(K1));
	//
	//     %% solve the matrix system
	//     duu = L\(P*FextReduced);
	//     du1 = U\duu;
	//     duu = L\(P*R);
	//     du2 = U\duu;
	//     du2 = -du2; % this is because the Residual is added to the RHS
	//
	//     dl = (a*du2 - A)/(b+a*du1);
	//
	//     du = -du2 + dl*du1;
	// end
	return
}

// TODO : NO NEED FUNCTION
// func timeSteppingParameters_Solid() {
// function td = timeSteppingParameters_Solid(tis, rho, dt)
//
// td = zeros(100,1);
//
// switch tis
//    case 0
//      alpf = 1.0;
//      alpm = 1.0;
//      gamm = 1.0;
//
//      td(1)  = alpm;
//      td(2)  = alpf;
//      td(3)  = alpm;
//      td(4)  = gamm;
//      td(7)  = alpf;
//
//             % velocity is used as the primary variable for
//             % solid dynamics problem
//             % td[10] is the multiplication when converting from
//             % displacement based formulation to velocity based formulation
//             % It is set to ONE for static problem
//
//       td(10) = 1.0;
//
//     case 2 % Backward-Euler
//
//        alpf = 1.0;
//        alpm = 1.0;
//        gamm = 1.0;
//
//        td(1)  = alpm;
//        td(2)  = alpf;
//        td(3)  = alpm;
//        td(4)  = gamm;
//
//        td(5)  = 1.0/dt/dt;
//        td(6)  = 1.0/dt;
//        td(7)  = alpf;
//
//        %displacement as the primary variable
//        %v_{n+1}  = td(10)*d_{n+1} + td(11)*d_n + td(12)*v_n + td(13)*a_n + td(14)*ddot_n;
//        %a_{n+1}  = td(15)*d_{n+1} + td(16)*d_n + td(17)*v_n + td(18)*a_n + td(19)*ddot_n;
//
//        td(10) = 1.0/dt;    % d_{n+1}
//        td(11) = -td(10);   % d_n
//        td(12) = 0.0;       % v_n
//        td(13) = 0.0;       % a_n
//        td(14) = 0.0;       % ddot_n
//
//        td(15) = 1.0/dt/dt; % d_{n+1}
//        td(16) = -td(15);   % d_n
//        td(17) = -1.0/dt;   % vn
//        td(18) = 0.0;       % an
//        td(19) = 0.0;       % ddot_n
//
//        %velocity as the primary variable
//        %d_{n+1}  = td(20)*v_{n+1} + td(21)*d_n + td(22)*v_n + td(23)*a_n + td(24)*ddot_n ;
//        %a_{n+1}  = td(25)*v_{n+1} + td(26)*d_n + td(27)*v_n + td(28)*a_n + td(29)*ddot_n ;
//
//        td(40) = dt;        % v_{n+1}
//        td(41) = 1.0;       % d_n
//        td(42) = 0.0;       % v_n
//        td(43) = 0.0;       % a_n
//        td(44) = 0.0;       % ddot_n
//
//        td(45) = 1.0/dt;    % v_{n+1}
//        td(46) = 0.0;       % d_n
//        td(47) = -td(45);   % v_n
//        td(48) = 0.0;       % a_n
//        td(49) = 0.0;       % ddot_n
//
//    case 3 % CH-aplha
//
//        alpm = (2.0-rho)/(rho+1.0);
//        alpf = 1.0/(rho+1.0);
//
//        gamm = 0.5 + alpm - alpf;
//        beta = 0.25*(1.0+alpm-alpf)*(1.0+alpm-alpf);
//
//        td(1)  = alpm;
//        td(2)  = alpf;
//        td(3)  = alpm;
//        td(4)  = gamm;
//
//        td(5)  = alpm/beta/dt/dt;
//        td(6)  = alpf*gamm/beta/dt;
//        td(7)  = alpf;
//
//        %displacement as the primary variable
//        %v_{n+1}  = td(10)*d_{n+1} + td(11)*d_n + td(12)*v_n + td(13)*a_n + td(14)*ddot_n;
//        %a_{n+1}  = td(15)*d_{n+1} + td(16)*d_n + td(17)*v_n + td(18)*a_n + td(19)*ddot_n;
//
//        td(10) = gamm/beta/dt;           % d_{n+1}
//        td(11) = -td(10);                % d_n
//        td(12) = 1.0-gamm/beta;          % v_n
//        td(13) = dt*(1.0-gamm/2.0/beta); % a_n
//        td(14) = 0.0;                    % ddot_n
//
//        td(15) = 1.0/beta/dt/dt;         % d_{n+1}
//        td(16) = -td(15);                % d_n
//        td(17) = -1.0/beta/dt;           % v_n
//        td(18) = 1.0-1.0/2.0/beta;       % a_n
//        td(19) = 0.0;                    % ddot_n
//
//        %velocity as the primary variable
//        %d_{n+1}  = td(20)*v_{n+1} + td(21)*d_n + td(22)*v_n + td(23)*a_n + td(24)*ddot_n ;
//        %a_{n+1}  = td(25)*v_{n+1} + td(26)*d_n + td(27)*v_n + td(28)*a_n + td(29)*ddot_n ;
//
//        td(40) = dt*beta/gamm;                     % v_{n+1}
//        td(41) = 1.0;                              % d_n
//        td(42) = dt*(gamm-beta)/gamm;              % v_n
//        td(43) = dt*dt*(gamm-2.0*beta)/(2.0*gamm); % a_n
//        td(44) = 0.0;                              % ddot_n
//
//        td(45) = 1.0/(gamm*dt);                    % v_{n+1}
//        td(46) = 0.0;                              % d_n
//        td(47) = -td(45);                          % v_n
//        td(48) = (gamm-1.0)/gamm;                  % a_n
//        td(49) = 0.0;                              % ddot_n
//
//    case 4 % JHW-alpha or KDP-alpha
//
//        alpf = 1.0/(1.0 + rho);
//        alpm = 0.5*(3.0 - rho)/(1.0 + rho);
//
//        gamm = 0.5 + alpm - alpf;
//
//        td(1)  = alpm;
//        td(2)  = alpf;
//        td(3)  = alpm;
//        td(4)  = gamm;
//
//        td(5)  = (alpm*alpm)/(alpf*gamm*gamm*dt*dt);
//        td(6)  = alpm/gamm/dt;
//        td(7)  = alpf;
//
//        %displacement as the primary variable
//        %v_{n+1}    = td(10)*d_{n+1} + td(11)*d_n + td(12)*v_n + td(13)*a_n + td(14)*ddot_n;
//        %a_{n+1}    = td(15)*d_{n+1} + td(16)*d_n + td(17)*v_n + td(18)*a_n + td(19)*ddot_n;
//        %ddot_{n+1} = td(20)*d_{n+1} + td(21)*d_n + td(22)*v_n + td(23)*a_n + td(24)*ddot_n;
//
//        td(10) = alpm/(alpf*gamm*dt);             % d_{n+1}
//        td(11) = -td(10);                         % d_n
//        td(12) = (alpf-1.0)/alpf;                 % v_n
//        td(13) = 0.0;                             % a_n
//        td(14) = (gamm-alpm)/alpf/gamm;           % ddot_n
//
//        td(15) = alpm/(alpf*gamm*gamm*dt*dt);     % d_{n+1}
//        td(16) = -td(15);                         % d_n
//        td(17) = -1.0/(alpf*gamm*dt);             % v_n
//        td(18) = (gamm-1.0)/gamm;                 % a_n
//        td(19) = (gamm-alpm)/(alpf*gamm*gamm*dt); % ddot_n
//
//        td(20) = 1.0/(gamm*dt);                   % d_{n+1}
//        td(21) = -td(20);                         % d_n
//        td(22) = 0.0;                             % v_n
//        td(23) = 0.0;                             % a_n
//        td(24) = (gamm-1.0)/gamm;                 % ddot_n
//
//        %velocity as the primary variable
//        %d_{n+1}    = td(20)*v_{n+1} + td(21)*d_n + td(22)*v_n + td(23)*a_n + td(24)*ddot_n ;
//        %a_{n+1}    = td(25)*v_{n+1} + td(26)*d_n + td(27)*v_n + td(28)*a_n + td(29)*ddot_n ;
//        %ddot_{n+1} = td(30)*v_{n+1} + td(31)*d_n + td(32)*v_n + td(33)*a_n + td(34)*ddot_n;
//
//        td(40) = alpf*gamm*dt/alpm;           % v_{n+1}
//        td(41) = 1.0;                         % d_n
//        td(42) = (1.0-alpf)*gamm*dt/alpm;     % v_n
//        td(43) = 0.0;                         % a_n
//        td(44) = (alpm-gamm)*dt/alpm;        % ddot_n
//
//        td(45) = 1.0/gamm/dt;                 % v_{n+1}
//        td(46) = 0.0;                         % d_n
//        td(47) = -td(45);                     % v_n
//        td(48) = 1.0-1.0/gamm;                % v_n
//        td(49) = 0.0;                         % ddot_n
//
//        td(50) = alpf/alpm;                   % v_{n+1}
//        td(51) = 0.0;                         % d_n
//        td(52) = (1.0-alpf)/alpm;             % v_n
//        td(53) = 0.0;                         % a_n
//        td(54) = (alpm-1.0)/alpm;             % ddot_n
//
//    case 5 %% Newmark-beta method
//
//        alpm = 1.0;
//        alpf = 1.0;
//
//        td(1)  = alpm;
//        td(2)  = alpf;
//        td(3)  = alpm;
//        td(4)  = gamm;
//
//        td(5)  = alpm/beta/dt/dt;
//        td(6)  = alpf*gamm/beta/dt;
//        td(7)  = alpf;
//
//        %displacement as the primary variable
//        %v_{n+1}  = td(10)*d_{n+1} + td(11)*d_n + td(12)*v_n + td(13)*a_n + td(14)*ddot_n;
//        %a_{n+1}  = td(15)*d_{n+1} + td(16)*d_n + td(17)*v_n + td(18)*a_n + td(19)*ddot_n;
//
//        td(10) = gamm/beta/dt;           % d_{n+1}
//        td(11) = -td(10);                % d_n
//        td(12) = 1.0-gamm/beta;          % v_n
//        td(13) = dt*(1.0-gamm/2.0/beta); % a_n
//        td(14) = 0.0;                    % ddot_n
//
//        td(15) = 1.0/beta/dt/dt;         % d_{n+1}
//        td(16) = -td(15);                % d_n
//        td(17) = -1.0/beta/dt;           % v_n
//        td(18) = 1.0-1.0/2.0/beta;       % a_n
//        td(19) = 0.0;                    % ddot_n
//
//        %velocity as the primary variable
//        %d_{n+1}  = td(20)*v_{n+1} + td(21)*d_n + td(22)*v_n + td(23)*a_n + td(24)*ddot_n ;
//        %a_{n+1}  = td(25)*v_{n+1} + td(26)*d_n + td(27)*v_n + td(28)*a_n + td(29)*ddot_n ;
//
//        td(40) = dt*beta/gamm;                     % v_{n+1}
//        td(41) = 1.0;                              % d_n
//        td(42) = dt*(gamm-beta)/gamm;              % v_n
//        td(43) = dt*dt*(gamm-2.0*beta)/(2.0*gamm); % a_n
//        td(44) = 0.0;                              % ddot_n
//
//        td(45) = 1.0/(gamm*dt);                    % v_{n+1}
//        td(46) = 0.0;                              % d_n
//        td(47) = -td(45);                          % v_n
//        td(48) = (gamm-1.0)/gamm;                  % a_n
//        td(49) = 0.0;                              % ddot_n
//
//     otherwise
//        printf('tis error\n');
//        printf('This time integration scheme is not available yet!');
// end
// }

func truss_2D_model1() {
	// function [Klocal, Flocal]=Truss_2D_model1(elemData, elemConn, e, coords, soln, bf)
	//
	// ndof = 2;
	//
	// finite = (int32(elemData(elemConn(e,1), 1)) == 1) ;
	// rho0 = elemData(elemConn(e,1),2);
	// A0   = elemData(elemConn(e,1),3);
	// E    = elemData(elemConn(e,1),4);
	//
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//
	// disp  = zeros(4,1);
	//
	// X1 = coords(elemConn(e,3),1);
	// Y1 = coords(elemConn(e,3),2);
	// X2 = coords(elemConn(e,4),1);
	// Y2 = coords(elemConn(e,4),2);
	//
	// disp(1) = soln(ndof*(elemConn(e,3)-1)+1);
	// disp(2) = soln(ndof*(elemConn(e,3)-1)+2);
	// disp(3) = soln(ndof*(elemConn(e,4)-1)+1);
	// disp(4) = soln(ndof*(elemConn(e,4)-1)+2);
	//
	//
	// % compute the orientation of the element
	//
	//     dx = X2 - X1;
	//     dy = Y2 - Y1;
	//
	//     L0 = math.Sqrt(dx*dx+dy*dy);
	//
	//     x1 = X1 + disp(1);
	//     y1 = Y1 + disp(2);
	//     x2 = X2 + disp(3);
	//     y2 = Y2 + disp(4);
	//
	//     dx = x2 - x1;
	//     dy = y2 - y1;
	//
	//     L = math.Sqrt(dx*dx+dy*dy);
	//
	//     B0 = zeros(4,1);
	//     B(1) = -dx;   B(2) = -dy;   B(3) =  dx;   B(4) =  dy;
	//
	//     H = zeros(4,4);
	//     H(1,1) =  1.0;    H(1,3) = -1.0;
	//     H(2,2) =  1.0;    H(2,4) = -1.0;
	//     H(3,1) = -1.0;    H(3,3) =  1.0;
	//     H(4,2) = -1.0;    H(4,4) =  1.0;
	//
	//     Klocal = zeros(4,4);
	//     Flocal = zeros(4,1);
	//
	//     strain = (L/L0 - 1.0) ;
	//
	//     % axial force
	//     N = (A0*E)*strain;
	//
	//     % residual
	//     Flocal = Flocal - N*B'/L;
	//
	//     % stiffness
	//     Klocal = Klocal + ((E*A0*B')*B)/L^3 ;
	//     Klocal = Klocal + ((N/L)*H);
	// end
}

// func truss_2D_model2() {
// function [Klocal, Flocal]=Truss_2D_model2(elemData, elemConn, e, coords, soln, bf)
//
// ndof = 2;
//
// finite = (int32(elemData(elemConn(e,1), 1)) == 1) ;
// rho0 = elemData(elemConn(e,1),2);
// A0   = elemData(elemConn(e,1),3);
// E    = elemData(elemConn(e,1),4);
//
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// disp  = zeros(4,1);
//
// X1 = coords(elemConn(e,3),1);
// Y1 = coords(elemConn(e,3),2);
// X2 = coords(elemConn(e,4),1);
// Y2 = coords(elemConn(e,4),2);
//
// disp(1) = soln(ndof*(elemConn(e,3)-1)+1);
// disp(2) = soln(ndof*(elemConn(e,3)-1)+2);
// disp(3) = soln(ndof*(elemConn(e,4)-1)+1);
// disp(4) = soln(ndof*(elemConn(e,4)-1)+2);
//
// % compute the orientation of the element
//
//     dx = X2 - X1;
//     dy = Y2 - Y1;
//
//     L0 = math.Sqrt(dx*dx+dy*dy);
//
//     x1 = X1 + disp(1);
//     y1 = Y1 + disp(2);
//     x2 = X2 + disp(3);
//     y2 = Y2 + disp(4);
//
//     dx = x2 - x1;
//     dy = y2 - y1;
//
//     L = math.Sqrt(dx*dx+dy*dy);
//
//     B0 = zeros(4,1);
//     B(1) = -dx;   B(2) = -dy;   B(3) =  dx;   B(4) =  dy;
//     B = B/L0/L0;
//
//     H = zeros(4,4);
//     H(1,1) =  1.0;    H(1,3) = -1.0;
//     H(2,2) =  1.0;    H(2,4) = -1.0;
//     H(3,1) = -1.0;    H(3,3) =  1.0;
//     H(4,2) = -1.0;    H(4,4) =  1.0;
//
//     Klocal = zeros(4,4);
//     Flocal = zeros(4,1);
//
//     strain = (L*L/L0/L0 - 1.0)/2.0 ;
//
//     % axial force
//     N = (A0*E)*strain;
//
//     % residual
//     Flocal = Flocal - N*L0*B';
//
//     % stiffness
//     Klocal = Klocal + ((E*A0*L0*B')*B) ;
//     Klocal = Klocal + ((N/L0)*H);
// }

// func truss_3D_model2() {
// function [Klocal, Flocal]=Truss_3D_model2(elemData, elemConn, e, coords, soln, bf)
//
// ndof = 3;
//
// finite = (int32(elemData(1)) == 1) ;
// rho0 = elemData(2);
// A0   = elemData(3);
// E    = elemData(4);
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// disp  = zeros(6,1);
//
// % rotate nodal displacements and compute nodal positions on element axis
//
// X1 = coords(elemConn(e,3),1);
// Y1 = coords(elemConn(e,3),2);
// Z1 = coords(elemConn(e,3),3);
// X2 = coords(elemConn(e,4),1);
// Y2 = coords(elemConn(e,4),2);
// Z2 = coords(elemConn(e,4),3);
//
// disp(1) = soln(ndof*(elemConn(e,3)-1)+1);
// disp(2) = soln(ndof*(elemConn(e,3)-1)+2);
// disp(3) = soln(ndof*(elemConn(e,3)-1)+3);
// disp(4) = soln(ndof*(elemConn(e,4)-1)+1);
// disp(5) = soln(ndof*(elemConn(e,4)-1)+2);
// disp(6) = soln(ndof*(elemConn(e,4)-1)+3);
//
// % compute the orientation of the element
//
//     dx = X2 - X1;
//     dy = Y2 - Y1;
//     dz = Z2 - Z1;
//
//     L0 = math.Sqrt(dx*dx+dy*dy+dz*dz);
//
//     x1 = X1 + disp(1);
//     y1 = Y1 + disp(2);
//     z1 = Z1 + disp(3);
//
//     x2 = X2 + disp(4);
//     y2 = Y2 + disp(5);
//     z2 = Z2 + disp(6);
//
//     dx = x2 - x1;
//     dy = y2 - y1;
//     dz = z2 - z1;
//
//     L = math.Sqrt(dx*dx+dy*dy+dz*dz);
//
//     B0 = zeros(6,1);
//     B(1) = -dx;   B(2) = -dy;   B(3) = -dz;
//     B(4) =  dx;   B(5) =  dy;   B(6) =  dz;
//     B = B/L0/L0;
//
//     H = zeros(6,6);
//     H(1,1) =  1.0;    H(1,4) = -1.0;
//     H(2,2) =  1.0;    H(2,5) = -1.0;
//     H(3,3) =  1.0;    H(3,6) = -1.0;
//
//     H(4,1) = -1.0;    H(4,4) =  1.0;
//     H(5,2) = -1.0;    H(5,5) =  1.0;
//     H(6,3) = -1.0;    H(6,6) =  1.0;
//
//     Klocal = zeros(6,6);
//     Flocal = zeros(6,1);
//
//     strain = (L*L/L0/L0 - 1.0)/2.0 ;
//
//     % axial force
//     N = (A0*E)*strain;
//
//     % residual
//     Flocal = Flocal - N*L0*B';
//
//     % stiffness
//     Klocal = Klocal + ((E*A0*L0*B')*B) ;
//     Klocal = Klocal + ((N/L0)*H);
//
// end
// }
