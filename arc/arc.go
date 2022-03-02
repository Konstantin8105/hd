package main

import "math"

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
	// [ndim, ndof, nnode, nelem, coords, elemConn, elemData, LM, neq, assy4r, dof_force, Fext, maxloadSteps, loadincr, outputlist] = processfile(fname)
	//
	// disp = zeros(neq,1);
	//
	// dispPrev  = disp;
	// dispPrev2 = disp;
	// dispPrev3 = disp;
	// dispPrev4 = disp;
	//
	// Kglobal = zeros(neq,neq);
	// Rglobal = zeros(neq,1);
	//
	// bf=[0.0 0.0];
	//
	// Ds = loadincr;
	// DsPrev = Ds;
	// DsMax = Ds;
	// DsMin = Ds;
	//
	// loadfactor      = loadincr;
	// loadfactorPrev2 = 0.0;
	// loadfactorPrev  = 0.0;
	//
	// converged = false;
	// convergedPrev = false;
	//
	// loadStepConverged = 0;
	// output = [disp(outputlist)];
	// llist = [0.0];
	//
	// dispFull = [disp];
	//
	// for  loadStep=1:maxloadSteps
	//     fprintf("load step = %d \n", loadStep);
	//
	//     if(loadStep > 1)
	//       Ds
	//       DsPrev
	//       DsFactor1 = Ds/DsPrev
	//       disp     = (1.0+DsFactor1)*dispPrev - DsFactor1*dispPrev2;
	//       loadfactor = (1.0+DsFactor1)*loadfactorPrev - DsFactor1*loadfactorPrev2;
	//     end
	//
	//     Du = disp - dispPrev;
	//     Dl = loadfactor - loadfactorPrev;
	//
	//     convergedPrev = converged;
	//     converged = false;
	//
	//     for iter = 1:10
	//         Kglobal(1:end,1:end) = 0.0;
	//         Rglobal(1:end) = 0.0;
	//
	//         if(ndim == 2)
	//           if(ndof == 2) % Truss element
	//             for e = 1:nelem
	//                 [Klocal, Flocal] = Truss_2D_model1(elemData, elemConn, e, coords, disp, bf);
	//                 Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
	//                 Rglobal = Assembly_Vector(Rglobal,Flocal,LM,e);
	//             end
	//           else % Beam element
	//             for e = 1:nelem
	//                 [Klocal, Flocal] = GeomExactBeam_2D(elemData, elemConn, e, coords, disp, bf);
	//                 Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
	//                 Rglobal = Assembly_Vector(Rglobal,Flocal,LM,e);
	//             end
	//           end
	//         else
	//           if(ndof == 3) % Truss element
	//             for e = 1:nelem
	//                 [Klocal, Flocal] = Truss_3D_model2(elemData, elemConn, e, coords, disp, bf);
	//                 Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
	//                 Rglobal = Assembly_Vector(Rglobal,Flocal,LM,e);
	//             end
	//           end
	//         end
	//
	//         Rglobal = Rglobal + loadfactor*Fext;
	//
	// %        [converged, du, dl] = solve_arclength(loadStep, neq, iter, Kglobal, Rglobal, dof_force, Fext, assy4r, Du, Dl, Ds);
	//         [converged, du, dl] = solve_arclength_split(loadStep, neq, iter, Kglobal, Rglobal, dof_force, Fext, assy4r, Du, Dl, Ds);
	//
	//         if(converged)
	//           break;
	//         end
	//
	//         disp(assy4r) = disp(assy4r) + du;
	//         loadfactor = loadfactor + dl;
	//
	//         Du(assy4r) = Du(assy4r) + du;
	//         Dl = Dl + dl;
	//     end
	//
	//     if (converged)
	// %      disp
	//       if(loadStep == 1)
	//          Ds = math.Sqrt(Du'*Du + loadfactor*loadfactor*Fext'*Fext);
	//          DsMax = Ds;
	//          DsMin = Ds/1024.0;
	//       end
	//
	//       loadfactorPrev2 = loadfactorPrev;
	//       loadfactorPrev  = loadfactor;
	//       dispPrev2 = dispPrev;
	//       dispPrev  = disp;
	//
	//       DsPrev = Ds;
	//       if(convergedPrev)
	//         Ds = min(max(2.0*Ds, DsMin), DsMax);
	//       end
	//
	//       dispFull = [dispFull; disp];
	//       output = [output disp(outputlist)];
	//       llist = [llist; loadfactor];
	//
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
	//       loadStepConverged = loadStepConverged + 1;
	//     else
	//       if(convergedPrev)
	//         Ds = max(Ds*0.5, DsMin);
	//       else
	//         Ds = max(Ds*0.25, DsMin);
	//       end
	//     end
	//
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

func geomExactBeam_2D() {
	// function [Klocal, Flocal]=GeomExactBeam_2D(elmDat, IEN, e, XX, soln, bf)
	// %%% Shape Function Routine for a 1D Lagrange polynomials
	// af = 1.0;
	//
	// p = 1;
	// nlocal = 2;
	// ndof   = 3;
	// nsize  = 6;
	//
	// rho  = elmDat(2);
	// A    = elmDat(3);
	// I    = elmDat(4);
	// E    = elmDat(5);
	// nu   = elmDat(6);
	// kappa= elmDat(7);
	//
	// G  = E/2.0/(1.0+nu);
	// EA = E*A;
	// EI = E*I;
	// GA = G*A*kappa;
	//
	//
	// Klocal=zeros(nsize,nsize); % Local stiffness matrix
	// Flocal=zeros(nsize,1);   % Local load vector
	//
	// x0(1) = XX(IEN(e,3),1);
	// y0(1) = XX(IEN(e,3),2);
	// x0(2) = XX(IEN(e,4),1);
	// y0(2) = XX(IEN(e,4),2);
	//
	// dx = x0(2) - x0(1);
	// dy = y0(2) - y0(1);
	// h0 = math.Sqrt(dx*dx+dy*dy);
	//
	// cth0 = dx/h0;
	// sth0 = dy/h0;
	//
	// RotMat=zeros(6,6);
	//
	// RotMat(1,1) =  cth0; RotMat(1,2) = -sth0;
	// RotMat(2,1) =  sth0; RotMat(2,2) =  cth0;
	// RotMat(3,3) =  1.0;
	// RotMat(4,4) =  cth0; RotMat(4,5) = -sth0;
	// RotMat(5,4) =  sth0; RotMat(5,5) =  cth0;
	// RotMat(6,6) =  1.0;
	//
	//
	// uxn=[0.0,0.0];
	// uzn=[0.0,0.0];
	// btn=[0.0,0.0];
	//
	// res = zeros(3,1);
	// B=zeros(3,6);
	// D=zeros(3,3);
	//
	// uxn(1) = soln(ndof*(IEN(e,3)-1)+1);
	// uzn(1) = soln(ndof*(IEN(e,3)-1)+2);
	// btn(1) = soln(ndof*(IEN(e,3)-1)+3);
	//
	// uxn(2) = soln(ndof*(IEN(e,4)-1)+1);
	// uzn(2) = soln(ndof*(IEN(e,4)-1)+2);
	// btn(2) = soln(ndof*(IEN(e,4)-1)+3);
	//
	// dummy = RotMat'*[uxn(1); uzn(1); btn(1); uxn(2); uzn(2); btn(2)];
	// uxn(1) = dummy(1);
	// uzn(1) = dummy(2);
	// btn(1) = dummy(3);
	// uxn(2) = dummy(4);
	// uzn(2) = dummy(5);
	// btn(2) = dummy(6);
	//
	//
	// nGP = 1;
	// [gpvec, gwvec] = get_Gauss_points(nGP);
	//
	// for gp = 1:nGP
	//     [N,dN_dx,d2N_dx2,J,xcoord]=shape_functions_Lagrange_1D(IEN(e,3:end), XX, p, gpvec(gp));
	//
	//     ux = 0.0;uz =0.0; bt = 0.0;
	//     dux = 0.0;duz = 0.0; dbt = 0.0;
	//
	//     for ii=1:nlocal
	//         ux  = ux  + uxn(ii) * N(ii);
	//         uz  = uz  + uzn(ii) * N(ii);
	//         bt  = bt  + btn(ii) * N(ii);
	//         dux = dux + uxn(ii) * dN_dx(ii);
	//         duz = duz + uzn(ii) * dN_dx(ii);
	//         dbt = dbt + btn(ii) * dN_dx(ii);
	//     end
	//
	//     sbt = sin(bt);
	//     cbt = cos(bt);
	//
	//     %compute average normal strain, shear strain and curvature
	//
	//     fact = (1.0+dux)*cbt - duz*sbt;
	//
	//     E = dux + 0.5*(dux*dux + duz*duz);
	//     G = (1.0+dux)*sbt + duz*cbt;
	//     K = dbt * fact;
	//
	//     % compute material response (elastic)
	//
	//     NF = EA * E; % normal force
	//     SF = GA * G; % shear force
	//     BM = EI * K; % bending moment
	//
	//     % multiply with volume element
	//
	//     dvol  = J*gwvec(gp);
	//     fact1 = dvol * af;
	//     NF    = NF * dvol;
	//     SF    = SF * dvol;
	//     BM    = BM * dvol;
	//     EAdv  = EA * fact1;
	//     GAdv  = GA * fact1;
	//     EIdv  = EI * fact1;
	//
	//     B(1,1) = (1.0+dux) * dN_dx(1);
	//     B(1,2) = duz * dN_dx(1);
	//     B(1,3) = 0.0;
	//
	//     B(2,1) = sbt * dN_dx(1);
	//     B(2,2) = cbt * dN_dx(1);
	//     B(2,3) = fact * N(1);
	//
	//     B(3,1) = dbt*cbt * dN_dx(1);
	//     B(3,2) = - dbt*sbt * dN_dx(1);
	//     B(3,3) = fact * dN_dx(1) - G*dbt * N(1);
	//
	//     B(1,4) = (1.0+dux) * dN_dx(2);
	//     B(1,5) = duz * dN_dx(2);
	//     B(1,6) = 0.0;
	//
	//     B(2,4) = sbt * dN_dx(2);
	//     B(2,5) = cbt * dN_dx(2);
	//     B(2,6) = fact * N(2);
	//
	//     B(3,4) = dbt*cbt * dN_dx(2);
	//     B(3,5) = - dbt*sbt * dN_dx(2);
	//     B(3,6) = fact * dN_dx(2) - G*dbt * N(2);
	//
	//     D(1,1) = EAdv;
	//     D(2,2) = GAdv;
	//     D(3,3) = EIdv;
	//
	//     res(1) = NF;
	//     res(2) = SF;
	//     res(3) = BM;
	//
	//     Klocal = Klocal + ( B'*D*B );
	//     Flocal = Flocal - ( B'*res );
	//
	//     fact1 = (+ SF * cbt - BM * dbt * sbt) * af;
	//     fact2 = (- SF * sbt - BM * dbt * cbt) * af;
	//
	//     for ii=1:nlocal
	//         TI   =  3*(ii-1)+1;
	//         TIp1 =  TI+1;
	//         TIp2 =  TI+2;
	//
	//         for jj=1:nlocal
	//             TJ   = 3*(jj-1)+1;
	//             TJp1 = TJ+1;
	//             TJp2 = TJ+2;
	//
	//             Klocal(TI,TJ)     =  Klocal(TI,TJ)     + dN_dx(ii)*NF*dN_dx(jj) * af;
	//             Klocal(TIp1,TJp1) =  Klocal(TIp1,TJp1) + dN_dx(ii)*NF*dN_dx(jj) * af;
	//
	//             fact3 =  dN_dx(ii)*BM*cbt*dN_dx(jj) * af;
	//             fact4 = -dN_dx(ii)*BM*sbt*dN_dx(jj) * af;
	//
	//             Klocal(TI,TJp2)   = Klocal(TI,TJp2)    + (fact3 + dN_dx(ii)*fact1*N(jj) );
	//             Klocal(TIp1,TJp2) = Klocal(TIp1,TJp2)  + (fact4 + dN_dx(ii)*fact2*N(jj) );
	//             Klocal(TIp2,TJ)   = Klocal(TIp2,TJ)    + (fact3 + N(ii)*fact1*dN_dx(jj) );
	//             Klocal(TIp2,TJp1) = Klocal(TIp2,TJp1)  + (fact4 + N(ii)*fact2*dN_dx(jj) );
	//             Klocal(TIp2,TJp2) = Klocal(TIp2,TJp2)  + (N(ii)*(-SF*G-BM*dbt*fact)*N(jj) - dN_dx(ii)*BM*G*N(jj) - N(ii)*BM*G*dN_dx(jj)) * af;
	//         end
	//     end
	// end
	//
	// h = XX(IEN(e,4)) - XX(IEN(e,3));
	//
	// Flocal(1) = Flocal(1) + 0.5*h*bf(1);
	// Flocal(4) = Flocal(4) + 0.5*h*bf(1);
	//
	// Flocal(2) = Flocal(2) + 0.5*h*bf(2);
	// Flocal(5) = Flocal(5) + 0.5*h*bf(2);
	//
	// Flocal = RotMat*Flocal;
	// Klocal = RotMat*Klocal*RotMat';
	//
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

func processfile() {
	// function [ndim, ndof, nnode, nelem, coords, elemConn, elemData, LM, neq, assy4r, dof_force, Fext, maxloadSteps, loadincr, outputlist] = processfile(fname)
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
	//
	// % ndof
	//
	// line=fgets(fid);
	// linestr = strsplit(line, ",");
	// ndof = int32(str2num(linestr{1,2}))
	//
	// % nodes
	//
	// line=fgets(fid);
	// linestr = strsplit(line, ",");
	// nnode = int32(str2num(linestr{1,2}))
	//
	// nperelem = 2;
	// nsize = nperelem*ndof;
	// neq = nnode*ndof;
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
	//
	//
	// % element data
	//
	// line=fgets(fid);
	// linestr = strsplit(line, ",");
	// nelemData = int32(str2num(linestr{1,2}))
	// elemData = zeros(nelemData, 10);
	// for i=1:nelemData
	//     line = fgets(fid);
	//     linestr = strsplit(line, ",");
	//     for j=1:10
	//       elemData(i,j) = double(str2num(linestr{1,j+1}));
	//     end
	// end
	//
	// % elements
	//
	// line=fgets(fid);
	// linestr = strsplit(line, ",");
	// nelem = int32(str2num(linestr{1,2}))
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
	//
	// % Dirichlet boundary conditions
	//
	// line=fgets(fid);
	// linestr = strsplit(line, ",");
	// nDBC    = int32(str2num(linestr{1,2}))
	//
	// %dbclist = zeros(nDBC, 3);
	// dbcnodes = zeros(nDBC, 1, "int32");
	// for i=1:nDBC
	//     line = fgets(fid);
	//     linestr = strsplit(line, ",");
	//     n1 = int32(str2num(linestr{1,1}));
	//     n2 = int32(str2num(linestr{1,2}));
	// %    dbclist(i,3) = double(str2num(linestr{1,3}));
	//     dbcnodes(i) = (n1-1)*ndof+n2;
	// end
	//
	// assy4r = setdiff([1:neq], dbcnodes)';
	//
	// % Force boundary conditions
	//
	// line=fgets(fid);
	// linestr = strsplit(line, ",");
	// nFBC    = int32(str2num(linestr{1,2}))
	// fbclist = zeros(nFBC, 3);
	// dof_force = zeros(nFBC, 1, "int32");
	// Fext = zeros(neq, 1);
	// for i=1:nFBC
	//     line = fgets(fid);
	//     linestr = strsplit(line, ",");
	//     n1 = int32(str2num(linestr{1,1}));
	//     n2 = int32(str2num(linestr{1,2}));
	//     ind = (n1-1)*ndof + n2;
	//     dof_force(i) = ind;
	//     Fext(ind) = double(str2num(linestr{1,3}));
	// end
	//
	// % for output
	//
	// line=fgets(fid);
	// linestr = strsplit(line, ",");
	// nOutput = int32(str2num(linestr{1,2}))
	// outputlist = zeros(nOutput, 1, "int32");
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
	//
	// line=fgets(fid)
	// linestr = strsplit(line, ",");
	// maxloadSteps = int32(str2num(linestr{1,1}));
	//
	// line=fgets(fid)
	// linestr = strsplit(line, ",");
	// loadincr = double(str2num(linestr{1,1}));
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
}

func Shape_functions_Lagrange_1D(NodeNums []int, XX [][]float64, p int, xi float64) (
	N, dN_dx, d2N_dx2 []float64, J, x float64,
) {
	// function [N,dN_dx,d2N_dx2,J,x]=shape_functions_Lagrange_1D(NodeNums, XX, p, xi)
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
		Jx = Jx + dN_dxi[kk]*XX[NodeNums[kk]][0] //1]
		Jy = Jy + dN_dxi[kk]*XX[NodeNums[kk]][1] // 2]
		x = x + N[kk]*XX[NodeNums[kk]][0]        // TODO not clear?
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

func solve_arclength_split() {
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
	// function [Klocal, Flocal]=Truss_2D_model1(elmDat, IEN, e, XX, soln, bf)
	//
	// ndof = 2;
	//
	// finite = (int32(elmDat(IEN(e,1), 1)) == 1) ;
	// rho0 = elmDat(IEN(e,1),2);
	// A0   = elmDat(IEN(e,1),3);
	// E    = elmDat(IEN(e,1),4);
	//
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//
	// disp  = zeros(4,1);
	//
	// X1 = XX(IEN(e,3),1);
	// Y1 = XX(IEN(e,3),2);
	// X2 = XX(IEN(e,4),1);
	// Y2 = XX(IEN(e,4),2);
	//
	// disp(1) = soln(ndof*(IEN(e,3)-1)+1);
	// disp(2) = soln(ndof*(IEN(e,3)-1)+2);
	// disp(3) = soln(ndof*(IEN(e,4)-1)+1);
	// disp(4) = soln(ndof*(IEN(e,4)-1)+2);
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
// function [Klocal, Flocal]=Truss_2D_model2(elmDat, IEN, e, XX, soln, bf)
//
// ndof = 2;
//
// finite = (int32(elmDat(IEN(e,1), 1)) == 1) ;
// rho0 = elmDat(IEN(e,1),2);
// A0   = elmDat(IEN(e,1),3);
// E    = elmDat(IEN(e,1),4);
//
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// disp  = zeros(4,1);
//
// X1 = XX(IEN(e,3),1);
// Y1 = XX(IEN(e,3),2);
// X2 = XX(IEN(e,4),1);
// Y2 = XX(IEN(e,4),2);
//
// disp(1) = soln(ndof*(IEN(e,3)-1)+1);
// disp(2) = soln(ndof*(IEN(e,3)-1)+2);
// disp(3) = soln(ndof*(IEN(e,4)-1)+1);
// disp(4) = soln(ndof*(IEN(e,4)-1)+2);
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
// function [Klocal, Flocal]=Truss_3D_model2(elmDat, IEN, e, XX, soln, bf)
//
// ndof = 3;
//
// finite = (int32(elmDat(1)) == 1) ;
// rho0 = elmDat(2);
// A0   = elmDat(3);
// E    = elmDat(4);
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// disp  = zeros(6,1);
//
// % rotate nodal displacements and compute nodal positions on element axis
//
// X1 = XX(IEN(e,3),1);
// Y1 = XX(IEN(e,3),2);
// Z1 = XX(IEN(e,3),3);
// X2 = XX(IEN(e,4),1);
// Y2 = XX(IEN(e,4),2);
// Z2 = XX(IEN(e,4),3);
//
// disp(1) = soln(ndof*(IEN(e,3)-1)+1);
// disp(2) = soln(ndof*(IEN(e,3)-1)+2);
// disp(3) = soln(ndof*(IEN(e,3)-1)+3);
// disp(4) = soln(ndof*(IEN(e,4)-1)+1);
// disp(5) = soln(ndof*(IEN(e,4)-1)+2);
// disp(6) = soln(ndof*(IEN(e,4)-1)+3);
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
