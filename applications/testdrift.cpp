#include "../fem/fem_header.h"
#include "../mesh/mesh_header.h"
#include "../general/array.hpp"
#include "../timedependent/time_header.h"
#include <cmath>

double sg[2] = {-sqrt(3)/3,sqrt(3)/3};
double k(Array<double> &xy) { return 1.0; }
double kn(Array<double> &xy) { return 0.001; }
double kp(Array<double> &xy) { return 0.001; }
double u(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return x*x*x - y*y*y + 1.0; }
double un(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return 6.0*x + x*y; }
double up(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return 6.0*y + x*y; }
double f(double n, double p) { return p - n; }
double fn(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return 54*x*x + 9*x*x*y - 9*x*y*y - 36*x*y; }
double fp(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return 54*y*y + 9*x*y*y - 9*x*x*y - 36*x*y; }
double ub(Array<double> &xy) { double x=xy[0]; double y=xy[1]; 
  if( fabs(x) < 1.0e-7 )
    return 1 - y*y*y;
  else if( fabs(y) < 1.0e-7 )
    return 1 + x*x*x;
  else if( fabs(x-1.0) < 1.0e-7 )
    return 2.0 - y*y*y;
  else if( fabs(y-1.0) < 1.0e-7 )
    return x*x*x;
  else
    return 0.0;
}
double unb(Array<double> &xy) { double x=xy[0]; double y=xy[1];
  if( fabs(x) < 1.0e-7 )
    return 0.0;
  else if( fabs(y) < 1.0e-7 )
    return 6.0*x;
  else if( fabs(x-1.0) < 1.0e-7 )
    return 6.0 + y;
  else if( fabs(y-1.0) < 1.0e-7 )
    return 7.0*x;
  else
  return 0.0;
}
double upb(Array<double> &xy) { double x=xy[0]; double y=xy[1];
  if( fabs(x) < 1.0e-7 )
    return 6.0*y;
  else if( fabs(y) < 1.0e-7 )
    return 0.0;
  else if( fabs(x-1.0) < 1.0e-7 )
    return 7.0*y;
  else if( fabs(y-1.0) < 1.0e-7 )
    return 6.0 + x;
  else
  return 0.0;
}


int main(int argc, const char * argv[])
{
  int Nx = 40;
  int Ny = 40;
  double hx = 1.0/Nx;
  double hy = 1.0/Ny;

  Array<int> n(2);
  n[0] = Nx;
  n[1] = Ny;

  Array<double> L(4);
  L[0] = 0.0;
  L[1] = 1.0;
  L[2] = 0.0;
  L[3] = 1.0;

  Mesh *tmesh;
  tmesh = new TMesh<2> (n, L, Element::TRIANGLE);

  int NV = tmesh->GetNV();
  int NE = tmesh->GetNE();

  Array<double *> edata(1);
  edata[0] = new double[NV];

  Array<double *> adata(1);
  adata[0] = new double[NV];

  Array<double *> fdata(1);
  fdata[0] = new double[NV];

  Array<double> coord(2);

  Data *tdata = new Data();
  Data *data = new Data();

  Vector Dval(NV);
  Vector uu(NV);
  Vector sol(NV);
  Vector soln(NV);
  Vector solp(NV);
  Vector ssol(NV);
  Vector ssoln(NV);
  Vector ssolp(NV);
  Dval = 0.0;
  sol = 0.0;
  soln = 0.0;
  solp = 0.0;
  ssol = 0.0;
  ssoln = 0.0;
  ssolp = 0.0;

  Array<double *> dval(1);
  dval[0] = new double[NV];
  
  Array<bool *> bdry_is_dirichlet(1);
  bdry_is_dirichlet[0] = new bool[4];
  for(int i=0; i<4; i++)
    bdry_is_dirichlet[0][i] = true;

  for(int i=0; i<NV; i++)
    {
      coord[0] = tmesh->GetVertex(i, 0);
      coord[1] = tmesh->GetVertex(i, 1);
      edata[0][i] = k(coord);
      fdata[0][i] = f( unb(coord), upb(coord) );
      uu(i) = u(coord);
    }
  
  LinearFEM2D *UW;
  LinearFUPG2D *UWn;
  LinearFUPG2Dp *UWp;
  
  // Iterative Solution
  int NN = 200;
  double tol = 0.0;
  for(int I=1; I<NN; I++)
    {
      ssol = sol;
      ssoln = soln;
      ssolp = solp;
      data->SetEllipticData(edata);
      data->SetForceData(fdata);
      
      UW = new LinearFEM2D(tmesh, data);
      UW->Assemble();
      UW->ProjectAFunction(ub, Dval);
      
      for(int i=0; i<NV; i++)
	dval[0][i] = Dval(i);
      
      UW->SetDirichletBoundary(bdry_is_dirichlet, dval);
      
      //sol = Dval;
      UW->Solve(sol, "pcg", 1000, 0, 1.0e-16, 1.0e-32);

      // ======= solve the second equation =================
      for(int i=0; i<NV; i++)
	{
	  coord[0] = tmesh->GetVertex(i, 0);
	  coord[1] = tmesh->GetVertex(i, 1);
	  edata[0][i] = kn(coord);
	  adata[0][i] = sol(i);
	  fdata[0][i] = fn(coord);
	}

      tdata->SetEllipticData(edata);
      tdata->SetAdvectionData(adata);
      tdata->SetForceData(fdata);

      UWn = new LinearFUPG2D(tmesh, tdata);
      UWn->Assemble();      
      UWn->ProjectAFunction(unb, Dval);

      for(int i=0; i<NV; i++)
	dval[0][i] = Dval(i);
      
      UWn->SetDirichletBoundary(bdry_is_dirichlet, dval);
      
      //soln = Dval;
      UWn->Solve(soln, "gmres", 100, 0, 1.0e-16, 1.0e-32);
            
      // ======== solve the third equation ===================
      for(int i=0; i<NV; i++)
	{
	  coord[0] = tmesh->GetVertex(i, 0);
	  coord[1] = tmesh->GetVertex(i, 1);
	  edata[0][i] = kp(coord);
	  adata[0][i] = sol(i);
	  fdata[0][i] = fp(coord);
	}

      tdata->SetEllipticData(edata);
      tdata->SetAdvectionData(adata);
      tdata->SetForceData(fdata);
      
      UWp = new LinearFUPG2Dp(tmesh, tdata);
      UWp->Assemble();
      UWp->ProjectAFunction(upb, Dval);
      
      for(int i=0; i<NV; i++)
	dval[0][i] = Dval(i);
      
      UWp->SetDirichletBoundary(bdry_is_dirichlet, dval);
      
      //solp = Dval;
      UWp->Solve(solp, "gmres", 100, 0, 1.0e-16, 1.0e-32);
      
      for(int i=0; i<NV; i++)
	{
	  coord[0] = tmesh->GetVertex(i, 0);
	  coord[1] = tmesh->GetVertex(i, 1);
	  edata[0][i] = k(coord);
	  fdata[0][i] = f( soln(i), solp(i) );
	}
      
      tol = 0.0;
      for(int i=0; i<NV; i++)
	{
	  if( tol<fabs(sol(i)-ssol(i)) )
	    tol = fabs(sol(i)-ssol(i));
	  if( tol<fabs(soln(i)-ssoln(i)) )
	    tol = fabs(soln(i)-ssoln(i));
	  if( tol<fabs(solp(i)-ssolp(i)) )
	    tol = fabs(solp(i)-ssolp(i));
	}
      cout<<I<<"th tol: "<<tol<<endl;
      //sol.Print();
      if(tol<0.45)
	{
	  cout<<"Good: Reached the solution before end of the iterations."<<endl;
	  break;
	}
      else
	{
	  delete UWp;
	  delete UWn;
	  delete UW;
	} 
      
    }
  
  tmesh->ParaviewPrint(Element::TRIANGLE);
  tmesh->ParaviewPrint(sol, Element::TRIANGLE);  

 
  //============== L2 Error ==============================================
  double err, uh, tr = 0.0;
  double xt[2] = {0.5 - 0.5/sqrt(3), 0.5 + 0.5/sqrt(3)};
  double yt[2] = {0.5 - 0.5/sqrt(3), 0.5 + 0.5/sqrt(3)};
  double xtr[2] = {-1.0/sqrt(3), 1.0/sqrt(3)};
  double ytr[2] = {-1.0/sqrt(3), 1.0/sqrt(3)};
  int indd[4];

  err = 0.0;
  for(int j=0; j<Ny; j++)
    for(int i=0; i<Nx; i++)
      {
	indd[0] = i + j*(Nx+1);
	indd[1] = indd[0] + 1;
	indd[2] = indd[1] + Nx + 1;
	indd[3] = indd[2] - 1;
	for(int jj=0; jj<2; jj++)
	  for(int ii=0; ii<2; ii++)
	    {
	      uh = sol(indd[0])*(1 - xt[ii])*(1 - yt[jj]) + sol(indd[1])*xt[ii]*(1 - yt[jj]) + sol(indd[2])*xt[ii]*yt[jj] + sol(indd[3])*(1 - xt[ii])*yt[jj];
	      
	      coord[0] = (i+0.5) * hx + 0.5*hx*xtr[ii];
	      coord[1] = (j+0.5) * hy + 0.5*hy*ytr[jj];
	      
	      tr = u(coord);
	      err += (uh - tr)*(uh - tr);	    
	    }
      } 
  cout<<"The error for u is: "<<sqrt(0.25*hx*hy*err)<<endl;
  
  err = 0.0;
  for(int j=0; j<Ny; j++)
    for(int i=0; i<Nx; i++)
      {
	indd[0] = i + j*(Nx+1);
	indd[1] = indd[0] + 1;
	indd[2] = indd[1] + Nx + 1;
	indd[3] = indd[2] - 1;
	for(int jj=0; jj<2; jj++)
	  for(int ii=0; ii<2; ii++)
	    {
	      uh = soln(indd[0])*(1 - xt[ii])*(1 - yt[jj]) + soln(indd[1])*xt[ii]*(1 - yt[jj]) + soln(indd[2])*xt[ii]*yt[jj] + soln(indd[3])*(1 - xt[ii])*yt[jj];
	      
	      coord[0] = (i+0.5) * hx + 0.5*hx*xtr[ii];
	      coord[1] = (j+0.5) * hy + 0.5*hy*ytr[jj];
	      
	      tr = un(coord);
	      err += (uh - tr)*(uh - tr);	    
	    }
      } 
  cout<<"The error for un is: "<<sqrt(0.25*hx*hy*err)<<endl;
  
  err = 0.0;
  for(int j=0; j<Ny; j++)
    for(int i=0; i<Nx; i++)
      {
	indd[0] = i + j*(Nx+1);
	indd[1] = indd[0] + 1;
	indd[2] = indd[1] + Nx + 1;
	indd[3] = indd[2] - 1;
	for(int jj=0; jj<2; jj++)
	  for(int ii=0; ii<2; ii++)
	    {
	      uh = solp(indd[0])*(1 - xt[ii])*(1 - yt[jj]) + solp(indd[1])*xt[ii]*(1 - yt[jj]) + solp(indd[2])*xt[ii]*yt[jj] + solp(indd[3])*(1 - xt[ii])*yt[jj];
	      
	      coord[0] = (i+0.5) * hx + 0.5*hx*xtr[ii];
	      coord[1] = (j+0.5) * hy + 0.5*hy*ytr[jj];
	      
	      tr = up(coord);
	      err += (uh - tr)*(uh - tr);	    
	    }
      } 
  cout<<"The error for up is: "<<sqrt(0.25*hx*hy*err)<<endl;
  

   // ============= compute flux =========================================
  Array<double *> flux(3);
  for(int i=0; i<3; i++) flux[i] = new double[NE];
  
  Array<double> fluxcv(NV);
  Array<double> fluxncv(NV);
  UWn->ComputeConservativeFlux(sol, soln, flux, fluxcv, fluxncv);

  ofstream fileout("n.out");
  for(int i=0; i<NV; i++) 
    fileout << i << " \t" << fluxcv[i] << "\t" << fluxncv[i] << endl;
  fileout.close();

  UWp->ComputeConservativeFlux(sol, solp, flux, fluxcv, fluxncv);
  ofstream fileout1("p.out");
  for(int i=0; i<NV; i++) 
    fileout1 << i << "\t" << fluxcv[i] << "\t" << fluxncv[i] << endl;
  fileout1.close();


  delete UWp;
  delete UWn;
  delete UW;

  delete []bdry_is_dirichlet[0];
  delete []dval[0];
  for(int i=0; i<3; i++)
    delete []flux[i];

  delete []adata[0];
  delete []edata[0];
  delete []fdata[0];
  delete data;
  delete tdata;
  delete tmesh;

  return 0;
}
