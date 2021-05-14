#include "../fem/fem_header.h"
#include "../mesh/mesh_header.h"
#include "../general/array.hpp"
#include "../timedependent/time_header.h"
#include <fstream>
#include <cmath>

double sg[2] = {-sqrt(3)/3,sqrt(3)/3};
/*
// ============== case 1 ========================
double k(Array<double> &xy) { double x=xy[0]; return 1.0 + 2.0*x; }
double vx(Array<double> &xy) { double x=xy[0]; return x + x*x; }
double vy(Array<double> &xy) { double y=xy[1]; return y + y*y; }
double f(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return 2*x*(5-2*x*x)*(y-y*y) + 2*(x-x*x)*(1+2*x+y-2*y*y*y); }
double u(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return (x-x*x)*(y-y*y); }
double unb(Array<double> &xy) { return 0.0; }
double uux(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return (1-2*x)*(y-y*y); }
double uuy(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return (x-x*x)*(1-2*y); }


// ============== case 2 ========================
double k(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return exp( -(x-x*x)*(y-y*y) ); }
double vx(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return exp( -(x-x*x)*(y-y*y) ); }
double vy(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return exp( -(x-x*x)*(y-y*y) ); }
double u(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return exp( (x-x*x)*(y-y*y) ) - 1.0; }
double f(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return 2*(x-x*x+y-y*y) + (1-2*x)*(y-y*y)*exp( -(x-x*x)*(y-y*y) ) + (x-x*x)*(1-2*y)*exp( -(x-x*x)*(y-y*y) ); }
double unb(Array<double> &xy) { return 0.0; }
double uux(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return (1-2*x)*(y-y*y)*exp( (x-x*x)*(y-y*y) ); }
double uuy(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return (x-x*x)*(1-2*y)*exp( (x-x*x)*(y-y*y) ); }
*/

// ============== case 3 ========================
double pi = 4.0*atan(1.0);
double k(Array<double> &xy) { double y=xy[1]; return y / (pi*pi); }
double vx(Array<double> &xy) { double y=xy[1]; return y / pi; }
double vy(Array<double> &xy) { double x=xy[0]; return x / pi; }
double u(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return sin(pi*x)*sin(pi*y); }
double f(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return 2.0*y*sin(pi*x)*sin(pi*y) - sin(pi*x)*cos(pi*y) / pi + y*cos(pi*x)*sin(pi*y) + x*sin(pi*x)*cos(pi*y); }
double unb(Array<double> &xy) { return 0.0; }
double uux(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return pi*cos(pi*x)*sin(pi*y); }
double uuy(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return pi*sin(pi*x)*cos(pi*y); }

/*
// ============== case 4 ========================
double pi = 4.0*atan(1.0);
double k(Array<double> &xy) { return 1.0; }
double vx(Array<double> &xy) { return 1.0; }
double vy(Array<double> &xy) { return 1.0; }
double u(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return sin(pi*x)*(1.0 - cos(2.0*pi*y)); }
double uux(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return pi*cos(pi*x)*(1.0 - cos(2.0*pi*y)); }
double uuy(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return 2.0*pi*sin(pi*x)*sin(2.0*pi*y); }
double f(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return pi*pi*sin(pi*x)*(1.0 - 5.0*cos(2.0*pi*y)) + vx(xy)*uux(xy) + vy(xy)*uuy(xy); }
double unb(Array<double> &xy) { return 0.0; }
*/


int main(int argc, const char * argv[])
{
  int Nx = 100;
  int Ny = 100;
  double hx= 1.0/Nx;
  double hy = 1.0/Ny;

  Array<int> n(2);
  n[0] = Nx;
  n[1] = Ny;

  Array<double> L(4);
  L[0] = 0.0;
  L[1] = 1.0;
  L[2] = 0.0;
  L[3] = 1.0;

  cout<< "hello " << exp(1000*0.999-1000) - exp(-1000)  <<endl;

  Mesh *mesh;
  mesh = new TMesh<2> (n, L, Element::TRIANGLE);

  int NV = mesh->GetNV();

  Array<double *> edata(1);
  edata[0] = new double[NV];

  Array<double *> adata(2);
  adata[0] = new double[NV];
  adata[1] = new double[NV];

  Array<double *> fdata(1);
  fdata[0] = new double[NV];
  
  Array<double> coord(2);
  for(int i=0; i<NV; i++)
    {
      coord[0] = mesh->GetVertex(i, 0);
      coord[1] = mesh->GetVertex(i, 1);
      edata[0][i] = k(coord);
      adata[0][i] = vx(coord);
      adata[1][i] = vy(coord);
      fdata[0][i] = f(coord);
    }

  Data *data;
  data = new Data();
  data->SetEllipticData(edata);
  data->SetConservativeAdvectionData(adata);
  data->SetForceData(fdata);

  LinearFEM2D *UW = new LinearFEM2D(mesh, data);
  UW->Assemble();

  Vector Dval(NV);
  Dval = 0.0;
  UW->ProjectAFunction(unb, Dval);

  Array<double *> dval(1);
  dval[0] = new double[NV];
    for(int i=0; i<NV; i++)
      dval[0][i] = Dval(i);

  Array<bool *> bdry_is_dirichlet(1);
  bdry_is_dirichlet[0] = new bool[4];
  for(int i=0; i<4; i++)
    bdry_is_dirichlet[0][i] = true;

  UW->SetDirichletBoundary(bdry_is_dirichlet, dval);

  Vector sol(NV);
  sol = Dval;
  char solv[20] = "multifrontal";
  UW->Solve(sol, solv, 1000, 1, 1.0e-12, 1.0e-24);
 
  // ============= compute flux =========================================
  Array<double *> flux(3);
  for(int i=0; i<3; i++) flux[i] = new double[mesh->GetNE()];

  //UW->ComputeConservativeFlux(sol, flux);
  UW->ComputeAdvectionDiffusionConservativeFlux(sol, flux);

  mesh->ParaviewPrint(Element::TRIANGLE);
  mesh->ParaviewPrint(sol, Element::TRIANGLE);  

  //============== L2 Error ==============================================
  double err, uh, tr = 0.0;
  double xt[2] = {0.5 - 0.5/sqrt(3), 0.5 + 0.5/sqrt(3)};
  double yt[2] = {0.5 - 0.5/sqrt(3), 0.5 + 0.5/sqrt(3)};
  double xtr[2] = {-1.0/sqrt(3), 1.0/sqrt(3)};
  double ytr[2] = {-1.0/sqrt(3), 1.0/sqrt(3)};
  int indd[4];

  err = 0.0;
  double err2 = 0.0;
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

	      uh = ( sol(indd[0])*(yt[jj] - 1) + sol(indd[1])*(1 - yt[jj]) + sol(indd[2])*yt[jj] - sol(indd[3])*yt[jj] ) / hx;
	      uh -= uux(coord);
	      uh = uh*uh;
	      tr = ( sol(indd[0])*(xt[ii] - 1) - sol(indd[1])*xt[ii] + sol(indd[2])*xt[ii] + sol(indd[3])*(1 - xt[ii]) ) / hy;
	      tr -= uuy(coord);
	      tr = tr*tr;
	      err2 += uh + tr;	    
	    }
      } 
  cout<<"The L2 error is: "<<sqrt(0.25*hx*hy*err)<<endl;
  cout<<"The semi-norm error is: "<<sqrt(0.25*hx*hy*err2)<<endl;
  
  // ================ Output Data =========================================
  int tt = 0;
  ofstream fileout("2d.z");
  fileout<<"! nx "<<Nx+1<<" ny "<<Ny+1<<" xmin 0 xmax 1 ymin 0 ymax 1"<<endl;
  for(int j=0; j<=Ny; j++)
    {
      for(int i=0; i<=Nx; i++)
	{
	  tt = i + j*(Nx + 1);
	  fileout<<setprecision(4)<<sol(tt)<<" ";
	}
      fileout<<endl;
    }  
  fileout.close(); 
  

  delete []bdry_is_dirichlet[0];
  delete []dval[0];
  for(int i=0; i<3; i++)
    delete []flux[i];
  delete []edata[0];
  delete []adata[1];
  delete []adata[0];
  delete []fdata[0];
  delete data;
  delete mesh;
  delete UW;

  return 0;
}
