#include "../linalg/linalg_header.h"
#include "../fem/fem_header.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include "../linalg/linalg_header.h"
#include "../mesh/mesh_header.h"

using namespace std;

//======================================================

double evalk(double x)
{ 

  return 1.0;
  //  double pi = 4.0*atan(1.0);
  //return 1.0/(1.0-0.9*sin(6.0*pi*x)); 

  //return 1.0/(1+x*x); 
}

//======================================================

double evalf(double x)
{
  return 1.0;
}

//======================================================

double exactsol(double x)
{
  
  // double pi = 4.0*atan(1.0);
  //if(x<0.5)
  //return x+x*x*x/3;
  //else
  //return 13./24.;
  // return x - x*x/2 + 3.0/(20.0*pi)*(cos(6.0*pi*x)-1.0)+1.0/(40*pi*pi)*(sin(6.0*pi*x)-6*pi*x*cos(6*pi*x));

  return -x*x*0.5+x*0.5;
}

//======================================================

int main()
{
  int Nx;
  cout << "input Nx : " << endl;
  cin >> Nx;
  double dx = 1.0/Nx;

  Mesh *mesh;
  mesh = new TMesh<1>(Nx, 0.0, 1.0);
  
  int nnx = Nx+1;
  double *grid = new double[Nx];
  for(int i=0; i<Nx; i++)
    grid[i] = dx;
  
  Array<double *> linellipticdata(1);
  Array<double *> linforcedata(1);
  for(int i=0; i<2; i++)
    {
      linellipticdata[i] = new double[nnx];
      linforcedata[i] = new double[nnx];
    }
  
  double coord = 0.0;
  for(int i=0; i<Nx; i++)
    {
      linellipticdata[0][i] = evalk(coord);
      linforcedata[0][i] = evalf(coord);
      coord += grid[i];
    }
  
  linellipticdata[0][Nx] = evalk(coord);
  linforcedata[0][Nx] = evalf(coord);
  
  coord = 0.0;
  double h; 
  for(int i=0; i<=Nx; i++)
    {
      h = grid[i];      
      linellipticdata[0][i] = evalk(coord);
      linforcedata[0][i] = evalf(coord);      
      coord += h;
    }
  
  
  Data *data;
  data = new Data();
  data->SetEllipticData(linellipticdata);
  data->SetForceData(linforcedata);
  
  LinearFEM1D *linearfem = new LinearFEM1D(mesh, data);  
  linearfem->Assemble();
  Array<bool> bdry_is_dirichlet(2);
  for (int i=0; i<2; i++)
    bdry_is_dirichlet[i] = true;


  Vector Dval(Nx+1);
  Dval = 0.0;
  Dval(0) = 1.0;

  linearfem->SetDirichletBoundary(bdry_is_dirichlet, Dval);

  Vector cubicsol;
  linearfem->Solve(cubicsol);
  
  ofstream exact("exact.out");
  ofstream femsol("femsol.out");

  coord = 0.0;
  for(int i=0; i<=Nx; i++)
    {
      exact << coord << " " << exactsol(coord) << endl;
      femsol << coord << " " << cubicsol(i) << endl; 
      coord += grid[i];
    }
  
  exact.close();
  femsol.close();
  
  delete mesh;
}

