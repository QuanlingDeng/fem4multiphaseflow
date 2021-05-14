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
  return x-0.5;
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

  return -x*x*x/6+x*x/4-1/24.0;
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
  
  int nnx = 3*Nx+1;
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
  for(int i=0; i<Nx; i++)
    {
      h = grid[i]/3.0;
      coord += h;
      int k = Nx+1+2*i;
      
      linellipticdata[0][k] = evalk(coord);
      linforcedata[0][k] = evalf(coord);
      
      coord += h;
      k++;
      
      linellipticdata[0][k] = evalk(coord);
      linforcedata[0][k] = evalf(coord);
      
      coord += h;
    }
  
  
  Data *data;
  data = new Data();
  data->SetEllipticData(linellipticdata);
  data->SetForceData(linforcedata);
  
  CubicFEM1D *cubicfem = new CubicFEM1D(mesh, data);  
  cubicfem->Assemble();

  Vector cubicsol;
  cubicfem->Solve(cubicsol);
  
  ofstream exact("exact.out");
  ofstream femsol("femsol.out");

  coord = 0.0;
  for(int i=0; i<Nx; i++)
    {
      exact << coord << " " << exactsol(coord) << endl;
      femsol << coord << " " << cubicsol(i) << endl; 
      coord += grid[i];
    }
  
  exact << coord << " " << exactsol(coord) << endl;
  femsol << coord << " " << cubicsol(Nx) << endl;
  
  coord = 0.0;
  for(int i=0; i<Nx; i++)
    {
      double h = grid[i]/3.0;
      coord += h;
      int k = Nx+1+2*i;
      
      exact << coord << " " << exactsol(coord) << endl;
      femsol << coord << " " << cubicsol(k) << endl;
      
      coord += h;
      k++;
      
      exact << coord << " " << exactsol(coord) << endl;
      femsol << coord << " " << cubicsol(k) << endl;

      coord += h;
    }
  
  exact.close();
  femsol.close();
  

  delete mesh;


}

