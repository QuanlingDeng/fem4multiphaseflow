#include "../linalg/linalg_header.h"
#include "../fem/fem_header.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include "../linalg/linalg_header.h"
#include "../mesh/mesh_header.h"

using namespace std;


/*
This is a test file for the quadrilateral/rectangular mesh
 */

//======================================================

double evalk(Array<double> &coord)
{ 
  return 1;  
}

//======================================================

double evalf(Array<double> &coord)
{
  return 0;
}

//======================================================

double evaldir(Array<double> &coord)
{
  return 1.0 - coord[0];
}

//======================================================

int main()
{

  Array<int> N(2);//Array of number of elements
  N[0] = 4;
  N[1] = 4;

  Array<double> L(4);
  L[0] = 0.0;
  L[1] = 1.0;
  L[2] = 0.0;
  L[3] = 1.0;

  double dx = 1.0/N[0];
  double dy = 1.0/N[1];

  Mesh *mesh;
  mesh = new TMesh<2>(N, L, Element::QUADRILATERAL);
  int nnn = (N[0]+1)*(N[1]+1);

  ofstream out;
  out.open("mesh.out");
  mesh->Print(out, Element::QUADRILATERAL);
  out.close();
  
  mesh->ParaviewPrint(Element::QUADRILATERAL);
   
  Array<double *> linellipticdata(1);
  Array<double *> linforcedata(1);

  linellipticdata[0] = new double[nnn];
  linforcedata[0] = new double[nnn];
  
  int index=0;
  Array<double> coord(2);

  coord[1] = 0.0;
  for(int j=0; j<=N[1]; j++)
    {
      coord[0] = 0.0;
      for(int i=0; i<=N[0]; i++)
	{
	  linellipticdata[0][index] = evalk(coord);
	  linforcedata[0][index] = evalf(coord);
	  coord[0] += dx;
	  index++;
	}
      coord[1] += dy;
    }
 
  Data *data;
  data = new Data();
  data->SetEllipticData(linellipticdata);
  data->SetForceData(linforcedata);


  FEM *fem = new BilinearFEM2D(mesh, data);
  fem->Assemble();

  Vector Dval;
  fem->ProjectAFunction(evaldir, Dval);
  
  int nbdry = mesh->GetNBdrs();
  Array<bool *> bdry_is_dirichlet(1);
  bdry_is_dirichlet[0] = new bool[nbdry];
  for(int i=0; i<4; i++)
    bdry_is_dirichlet[0][i] = false;  

  bdry_is_dirichlet[0][1] = true; 
  bdry_is_dirichlet[0][3] = true;  



  int NV = mesh->GetNV();
  Array<double *> dval(1);
  dval[0] = new double[NV];
  for(int i=0; i<NV; i++)
    {
      dval[0][i] = Dval(i);
    }

  fem->SetDirichletBoundary(bdry_is_dirichlet, dval);

  Vector sol;
  sol = Dval;
  fem->Solve(sol, "multifrontal", 10000, 1, 1.0e-50, 1.0e-80);

  out.open("sol.out");
  mesh->Print(out, sol);
  out.close();
  
  mesh->ParaviewPrint(sol, Element::QUADRILATERAL);

  delete mesh;
  delete []linellipticdata[0];
  delete []linforcedata[0];
  delete []bdry_is_dirichlet[0];
  delete []dval[0];
  delete data;
  delete fem;

  
}

