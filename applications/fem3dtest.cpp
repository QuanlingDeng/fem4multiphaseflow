#include "../linalg/linalg_header.h"
#include "../fem/fem_header.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include "../linalg/linalg_header.h"
#include "../mesh/mesh_header.h"

using namespace std;

//======================================================

double evalk(Array<double> &coord)
{ 
  return 1;  
}

//======================================================

double evalf(Array<double> &coord)
{
  return 2;
}

//======================================================

double evaldir(Array<double> &coord)
{
  return 1;
}

//======================================================

int main()
{

  Array<int> N(3);
  N[0] = 2;
  N[1] = 1;
  N[2] = 1;

  Array<double> L(6);
  L[0] = 0.0;
  L[1] = 1.0;
  L[2] = 0.0;
  L[3] = 1.0;
  L[4] = 0.0;
  L[5] = 1.0;

  double dx = 1.0/N[0];
  double dy = 1.0/N[1];
  double dz = 1.0/N[2];

  Mesh *mesh;
  mesh = new TMesh<3>(N, L, Element::TETRAHEDRAL);
  int nnn = (N[0]+1)*(N[1]+1)*(N[2]+1);

  ofstream out;
  out.open("mesh.out");
  mesh->Print(out, Element::TETRAHEDRAL);
  out.close();

  mesh->ParaviewPrint(Element::TETRAHEDRAL);
  Array<double *> linellipticdata(1);
  Array<double *> linforcedata(1);

  linellipticdata[0] = new double[nnn];
  linforcedata[0] = new double[nnn];
  
  int index=0;
  Array<double> coord(3);
  coord[2] = 0.0;
  for(int k=0; k<=N[2]; k++)
    {
      coord[1] = 0.0;
      for(int j=0; j<=N[1]; j++)
	{
	 coord[0] = 0;
	  for(int i=0; i<=N[0]; i++)
	    {
	      linellipticdata[0][index] = evalk(coord);
	      linforcedata[0][index] = evalf(coord);
	      coord[0] += dx;
	      index++;
	    }
	  coord[1] += dy;
	}
      coord[2] += dz;
    }
 
  Data *data;
  data = new Data();
  data->SetEllipticData(linellipticdata);
  data->SetForceData(linforcedata);
 
  int nbdry = mesh->GetNBdrs();
  cout << nbdry << endl;
  Array<bool> dirichlet(nbdry);
  
  for(int j=0; j<nbdry; j++)
    {
      dirichlet[j] = false;
    }

  dirichlet[0] = true;
  Vector Dval;
  FEM *fem = new LinearFEM3D(mesh, data);
  fem->Assemble();
  fem->ProjectAFunction(evaldir, Dval);
  fem->SetDirichletBoundary(dirichlet, Dval);

  Vector sol;
  sol = Dval;
  fem->Solve(sol, 100000, 1);
  mesh->ParaviewPrint(sol, Element::TETRAHEDRAL); 
  
 delete mesh;
 delete data;
 delete []linellipticdata[0];
 delete []linforcedata[0];
 delete fem;
 
  
}

