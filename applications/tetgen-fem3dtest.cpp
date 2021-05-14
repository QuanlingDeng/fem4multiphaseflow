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
  return  fabs(sin(6.28*coord[1]));  
}

//======================================================

double evalf(Array<double> &coord)
{
  return 0.0; //cos(coord[2]);
}

//======================================================

double evaldir(Array<double> &coord)
{
  if(coord[1] < 5)
    {
      return 1.0; //sin(6.28*coord[0]) + cos(3.14*coord[1]) + coord[2]*coord[2];
    }

  return 0.0;
}

//======================================================

double evalneu(Array<double> &coord)
{

  return 0;
  /*
  double a= -10.0;
  if (coord[2] > 3.0)
    a = 10.0;

    return a;*/
}

//======================================================

int main()
{
  
  char file1[256];
  char file2[256];
  char file3[256];

  sprintf(file1, "wing.node");
  sprintf(file2, "wing.ele");
  sprintf(file3, "wing.face");

  Array<ifstream *> files;
  files.SetSize(3);
  files[0] = new ifstream( file1 );
  files[1] = new ifstream( file2 );
  files[2] = new ifstream( file3 );

  Mesh *mesh;
  mesh = new TMesh<3>(files, Element::TETRAHEDRAL);

  for(int i=0; i<3; i++)
    {
      files[i]->close();
      delete files[i];
    }
  mesh->ParaviewPrint(Element::TETRAHEDRAL);

  int nv = mesh->GetNV();
  
  Array<double *> linellipticdata(1);
  Array<double *> linforcedata(1);
  Array<double *> linneumanndata(1);

  linellipticdata[0] = new double[nv];
  linforcedata[0] = new double[nv];
  linneumanndata[0] = new double[nv];
  
  int index=0;
  for(int k=0; k<nv; k++)
    {
      Array<double> coord(3);
      for(int i=0; i<3; i++)
	{
	  coord[i] = mesh->GetVertex(k, i);
	}

      linellipticdata[0][k] = evalk(coord);
      linforcedata[0][k] = evalf(coord);
      linneumanndata[0][k] = evalneu(coord);
    }

  Data *data;
  data = new Data();
  data->SetEllipticData(linellipticdata);
  data->SetForceData(linforcedata);
  data->SetNeumannData(linneumanndata);
 
  int nbdry = mesh->GetNBdrs();
  Array<bool> dirichlet(nbdry);
 
  for(int j=0; j<nbdry; j++)
    {
      dirichlet[j] = true;
    }

  /*
  dirichlet[3] = false;
  dirichlet[7] = false;

  dirichlet[9] = false;
  dirichlet[19] = false;
  dirichlet[0] = false;
  dirichlet[1] = false;
  dirichlet[2] = false;
  */
  
  cout << "fembuild" << endl;
  FEM *fem = new LinearFEM3D(mesh, data);
  cout << "assemble" << endl; 
  fem->Assemble();

  Vector Dval;
  fem->ProjectAFunction(evaldir, Dval);
  cout << "projected" << endl;
  fem->SetDirichletBoundary(dirichlet, Dval);
  cout << "Boundaryset" << endl;
  fem->Solve(Dval, 100000, 1);

  mesh->ParaviewPrint(Dval, Element::TETRAHEDRAL); 
 
 delete mesh;
 delete data;
 delete []linellipticdata[0];
 delete []linforcedata[0];
 delete fem;

     
}

