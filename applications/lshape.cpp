#include <iostream>
#include "../fem/fem_header.h"
#include "../mesh/mesh_header.h"
#include "../general/array.hpp"
#include <time.h>
using namespace std;

//============================================

void GetPolar(Array<double> &coord, Array<double> &polar)
{
  double r = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
  double phi = atan2(coord[1], coord[0]);
  polar[0] = r;
  polar[1] = phi;
}

//============================================

double dirichlet1(Array<double> &coord)
{
  Array<double> polar(2);
  GetPolar(coord, polar);
  double r = polar[0];
  double phi = polar[1];
  double pi = 4.0*atan(1.0);
  double E = 1e5;
  double nu = 0.3;
  double lambda = nu*E/((1+nu)*(1.0-2.0*nu));
  double mu = E/(2.0*(1.0+nu));
  double alph = 0.544483737;
  double omega = 3.0*pi/4.0;
  double C_1 = -cos((alph+1.0)*omega)/cos((alph-1)*omega);
  double C_2 = 2*(lambda+2*mu)/(lambda+mu);
  double ut = 1.0/(2.0*mu)*pow(r,alph)*((alph+1)*sin((alph+1)*phi)+(C_2+alph-1)
					*C_1*sin((alph-1)*phi));
  double ur = 1.0/(2.0*mu)*pow(r,alph)*(-(alph+1)*cos((alph+1)*phi)+(C_2-alph-1)
					*C_1*cos((alph-1)*phi)); 
 
 return ur*cos(phi)-ut*sin(phi);
}

//============================================

double dirichlet2(Array<double> &coord)
{
  Array<double> polar(2);
  GetPolar(coord, polar);
  double r = polar[0];
  double phi = polar[1];
  double pi = 4.0*atan(1.0);
  double E = 1e5;
  double nu = 0.3;
  double lambda = nu*E/((1+nu)*(1.0-2.0*nu));
  double mu = E/(2.0*(1.0+nu));
  double alph = 0.544483737;
  double omega = 3.0*pi/4.0;
  double C_1 = -cos((alph+1.0)*omega)/cos((alph-1)*omega);
  double C_2 = 2*(lambda+2*mu)/(lambda+mu);
  double ut = 1.0/(2.0*mu)*pow(r,alph)*((alph+1)*sin((alph+1)*phi)+(C_2+alph-1)
					*C_1*sin((alph-1)*phi));
  double ur = 1.0/(2.0*mu)*pow(r,alph)*(-(alph+1)*cos((alph+1)*phi)+(C_2-alph-1)
					*C_1*cos((alph-1)*phi)); 

  return ur*sin(phi)+ut*cos(phi);
}

//============================================

double neumann1(Array<double> &coord)
{
  return 0.0;
}

//==========================================

double neumann2(Array<double> &coord)
{
  return 0.0;
}

//==========================================

int main()
{
  
  int length = 256;
  char file1[length];
  char file2[length];
  char file3[length];
  
  sprintf(file1, "lshape.3.node");
  sprintf(file2, "lshape.3.ele");
  sprintf(file3, "lshape.3.poly");
  
  Array<ifstream *> files(3);
  files[0] = new ifstream(file1);
  files[1] = new ifstream(file2);
  files[2] = new ifstream(file3);
  
  Mesh *mesh = new TMesh<2>(files, Element::TRIANGLE);
  cout << "mesh created ... " << endl;
  mesh->ParaviewPrint(Element::TRIANGLE);

  for(int i=0; i<3; i++)
    delete files[i];
  
  Array<double *> edata(2);
  Array<double *> fdata(2);
  Array<double *> gdata(2);

  int nvert = mesh->GetNV();
  int dof = 2*nvert;

  Data *data = new Data();
  LinearFEMELAST2D *method = new LinearFEMELAST2D(mesh, data);

  Vector N1(nvert), N2(nvert);
  method->ProjectAFunction(neumann1, N1);
  method->ProjectAFunction(neumann2, N2);

  delete method;

  for(int i=0; i<2; i++)
    {
      edata[i] = new double[nvert];
      fdata[i] = new double[nvert];
      gdata[i] = new double[nvert];
    }

  double lambda, mu;
  double nu = 0.3;
  double E = 1e5;
  lambda = nu*E/((1+nu)*(1.0-2.0*nu));
  mu = E/(2.0*(1.0+nu));

  for(int i=0; i<nvert; i++)
    {
      edata[0][i] = lambda;
      edata[1][i] = mu;
      fdata[0][i] = 0.0;
      fdata[1][i] = 0.0;
      gdata[0][i] = N1(i);
      gdata[1][i] = N2(i);
    }

  data->SetEllipticData(edata);
  data->SetForceData(fdata);
  data->SetNeumannData(gdata);

  method = new LinearFEMELAST2D(mesh,data);
  cout << "method created ... " << endl;

  method->Assemble();
  
  //setup dirichlet boundary
  Vector D1(nvert), D2(nvert);

  method->ProjectAFunction(dirichlet1, D1);
  method->ProjectAFunction(dirichlet2, D2);

  Array<double *> dbdy(2);
  for(int i=0; i<2; i++)
    dbdy[i] = new double[nvert];

  for(int i=0; i<nvert; i++)
    {
      dbdy[0][i] = D1(i);
      dbdy[1][i] = D2(i);
    }
   
  Array<bool *> is_dirichlet(2);
  for(int i=0; i<2; i++)
    is_dirichlet[i] = new bool[mesh->GetNBdrs()];
 
  for(int i=0; i<mesh->GetNBdrs(); i++)
    {
      is_dirichlet[0][i] = true;
      is_dirichlet[1][i] = true;
    }

  method->SetDirichletBoundary(is_dirichlet, dbdy);
  Vector sol(dof);
  sol = 0;

  Array<int> ind;
  for(int i=0; i<mesh->GetNBE(); i++)
    {
      mesh->GetBdrElementVertices(i, ind);
      for(int j=0; j<ind.Size(); j++)
	{
	  sol(ind[j]) = D1(ind[j]);
	  sol(ind[j]+nvert) = D2(ind[j]);
	}
    }
  cout << "solving..." << endl;
  method->Solve(sol, "pcg");
  cout << "solved" << endl;
  ofstream out;
  method->PrintDeformedMesh(out, sol, 3000);  
  

  delete method;
  delete mesh;
  delete data;
  for(int i=0; i<2; i++)
    {
      delete []edata[i];
      delete []fdata[i];
      delete []gdata[i];
      delete []dbdy[i];
      delete []is_dirichlet[i];
    }
  

}


