#include "../fem/fem_header.h"
#include "../mesh/mesh_header.h"
#include "../general/array.hpp"
#include "../timedependent/time_header.h"
#include <fstream>
#include <cmath>

double b(Array<double> &coord)
{
  if(coord[0] < 0.001)
    return 1.0;
  else 
    return 0.0;
}

int main(int argc, const char * argv[])
{  
  int Nx = 15;
  int Ny = 15;

  Array<int> n(2);
  n[0] = Nx;
  n[1] = Ny;

  Array<double> L(4);
  L[0] = 0.0;
  L[1] = 1.0;
  L[2] = 0.0;
  L[3] = 1.0;

  Mesh *mesh = new TMesh<2> (n, L, Element::TRIANGLE);

  /*
  int length = 256;
  char file1[length];
  char file2[length];
  char file3[length];
  
  sprintf(file1, "unitsquare.2.node");
  sprintf(file2, "unitsquare.2.ele");
  sprintf(file3, "unitsquare.2.poly");

  Array<ifstream *> files(3);
  files[0] = new ifstream(file1);
  files[1] = new ifstream(file2);
  files[2] = new ifstream(file3);

  Mesh *mesh = new TMesh<2>(files, Element::TRIANGLE);
  cout << " mesh created" << flush;

  for(int i=0; i<3; i++)
    delete files[i];
  */


  int NV = mesh->GetNV();
  int NE = mesh->GetNE();

  Array<double *> edata(1);
  edata[0] = new double[NV];

  Array<double *> fdata(1);
  fdata[0] = new double[NV];
  
  for(int i=0; i<NV; i++)
    {
      fdata[0][i] = 0.0;
      edata[0][i] = 1.0;
    }
  
  Data *data = new Data();
  data->SetEllipticData(edata);
  data->SetForceData(fdata);
  cout << endl;
  cout << " data set " << flush;

  LinearFEM2D *method = new LinearFEM2D(mesh, data);
  cout << endl;
  cout << " FEM created " << flush;
  method->Assemble();
  cout << endl;
  cout << " assembled " << flush;

  Array<bool *> bdry_is_dirichlet(1);
  bdry_is_dirichlet[0] = new bool[4];

  for(int i=0; i<4; i++)
    {
      bdry_is_dirichlet[0][i] = true;
    }

  bdry_is_dirichlet[0][0] = false;
  bdry_is_dirichlet[0][2] = false;

  Vector D(NV);
  method->ProjectAFunction(b, D);

  Array<double *> bval(1);
  bval[0] = new double[NV];
  
  for(int i=0; i<NV; i++)
    {
      bval[0][i] = D(i);
    }

  method->SetDirichletBoundary(bdry_is_dirichlet, bval);
  cout << endl;
  cout << " dirichlet set " << flush;

  Vector sol(NV);
  cout << endl;
  cout << " solving... " << flush;

  method->Solve(sol, "pcg", 1000, 0, 1e-32, 1e-32);

  cout << endl;
  cout << " solved " << flush;
  cout << endl;

  Array<double *> flux(3);
  for(int i=0; i<3; i++)
    flux[i] = new double[NE];

  method->ComputeConservativeFlux(sol, flux);

  //  mesh->ParaviewPrint(sol, Element::TRIANGLE);
  //mesh->ParaviewPrint(Element::TRIANGLE);

  /*
  for(int i=0; i<NE; i++)
    cout << i << " " << flux[0][i] << " " << flux[1][i] << " " 
         << flux[2][i] << endl;
  */

  
  //time marching ================================
  
  ControlVolumes *cv = new ControlVolumes(mesh);
  
  Vector initsat(NV);
  method->ProjectAFunction(b, initsat);
  Array<double> initsol(NV);

  //  mesh->ParaviewPrint(initsat, Element::TRIANGLE);

  for(int i=0; i<NV; i++)
    initsol[i] = initsat(i);

  Array<double> sat(NV);
  
  double inittime = 0.0;
  double finaltime = 0.1;
  int numsteps = 10000;
  
  
  TimeMarchSolve(initsol, inittime, finaltime, numsteps, sat, cv, 
  		 mesh, flux);

  delete cv;
  

  for(int i=0; i<3; i++)
    delete []flux[i];
  

  delete []bdry_is_dirichlet[0];
  delete []bval[0];
  delete []edata[0];
  delete []fdata[0];
  delete data;
  delete mesh;
  delete method;

  return 0;
}


  /*
  //print data for gle plotting
  //======================================================================
  double *v;
  Array<double *> coord(4);
  Array<int> dofs;
  ofstream out;
  out.open("glesat.out");
  out << 3*mesh->GetNE() << endl;
  for(int i=0; i<mesh->GetNE(); i++)
    {
      mesh->GetElementVertices(i, dofs);
      cv->GetCoordinates(i, coord);
      for(int k=0; k<3; k++)
	{
	  v = mesh->GetVertex(dofs[k]);
	  if(k==0)
	    {
	      out << v[0] << " " << v[1] << " " << coord[0][0] << " " 
                  << coord[0][1] << " " << coord[3][0] << " " << coord[3][1] 
                  << " " << coord[2][0] << " " << coord[2][1] << " " 
                  << sat[dofs[k]] << endl;
	    }
	  else if(k==1)
	    {
	      out << v[0] << " " << v[1] << " " << coord[1][0] << " " 
                  << coord[1][1] << " "  << coord[3][0] << " " << coord[3][1] 
                  << " " << coord[0][0] << " " << coord[0][1] << " " 
                  << sat[dofs[k]] << endl;
	    }
	  else if(k==2)
	    {
	      out << v[0] << " " << v[1] << " " << coord[2][0] << " " 
                  << coord[2][1] << " " << coord[3][0] << " " << coord[3][1] 
                  << " " << coord[1][0] << " " << coord[1][1] << " " 
                  << sat[dofs[k]] << endl;
	    }
	}
    }
  out.close();
  //=========================================================================
*/
