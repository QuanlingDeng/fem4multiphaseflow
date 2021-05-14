#include "../fem/fem_header.h"
#include "../mesh/mesh_header.h"
#include "../general/array.hpp"


double b(Array<double> &coord)
{
  if(coord[0] < 0.5)
    return 1.0;
  else 
    return 0.0;
}

int main(int argc, const char * argv[])
{
  int Nx = 120;
  int Ny = 120;

  Array<int> n(2);
  n[0] = Nx;
  n[1] = Ny;

  Array<double> L(4);
  L[0] = 0.0;
  L[1] = 1.0;
  L[2] = 0.0;
  L[3] = 1.0;

  Mesh *mesh = new TMesh<2> (n, L, Element::TRIANGLE);

  int NV = mesh->GetNV();

  Array<double *> edata(1);
  edata[0] = new double[NV];

  Array<double *> fdata(1);
  fdata[0] = new double[NV];
  
  Array<double *> ndata(1);
  ndata[0] = new double[NV];

  Vector uu(NV);

  ifstream permfile("perm_stanford_11.dat");
  double temp;
  for(int i=0; i<NV; i++)
    {
      fdata[0][i] = 0.0;
      ndata[0][i] = 0.0;
      permfile >> temp; 
      temp = exp(0.5*temp);
      edata[0][i] = temp;
    }
  
  permfile.close();
  
  Data *data = new Data();
  data->SetEllipticData(edata);
  data->SetForceData(fdata);
  data->SetNeumannData(ndata);

  cout << " data set " << endl;

  FEM *UW = new LinearFEM2D(mesh, data);
  cout << " FEM created " << endl;
  UW->Assemble();
  cout << " assembled " << endl;

  Vector Dval(NV);
  Dval = 0.0;
  UW->ProjectAFunction(b, Dval);
  cout << " projected " << endl;

  Array<bool *> bdry_is_dirichlet(1);
  bdry_is_dirichlet[0] = new bool[4];
  for(int i=0; i<4; i++)
    bdry_is_dirichlet[0][i] = true;  
  bdry_is_dirichlet[0][0] = false;
  bdry_is_dirichlet[0][2] = false;

  Array<double *> bval(1);
  bval[0] = new double[NV];
  for(int i=0; i<NV; i++)
    bval[0][i] = Dval(i);

  UW->SetDirichletBoundary(bdry_is_dirichlet, bval);
  cout << " dirichlet set " << endl;

  Vector sol(NV);
  cout << " solving... " << endl;
  UW->Solve(sol, "pcg");
  cout << " solved " << endl;

  mesh->ParaviewPrint(Element::TRIANGLE);
  mesh->ParaviewPrint(sol, Element::TRIANGLE);  
  
  delete []bdry_is_dirichlet[0];
  delete []bval[0];
  delete []ndata[0];
  delete []edata[0];
  delete []fdata[0];
  delete data;
  delete mesh;
  delete UW;

  return 0;
}
