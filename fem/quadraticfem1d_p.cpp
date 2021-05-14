#include "mpi_fem_header.h"

//==============================================================================

QuadraticFEM1D_p::QuadraticFEM1D_p(Mesh_p *_mesh, Data *_data) :
FEM_p(_mesh, _data)
{
    strcpy(Type, "QuadraticFEM1D_p");
    NumGlobalDOF = mesh->GetNV() + mesh->GetNE();
    NumLocalDOF = 3;
    rhs.SetSize(NumGlobalDOF);
    mat = new SparseMatrix(NumGlobalDOF, 10);
}

//============================================================================== 

void QuadraticFEM1D_p::GetElementDOF(int i, Array<int> &ind)
{
  //  GetElementVertices(i, ind);
  ind[0] = i;
  ind[2] = i + 1;
  ind[1] = mesh->GetNE() + 1 + i;
}

//==============================================================================

void LinearFEM1D_p::GetBdrElementDOF(int i, Array<int> &ind)
{
  GetBdrElementVertices(i, ind);
}

//==============================================================================

void QuadraticFEM1D_p::ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  Array<double> locd(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locd[j] = data->GetNodalEllipticCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double h = coord[0][1] - coord[0][0];
  
  locmat[0][0] = (37.0*locd[0] + 36.0*locd[1] - 3.0*locd[2])/(30.0*h);
  locmat[0][1] = (-22.0*locd[0] - 16.0*locd[1] - 2.0*locd[2])/(15.0*h);
  locmat[0][2] = (7.0*locd[0] - 4.0*locd[1] + 7.0*locd[2])/(30.0*h);
  locmat[1][1] = (24.0*locd[0] + 32.0*locd[1] + 24.0*locd[2])/(15.0*h);
  locmat[1][2] = (-2.0*locd[0] - 16.0*locd[1] - 22.0*locd[2])/(15.0*h);
  locmat[2][2] = (-3.0*locd[0] + 36.0*locd[1] + 37.0*locd[2])/(30.0*h);
  locmat[1][0] = locmat[0][1];
  locmat[2][0] = locmat[0][2];
  locmat[2][1] = locmat[1][2];

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void QuadraticFEM1D_p::ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  ;
}

//==============================================================================

void QuadraticFEM1D_p::ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
{
  Array<double> locd(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locd[j] = data->GetNodalForceCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double h = coord[0][1] - coord[0][0];

  locrhs[0] = h*(4.0*locd[0] + 2.0*locd[1] - locd[2])/30.0;
  locrhs[1] = h*(locd[0] + 8.0*locd[1] + locd[2])/15.0;
  locrhs[2] = h*(-locd[0] + 2.0*locd[1] + 4.0*locd[2])/30.0;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

void QuadraticFEM1D_p::ProjectAFunction(double (*func)(Array<double> &), Vector &pval)
{
  Array<double> coord(1);
  for(int i=0; i<mesh->GetNV(); i++)
    {
      coord[0] = mesh->GetVertex(i, 0);
      pval(i) = func(coord);
    }
}
