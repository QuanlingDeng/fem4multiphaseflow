#include "mpi_fem_header.h"

//==============================================================================

LinearFEM1D_p::LinearFEM1D_p(Mesh_p *_mesh, Data *_data) :
FEM_p(_mesh, _data)
{
    strcpy(Type, "LinearFEM1D_p");
    NumGlobalDOF = mesh->GetNV();
    NumLocalDOF = 2;
    rhs.SetSize(NumGlobalDOF);
    mat = new SparseMatrix(NumGlobalDOF, 10);
}

//============================================================================== 

void LinearFEM1D_p::GetElementDOF(int i, Array<int> &ind)
{
  GetElementVertices(i, ind);
}

//==============================================================================

void LinearFEM1D_p::GetBdrElementDOF(int i, Array<int> &ind)
{
  GetBdrElementVertices(i, ind);
}

//==============================================================================

void LinearFEM1D_p::ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalEllipticCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double len = coord[0][1] - coord[0][0];
    
  locmat[0][0] = 1.0/2.0 * 1.0/len * ( locdat[1] + locdat[0] );
  locmat[1][0] = -1.0/2.0 * 1.0/len * ( locdat[1] + locdat[0] );
  locmat[0][1] = -1.0/2.0 * 1.0/len * ( locdat[1] + locdat[0] );
  locmat[1][1] = 1.0/2.0 * 1.0/len * ( locdat[1] + locdat[0] );

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void LinearFEM1D_p::ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  ;
}

//==============================================================================

void LinearFEM1D_p::ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalForceCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double len = coord[0][1] - coord[0][0];

  locrhs[0] = len * ( 1.0/6.0 * locdat[1] + 1.0/3.0* locdat[0] );
  locrhs[1] = len * ( 1.0/3.0 * locdat[1] + 1.0/6.0* locdat[0] );

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}


