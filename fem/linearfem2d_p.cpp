#include "mpi_fem_header.h"

//==============================================================================

LinearFEM2D_p::LinearFEM2D_p(Mesh_p *_mesh, Data *_data) :
FEM_p(_mesh, _data)
{
    strcpy(Type, "LinearFEM2D_p");
    NumGlobalDOF = mesh->GetNV();
    NumLocalDOF = 3;
    rhs.SetSize(NumGlobalDOF);
    mat = new SparseMatrix(NumGlobalDOF, 10);
}

//============================================================================== 

void LinearFEM2D_p::GetElementDOF(int i, Array<int> &ind)
{
  GetElementVertices(i, ind);
}

//============================================================================== 

void LinearFEM2D_p::GetBdrElementDOF(int i, Array<int> &ind)
{
  GetBdrElementVertices(i, ind);
}

//==============================================================================

void LinearFEM2D_p::ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalEllipticCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[1][2] - coord[1][0]) - (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);
  double t = ( locdat[0] + locdat[1] + locdat[2] ) / ( 6.0*detJ );

  locmat[0][0] = t * ( (coord[1][1] - coord[1][2]) * (coord[1][1] - coord[1][2]) + (coord[0][2] - coord[0][1]) * (coord[0][2] - coord[0][1]) );
  locmat[0][1] = t * ( (coord[1][1] - coord[1][2]) * (coord[1][2] - coord[1][0]) + (coord[0][2] - coord[0][1]) * (coord[0][0] - coord[0][2]) );
  locmat[0][2] = t * ( (coord[1][1] - coord[1][2]) * (coord[1][0] - coord[1][1]) + (coord[0][2] - coord[0][1]) * (coord[0][1] - coord[0][0]) );
  locmat[1][0] = locmat[0][1];
  locmat[1][1] = t * ( (coord[1][0] - coord[1][2]) * (coord[1][0] - coord[1][2]) + (coord[0][2] - coord[0][0]) * (coord[0][2] - coord[0][0]) );
  locmat[1][2] = t * ( (coord[1][2] - coord[1][0]) * (coord[1][0] - coord[1][1]) + (coord[0][0] - coord[0][2]) * (coord[0][1] - coord[0][0]) );
  locmat[2][0] = locmat[0][2];
  locmat[2][1] = locmat[1][2];
  locmat[2][2] =  t * ( (coord[1][0] - coord[1][1]) * (coord[1][0] - coord[1][1]) + (coord[0][1] - coord[0][0]) * (coord[0][1] - coord[0][0]) );
    
  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void LinearFEM2D_p::ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalReactionCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[1][2] - coord[1][0]) - (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);
  
  locmat[0][0] = detJ * ( 3.0*locdat[0] + locdat[1] + locdat[2] ) / 60.0;
  locmat[0][1] = detJ * ( 2.0*locdat[0] + 2.0*locdat[1] + locdat[2] ) / 120.0;
  locmat[0][2] = detJ * ( 2.0*locdat[0] + locdat[1] + 2.0*locdat[2] ) / 120.0;
  locmat[1][0] = locmat[0][1];
  locmat[1][1] = detJ * ( locdat[0] + 3.0*locdat[1] + locdat[2] ) / 60.0;
  locmat[1][2] = detJ * ( locdat[0] + 2.0*locdat[1] + 2.0*locdat[2] ) / 120.0;
  locmat[2][0] = locmat[0][2];
  locmat[2][1] = locmat[0][1];
  locmat[2][2] = detJ * ( locdat[0] + locdat[1] + 3.0*locdat[2] ) / 60.0;
}

//==============================================================================

void LinearFEM2D_p::ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalForceCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[1][2] - coord[1][0]) - (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);

  locrhs[0] = detJ * ( 2.0*locdat[0] + locdat[1] + locdat[2] ) / 24.0;
  locrhs[1] = detJ * ( locdat[0] + 2.0*locdat[1] + locdat[2] ) / 24.0;
  locrhs[2] = detJ * ( locdat[0] + locdat[1] + 2.0*locdat[2] ) / 24.0;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

void LinearFEM2D_p::ProjectAFunction(double (*func)(Array<double> &), Vector &pval)
{
  Array<double> coord(2);
  for(int i=0; i<mesh->GetNV(); i++)
    {
      coord[0] = mesh->GetVertex(i, 0);
      coord[1] = mesh->GetVertex(i, 1);      
      pval(i) = func(coord);
    }
}

void LinearFEM2D_p::ComputeNeumannLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalNeumannCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetBdrElement(i)->GetNVertices()];
  
  mesh->GetBdrElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[0][1] - coord[0][0]) + (coord[1][1] - coord[1][0]) * (coord[1][1] - coord[1][0]);

  detJ = sqrt(detJ);

  locrhs[0] = detJ * ( 2.0*locdat[0] + locdat[1] ) / 6.0;
  locrhs[1] = detJ * ( locdat[0] + 2.0*locdat[1] ) / 6.0;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}
