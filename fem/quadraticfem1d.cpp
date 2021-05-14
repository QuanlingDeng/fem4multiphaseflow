#include "fem_header.h"

//==============================================================================

QuadraticFEM1D::QuadraticFEM1D(Mesh *_mesh, Data *_data) :
FEM(_mesh, _data)
{
    strcpy(Type, "QuadraticFEM1D");
    NumGlobalDOF = mesh->GetNV() + mesh->GetNE();
    NumLocalDOF = 3;
    rhs.SetSize(NumGlobalDOF);
    mat = new SparseMatrix(NumGlobalDOF, 10);
}

//============================================================================== 

void QuadraticFEM1D::GetElementDOF(int i, Array<int> &ind)
{
  //  GetElementVertices(i, ind);
  ind[0] = i;
  ind[2] = i + 1;
  ind[1] = mesh->GetNE() + 1 + i;
}

//==============================================================================

void LinearFEM1D::GetBdrElementDOF(int i, Array<int> &ind)
{
  GetBdrElementVertices(i, ind);
}

//==============================================================================

void QuadraticFEM1D::ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
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

void QuadraticFEM1D::ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  ;
}

//==============================================================================

void QuadraticFEM1D::ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
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

void QuadraticFEM1D::SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
			       const Array<double *> &Dval)
{
  Element *el;  
  for(int j=0; j<bdry_is_dirichlet.Size(); j++)
    {
      for (int i=0; i<GetNBE(); i++)
	{
	  el = GetBdrElement(i);
	  int attr = el->GetAttribute();
	  int nv = el->GetNVertices();
	  if (bdry_is_dirichlet[j][attr])
	    {
	      Array<int> row;
	      GetBdrElementDOF(i, row);
	      
	      // This assumes row is returned as dofs associated with first
	      // direction followed by second and so on. In 2d triangular 
	      // case it contains row[0] = dofv1, row[1] = dofv2, 
	      // row[2] = dofv1 + NV, row[3] =  dofv2 +NV. See the function in 
	      // linearfemelast2d.cpp for an example.

	      Array<int> index(nv);
	      for(int kk=0; kk<nv; kk++) { index[kk] = row[kk]; }
	      int ind = 0;
	      for (int ii=nv*j; ii<nv*j+nv; ii++)
		{
		  mat->EliminateRowCol(row[ii], Dval[j][index[ind]], rhs);
		  ind++;
		}
	    }
	}
    }
  mat->Finalize();
  // mat->Print();
}
