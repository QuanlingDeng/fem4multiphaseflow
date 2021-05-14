#include "fem_header.h"

//==============================================================================

LinearFEM1D::LinearFEM1D(Mesh *_mesh, Data *_data) :
FEM(_mesh, _data)
{
    strcpy(Type, "LinearFEM1D");
    NumGlobalDOF = mesh->GetNV();
    NumLocalDOF = 2;
    rhs.SetSize(NumGlobalDOF);
    mat = new SparseMatrix(NumGlobalDOF, 10);
}

//============================================================================== 

void LinearFEM1D::GetElementDOF(int i, Array<int> &ind)
{
  GetElementVertices(i, ind);
}

//==============================================================================

void LinearFEM1D::GetBdrElementDOF(int i, Array<int> &ind)
{
  GetBdrElementVertices(i, ind);
}

//==============================================================================

void LinearFEM1D::ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
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

void LinearFEM1D::ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  ;
}

//==============================================================================

void LinearFEM1D::ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
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

void LinearFEM1D::SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
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

