#include "fem_header.h"

//==============================================================================

CubicFEM1D::CubicFEM1D(Mesh *_mesh, Data *_data) :
FEM(_mesh, _data)
{
    strcpy(Type, "CubicFEM1D");
    NumGlobalDOF = mesh->GetNV() + 2*mesh->GetNE();
    NumLocalDOF = 4;
    rhs.SetSize(NumGlobalDOF);
    mat = new SparseMatrix(NumGlobalDOF, 10);
}

//============================================================================== 

void CubicFEM1D::GetElementDOF(int i, Array<int> &ind)
{
  ind[0] = i;
  ind[1] = (mesh->GetNE()+1)+2*i;
  ind[2] = ind[1]+1;
  ind[3] = i+1;
}

//==============================================================================

void CubicFEM1D::ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  Array<double> k(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    k[j] = data->GetNodalEllipticCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double h = (coord[0][1] - coord[0][0])/3.0;


  //row 1=====================================================================
  locmat[0][0] = (409*k[3]-1455*k[2]+4539*k[1]+4795*k[0])/(6720*h);
  locmat[0][1] = -(251*k[3]-351*k[2]+1377*k[1]+2251*k[0])/(2240*h);
  locmat[0][2] = (289*k[3]+81*k[2]-189*k[1]+827*k[0])/(2240*h);
  locmat[0][3] = -(523*k[3]-159*k[2]-159*k[1]+523*k[0])/(6720*h);

  //row 2=====================================================================
  locmat[1][0] = locmat[0][1];
  locmat[1][1] = (657*k[3]+1377*k[2]+2835*k[1]+3195*k[0])/(2240*h);
  locmat[1][2] = -(1233*k[3]+1539*k[2]+1539*k[1]+1233*k[0])/(2240*h);
  locmat[1][3] = (827*k[3]-189*k[2]+81*k[1]+289*k[0])/(2240*h);


  //row3=====================================================================
  locmat[2][0] = locmat[0][2];
  locmat[2][1] = locmat[1][2];
  locmat[2][2] = (3195*k[3]+2835*k[2]+1377*k[1]+657*k[0])/(2240*h);
  locmat[2][3] = -(2251*k[3]+1377*k[2]-351*k[1]+251*k[0])/(2240*h);

  //row4=====================================================================
  locmat[3][0] = locmat[0][3];
  locmat[3][1] = locmat[1][3];
  locmat[3][2] = locmat[2][3];
  locmat[3][3] = (4795*k[3]+4539*k[2]-1455*k[1]+409*k[0])/(6720*h);
  

  for (int k=0; k<coord.Size(); k++)
    delete []coord[k];
}

//==============================================================================

void CubicFEM1D::ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  ;
}

//==============================================================================

void CubicFEM1D::ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
{
  Array<double> f(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    f[j] = data->GetNodalForceCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double h = (coord[0][1] - coord[0][0])/3.0;

  locrhs[0] =  ((19*f[3]-36*f[2]+99*f[1]+128*f[0])*h)/560;
  locrhs[1] =  -((36*f[3]+81*f[2]-648*f[1]-99*f[0])*h)/560;
  locrhs[2] = ((99*f[3]+648*f[2]-81*f[1]-36*f[0])*h)/560;
  locrhs[3] = ((128*f[3]+99*f[2]-36*f[1]+19*f[0])*h)/560;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

void CubicFEM1D::SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
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

