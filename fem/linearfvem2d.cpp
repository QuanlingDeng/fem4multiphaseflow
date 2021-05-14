#include "fem_header.h"
#include "../general/general_header.h"

LinearFVEM2D::LinearFVEM2D(Mesh *_mesh, Data *_data) : FEM(_mesh, _data)
{
  strcpy(Type, "LinearFVEM2D");
  NumGlobalDOF = mesh->GetNV();
  NumLocalDOF = 3;
  rhs.SetSize(NumGlobalDOF);
  mat = new SparseMatrix(NumGlobalDOF, 10);
}

//============================================================================== 

void LinearFVEM2D::GetElementDOF(int i, Array<int> &ind)
{
  GetElementVertices(i, ind);
}

//============================================================================== 

void LinearFVEM2D::GetBdrElementDOF(int i, Array<int> &ind)
{
  GetBdrElementVertices(i, ind);
}

//==============================================================================

void LinearFVEM2D::ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[1][2] - coord[1][0]) 
                 - (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);

  Array<double *> normals(3);
  Array<double> elengths;
  Array<double *> ecoord(4);
  Array<double> areas;
  for(int k=0; k<ecoord.Size(); k++) ecoord[k] = new double[2];
  for(int k=0; k<normals.Size(); k++) normals[k] = new double[2];
  
  DualMesh *dualmesh = new DualMesh(mesh);
  dualmesh->GetElemDualInfo(i, ecoord, normals, elengths, areas);  

  double dp0x, dp0y, dp1x, dp1y, dp2x, dp2y;
  dp0x = (coord[1][1]-coord[1][2])/detJ;
  dp0y = (coord[0][2]-coord[0][1])/detJ;
  dp1x = (coord[1][2]-coord[1][0])/detJ;
  dp1y = (coord[0][0]-coord[0][2])/detJ;
  dp2x = (coord[1][0]-coord[1][1])/detJ;
  dp2y = (coord[0][1]-coord[0][0])/detJ;

  int n = 4;
  Array<double *> pw(n);
  for(int k=0; k<n; k++) pw[k] = new double[2];
  GetStandLineQuadPW(n, pw);

  Array<double *> tem(3);
  for(int k=0; k<3; k++) tem[k] = new double[3];
  for(int k=0; k<3; k++)
    for(int j=0; j<3; j++)
      tem[k][j] = 0.0;
  
  double xa, xb, ya, yb, tx, ty, ttx, tty;

  Function *kk = data->GetEllipticFunction();
  for(int k=0; k<3; k++)
    {
      xa = ecoord[k][0];
      xb = ecoord[3][0];
      ya = ecoord[k][1];
      yb = ecoord[3][1];
	  
      for(int l=0; l<n; l++)
	{
	  tx = 0.5 * (xb - xa) * pw[l][0] + 0.5 * (xa + xb);
	  ty = 0.5 * (yb - ya) * pw[l][0] + 0.5 * (ya + yb);
	  
	  ttx = dp1x * (tx - coord[0][0]) + dp1y * (ty - coord[1][0]);
	  tty = dp2x * (tx - coord[0][0]) + dp2y * (ty - coord[1][0]);

	  double bcoord[2];
	  bcoord[0] = tx;
	  bcoord[1] = ty;
	  	  
	  tem[k][0] += pw[l][1] * kk->Eval(bcoord) * ( dp0x*normals[k][0] + dp0y*normals[k][1] ) * 0.5 * elengths[k];
	  tem[k][1] += pw[l][1] * kk->Eval(bcoord) * ( dp1x*normals[k][0] + dp1y*normals[k][1] ) * 0.5 * elengths[k];
	  tem[k][2] += pw[l][1] * kk->Eval(bcoord) * ( dp2x*normals[k][0] + dp2y*normals[k][1] ) * 0.5 * elengths[k];
	}
    }
  
  locmat[0][0] = tem[2][0] - tem[0][0];
  locmat[0][1] = tem[2][1] - tem[0][1];
  locmat[0][2] = tem[2][2] - tem[0][2];
      
  locmat[1][0] = tem[0][0] - tem[1][0];
  locmat[1][1] = tem[0][1] - tem[1][1];
  locmat[1][2] = tem[0][2] - tem[1][2];
  
  locmat[2][0] = tem[1][0] - tem[2][0];
  locmat[2][1] = tem[1][1] - tem[2][1];
  locmat[2][2] = tem[1][2] - tem[2][2];
          
  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
  for (int k=0; k<tem.Size(); k++) { delete []tem[k]; }
  for (int k=0; k<pw.Size(); k++) { delete []pw[k]; }
  for (int k=0; k<ecoord.Size(); k++) { delete []ecoord[k]; }
  for (int k=0; k<normals.Size(); k++) { delete []normals[k]; }
}

//==============================================================================

void LinearFVEM2D::ComputeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  ;
}

//==============================================================================

void LinearFVEM2D::ComputeConservativeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  ;
}

//==============================================================================

void LinearFVEM2D::ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  ;
}

//==============================================================================

void LinearFVEM2D::ComputeForceLocalSystem(int i, Array<int> &ind, 
					  Array<double> &locrhs)
{
  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);

  Array<double *> normals(3);
  Array<double> elengths;
  Array<double *> ecoord(4);
  Array<double> areas;
  for(int k=0; k<ecoord.Size(); k++) ecoord[k] = new double[2];
  for(int k=0; k<normals.Size(); k++) normals[k] = new double[2];

  DualMesh *dualmesh = new DualMesh(mesh);
  dualmesh->GetElemDualInfo(i, ecoord, normals, elengths, areas);  

  int np = ConvertN2PN(5);
  Array<double *> tpw(np);
  for(int k=0; k<np; k++) tpw[k] = new double[3];
  GetStandTriQuadPW(5, tpw);

  Array<double> x(3);
  Array<double> y(3);
  double t1, t2, t3, t4, detJ;

  locrhs[0] = 0.0;
  locrhs[1] = 0.0;
  locrhs[2] = 0.0;

  Function *ff = data->GetForceFunction();
  for(int l=0; l<6; l++)
    {
      if(l==0)
	{
	  x[0] = coord[0][0];
	  y[0] = coord[1][0];
	  x[1] = ecoord[0][0];
	  y[1] = ecoord[0][1];
	  x[2] = ecoord[2][0];
	  y[2] = ecoord[2][1];
	}
      else if(l==1)
	{
	  x[0] = ecoord[3][0];
	  y[0] = ecoord[3][1];
	  x[1] = ecoord[2][0];
	  y[1] = ecoord[2][1];
	  x[2] = ecoord[0][0];
	  y[2] = ecoord[0][1];
	}
      else if(l==2)
	{
	  x[0] = coord[0][1];
	  y[0] = coord[1][1];
	  x[1] = ecoord[1][0];
	  y[1] = ecoord[1][1];
	  x[2] = ecoord[0][0];
	  y[2] = ecoord[0][1];
	}
      else if(l==3)
	{
	  x[0] = ecoord[3][0];
	  y[0] = ecoord[3][1];
	  x[1] = ecoord[0][0];
	  y[1] = ecoord[0][1];
	  x[2] = ecoord[1][0];
	  y[2] = ecoord[1][1];
	}
      else if(l==4)
	{
	  x[0] = coord[0][2];
	  y[0] = coord[1][2];
	  x[1] = ecoord[2][0];
	  y[1] = ecoord[2][1];
	  x[2] = ecoord[1][0];
	  y[2] = ecoord[1][1];
	}
      else
	{
	  x[0] = ecoord[3][0];
	  y[0] = ecoord[3][1];
	  x[1] = ecoord[1][0];
	  y[1] = ecoord[1][1];
	  x[2] = ecoord[2][0];
	  y[2] = ecoord[2][1];
	}
      
      t1 = x[1] - x[0];
      t2 = x[2] - x[0];
      t3 = y[1] - y[0];
      t4 = y[2] - y[0];
      detJ = t1 * t4 - t2 * t3;
      
      for(int k=0; k<np; k++)
	{
	  double xx = tpw[k][0];
	  double yy = tpw[k][1];
	  double gpx = x[0] + t1*xx + t2*yy;
	  double gpy = y[0] + t3*xx + t4*yy;

	  double bcoord[2];
	  bcoord[0] = gpx;
	  bcoord[1] = gpy;
	  locrhs[l/2] += tpw[k][2] * ff->Eval(bcoord) * detJ;
	} 
    }

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
  for (int k=0; k<tpw.Size(); k++) { delete []tpw[k]; }
  for (int k=0; k<ecoord.Size(); k++) { delete []ecoord[k]; }
  for (int k=0; k<normals.Size(); k++) { delete []normals[k]; }
}

//==============================================================================

void LinearFVEM2D::ProjectAFunction(double (*func)(Array<double> &), Vector &pval)
{
  Array<double> coord(2);
  for(int i=0; i<mesh->GetNV(); i++)
    {
      coord[0] = mesh->GetVertex(i, 0);
      coord[1] = mesh->GetVertex(i, 1);      
      pval(i) = func(coord);
    }
}

//==============================================================================

void LinearFVEM2D::SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
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
