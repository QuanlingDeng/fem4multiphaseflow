#include "fem_header.h"
#include "../general/general_header.h"

CubicFEM2D::CubicFEM2D(Mesh *_mesh, Data *_data) : FEM(_mesh, _data)
{
  strcpy(Type, "CubicFEM2D");
  NumGlobalDOF = mesh->GetNV() + 2 * mesh->GetNEdges() + mesh->GetNE();
  NumLocalDOF = 10;
  rhs.SetSize(NumGlobalDOF);
  mat = new SparseMatrix(NumGlobalDOF, 60);
}

//============================================================================== 

void CubicFEM2D::GetElementDOF(int i, Array<int> &ind)
{
  GetElementVertices(i, ind);

  Array<int> edges;
  mesh->GetElementEdges(i, edges);

  if(ind[0]<ind[1] && ind[0]<ind[2])
    {
      if(ind[1] < ind[2])
	{
	  ind.Append( mesh->GetNV() + 2*edges[0] );
	  ind.Append( mesh->GetNV() + 2*edges[0] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[1] );
	  ind.Append( mesh->GetNV() + 2*edges[1] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[2] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[2] );
	}
      else
	{
	  ind.Append( mesh->GetNV() + 2*edges[0] );
	  ind.Append( mesh->GetNV() + 2*edges[0] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[1] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[1] );
	  ind.Append( mesh->GetNV() + 2*edges[2] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[2] );
	}
    }
  else if(ind[1]<ind[0] && ind[1]<ind[2])
    {
      if(ind[0] < ind[2])
	{
	  ind.Append( mesh->GetNV() + 2*edges[0] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[0] );
	  ind.Append( mesh->GetNV() + 2*edges[1] );
	  ind.Append( mesh->GetNV() + 2*edges[1] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[2] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[2] );
	}
      else
	{
	  ind.Append( mesh->GetNV() + 2*edges[0] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[0] );
	  ind.Append( mesh->GetNV() + 2*edges[1] );
	  ind.Append( mesh->GetNV() + 2*edges[1] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[2] );
	  ind.Append( mesh->GetNV() + 2*edges[2] + 1);
	}
    }
  else
    {
      if(ind[0] < ind[1])
	{
	  ind.Append( mesh->GetNV() + 2*edges[0] );
	  ind.Append( mesh->GetNV() + 2*edges[0] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[1] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[1] );
	  ind.Append( mesh->GetNV() + 2*edges[2] );
	  ind.Append( mesh->GetNV() + 2*edges[2] + 1);
	}
      else
	{
	  ind.Append( mesh->GetNV() + 2*edges[0] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[0] );
	  ind.Append( mesh->GetNV() + 2*edges[1] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[1] );
	  ind.Append( mesh->GetNV() + 2*edges[2] );
	  ind.Append( mesh->GetNV() + 2*edges[2] + 1);
	}
    }
  
  ind.Append( mesh->GetNV() + 2 * mesh->GetNEdges() + i );
  //for(int k=0; k<ind.Size(); k++) cout<<ind[k] <<endl; cout<<endl;
}

//============================================================================== 

void CubicFEM2D::GetBdrElementDOF(int i, Array<int> &ind)
{
  GetBdrElementVertices(i, ind);

  if(ind[0] < ind[1])
    {
      ind.Append( mesh->GetNV() + 2 * mesh->GetBdrElementEdgeIndex(i) );
      ind.Append( mesh->GetNV() + 2 * mesh->GetBdrElementEdgeIndex(i) + 1);
    }
  else
    {
      ind.Append( mesh->GetNV() + 2 * mesh->GetBdrElementEdgeIndex(i) + 1);
      ind.Append( mesh->GetNV() + 2 * mesh->GetBdrElementEdgeIndex(i) );
    }

  //for(int k=0; k<ind.Size(); k++) cout<<ind[k] <<endl; cout<<endl;
}

//==============================================================================

void CubicFEM2D::ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  for(int k=0; k<locmat.Size(); k++)
    for(int j=0; j<locmat.Size(); j++)
      locmat[k][j] = 0.0;
	  
  int np = ConvertN2PN(6);
  Array<double *> tpw(np);
  for(int j=0; j<np; j++) tpw[j] = new double[3];
  GetStandTriQuadPW(6, tpw);

  Array<double *> c(mesh->GetDim()); // c for coordinate
  for (int j=0; j<c.Size(); j++)
    c[j] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, c);

  double t1 = c[0][1] - c[0][0];
  double t2 = c[0][2] - c[0][0];
  double t3 = c[1][1] - c[1][0];
  double t4 = c[1][2] - c[1][0];
  double detJ = t1 * t4 - t2 * t3;
    
  // Jacobian Inverse Transpose * detJ
  Array<double *> jit(2);
  jit[0] = new double[2];
  jit[1] = new double[2];
  jit[0][0] = c[1][2] - c[1][0];
  jit[0][1] = c[1][0] - c[1][1];
  jit[1][0] = c[0][0] - c[0][2];
  jit[1][1] = c[0][1] - c[0][0];

  double gpx, gpy, xx, yy, dotprod, vm[2], vn[2];

  BasisFunctions *phi = new BasisFunctions(3);
  Function *kk = data->GetEllipticFunction();
  
  for(int m=0; m<locmat.Size(); m++)
    {
      for(int n=m; n<locmat.Size(); n++)
	{
	  for(int k=0; k<np; k++)
	    {
	      gpx = tpw[k][0];
	      gpy = tpw[k][1];
	      xx = c[0][0] + t1*gpx + t2*gpy;
	      yy = c[1][0] + t3*gpx + t4*gpy;

	      vn[0] = jit[0][0] * phi->GradBF2D(n, 0, gpx, gpy) + jit[0][1] * phi->GradBF2D(n, 1, gpx, gpy);
	      vn[1] = jit[1][0] * phi->GradBF2D(n, 0, gpx, gpy) + jit[1][1] * phi->GradBF2D(n, 1, gpx, gpy);

	      vm[0] = jit[0][0] * phi->GradBF2D(m, 0, gpx, gpy) + jit[0][1] * phi->GradBF2D(m, 1, gpx, gpy);
	      vm[1] = jit[1][0] * phi->GradBF2D(m, 0, gpx, gpy) + jit[1][1] * phi->GradBF2D(m, 1, gpx, gpy);

	      dotprod = vn[0] * vm[0] + vn[1] * vm[1];

	      double coord[2];
	      coord[0] = xx;
	      coord[1] = yy;
	      locmat[m][n] += tpw[k][2] * kk->Eval(coord) * dotprod / detJ;
	    } 	  	  
	}
    }
  
  for(int m=0; m<locmat.Size(); m++)
    {
      for(int n=m+1; n<locmat.Size(); n++)
	{
	  locmat[n][m] = locmat[m][n];
	}
    }

  for (int k=0; k<jit.Size(); k++) { delete []jit[k]; }
  for (int k=0; k<c.Size(); k++) { delete []c[k]; }
}

//==============================================================================

void CubicFEM2D::ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  ;
}

//==============================================================================

void CubicFEM2D::ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
{
  for(int j=0; j<locrhs.Size(); j++) locrhs[j] = 0.0;

  int np = ConvertN2PN(6);
  Array<double *> tpw(np);
  for(int j=0; j<np; j++) tpw[j] = new double[3];
  GetStandTriQuadPW(6, tpw);

  Array<double *> c(mesh->GetDim()); // c for coordinate
  for (int j=0; j<c.Size(); j++)
    c[j] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, c);

  double t1 = c[0][1] - c[0][0];
  double t2 = c[0][2] - c[0][0];
  double t3 = c[1][1] - c[1][0];
  double t4 = c[1][2] - c[1][0];
  double detJ = t1 * t4 - t2 * t3;
  
  double gpx, gpy, xx, yy;

  BasisFunctions *phi = new BasisFunctions(3);
  Function *ff = data->GetForceFunction();
  for(int m=0; m<locrhs.Size(); m++)
    {
      for(int k=0; k<np; k++)
	{
	  gpx = tpw[k][0];
	  gpy = tpw[k][1];
	  xx = c[0][0] + t1*gpx + t2*gpy;
	  yy = c[1][0] + t3*gpx + t4*gpy;

	  double coord[2];
	  coord[0] = xx;
	  coord[1] = yy;
	  locrhs[m] += tpw[k][2] * ff->Eval(coord) * phi->BF2D(m, gpx, gpy) * detJ;
	}
    }
  
  for (int k=0; k<c.Size(); k++) { delete []c[k]; }
}

//==============================================================================

void CubicFEM2D::SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
			       const Array<double *> &Dval)
{
  Element *el;  
  for(int j=0; j<bdry_is_dirichlet.Size(); j++)
    {
      for (int i=0; i<GetNBE(); i++)
	{
	  el = GetBdrElement(i);
	  int attr = el->GetAttribute();
	  int nv = el->GetNVertices() + 2;
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
  //mat->Print();
  //rhs.Print();
}

//==============================================================================

void CubicFEM2D::ProjectAFunction(double (*func)(Array<double> &), Vector &pval)
{
  ;
}

//==============================================================================

double CubicFEM2D::ComputeL2Error(Function *exactsol, Vector &femsol, double time)
{
  int N = 6; // accuracy of Gaussian Quadrature
  int np = ConvertN2PN(N);
  Array<double *> pw(np);
  for(int i=0; i<np; i++) pw[i] = new double[3];
  GetStandTriQuadPW(N, pw);  
  
  BasisFunctions *phi = new BasisFunctions(3);
  double xx, yy, err = 0.0; 
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      GetElementDOF(i, ind);

      Array<double> x(3);
      Array<double> y(3);
      for(int k=0; k<3; k++)
	{
	  x[k] = mesh->GetVertex(ind[k], 0);
	  y[k] = mesh->GetVertex(ind[k], 1);	  
	}

      double t1 = x[1] - x[0];
      double t2 = x[2] - x[0];
      double t3 = y[1] - y[0];
      double t4 = y[2] - y[0];
      double detJ = t1 * t4 - t2 * t3;

      double coord[3];
      coord[2] = time;
      for(int k=0; k<np; k++)
	{
	  xx = pw[k][0];
	  yy = pw[k][1];
	  coord[0] = x[0] + t1*xx + t2*yy;
	  coord[1] = y[0] + t3*xx + t4*yy;

	  double uh = 0.0;
	  for(int j=0; j<10; j++) uh += femsol(ind[j]) * phi->BF2D(j, xx, yy);
	  uh -= exactsol->Eval(coord);
	  err += pw[k][2] * uh * uh * detJ;
	}
    }  
  return sqrt(err);
}

//==============================================================================

void CubicFEM2D::ComputePPL2Error(Array<double> &errs, Function *exactsol, Vector &femsol, Array<double *> &ppsol, double time)
{
  errs[0] = 0.0;
  errs[1] = 0.0;
  errs[2] = 0.0;
  
  int N = 6; // accuracy of Gaussian Quadrature
  int np = ConvertN2PN(N);
  Array<double *> pw(np);
  for(int i=0; i<np; i++) pw[i] = new double[3];
  GetStandTriQuadPW(N, pw);  
  
  BasisFunctions *phi = new BasisFunctions(3);
  double xx, yy; 
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      GetElementDOF(i, ind);

      Array<double> x(3);
      Array<double> y(3);
      for(int k=0; k<3; k++)
	{
	  x[k] = mesh->GetVertex(ind[k], 0);
	  y[k] = mesh->GetVertex(ind[k], 1);	  
	}

      double t1 = x[1] - x[0];
      double t2 = x[2] - x[0];
      double t3 = y[1] - y[0];
      double t4 = y[2] - y[0];
      double detJ = t1 * t4 - t2 * t3;

      double coord[3];
      coord[2] = time;
      for(int k=0; k<np; k++)
	{
	  xx = pw[k][0];
	  yy = pw[k][1];
	  coord[0] = x[0] + t1*xx + t2*yy;
	  coord[1] = y[0] + t3*xx + t4*yy;

	  double uu = exactsol->Eval(coord);
	  double femuh = 0.0;
	  for(int j=0; j<10; j++) femuh += femsol(ind[j]) * phi->BF2D(j, xx, yy);
	  double ppuh = 0.0;
	  for(int j=0; j<10; j++) ppuh += ppsol[j][i] * phi->BF2D(j, xx, yy);
	  
	  double tem = uu - femuh;
	  errs[0] += pw[k][2] * tem * tem * detJ;
	  tem = uu - ppuh;
	  errs[1] += pw[k][2] * tem * tem * detJ;
	  tem = ppuh - femuh;
	  errs[2] += pw[k][2] * tem * tem * detJ;
	}
    }

  errs[0] = sqrt(errs[0]);
  errs[1] = sqrt(errs[1]);
  errs[2] = sqrt(errs[2]);
}

//==============================================================================

double CubicFEM2D::ComputeH1Error(Array<Function *> &exactgrad, Vector &femsol)
{
  int N = 6; // accuracy of Gaussian Quadrature
  int np = ConvertN2PN(N);
  Array<double *> pw(np);
  for(int i=0; i<np; i++) pw[i] = new double[3];
  GetStandTriQuadPW(N, pw);  
  
  Array<double *> jit(2);
  jit[0] = new double[2];
  jit[1] = new double[2];

  double xx, yy, err = 0.0;
  BasisFunctions *phi = new BasisFunctions(3);
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      GetElementDOF(i, ind);

      Array<double> x(3);
      Array<double> y(3);
      for(int k=0; k<3; k++)
	{
	  x[k] = mesh->GetVertex(ind[k], 0);
	  y[k] = mesh->GetVertex(ind[k], 1);	  
	}

      double t1 = x[1] - x[0];
      double t2 = x[2] - x[0];
      double t3 = y[1] - y[0];
      double t4 = y[2] - y[0];
      double detJ = t1 * t4 - t2 * t3;

      jit[0][0] = y[2] - y[0];
      jit[0][1] = y[0] - y[1];
      jit[1][0] = x[0] - x[2];
      jit[1][1] = x[1] - x[0];

      for(int k=0; k<np; k++)
	{
	  xx = pw[k][0];
	  yy = pw[k][1];
	  double coord[2];
	  coord[0] = x[0] + t1*xx + t2*yy;
	  coord[1] = y[0] + t3*xx + t4*yy;

	  double uhx = 0.0;
	  double uhy = 0.0;
	  for(int j=0; j<10; j++)
	    {
	      uhx += femsol(ind[j]) * ( jit[0][0] * phi->GradBF2D(j, 0, xx, yy) + jit[0][1] * phi->GradBF2D(j, 1, xx, yy) );
	      uhy += femsol(ind[j]) * ( jit[1][0] * phi->GradBF2D(j, 0, xx, yy) + jit[1][1] * phi->GradBF2D(j, 1, xx, yy) );
	    }

	  uhx = exactgrad[0]->Eval(coord) - uhx / detJ;
	  uhy = exactgrad[1]->Eval(coord) - uhy / detJ;
	  err += pw[k][2] * (uhx*uhx + uhy*uhy) * detJ;
	}
    }
  return sqrt(err);
}

//==============================================================================

void CubicFEM2D::ComputePPH1Error(Array<double> &errs, Array<Function *> &exactgrad, Vector &femsol, Array<double *> &ppsol)
{
  errs[0] = 0.0;
  errs[1] = 0.0;
  errs[2] = 0.0;
  
  int N = 6; // accuracy of Gaussian Quadrature
  int np = ConvertN2PN(N);
  Array<double *> pw(np);
  for(int i=0; i<np; i++) pw[i] = new double[3];
  GetStandTriQuadPW(N, pw);  
  
  Array<double *> jit(2);
  jit[0] = new double[2];
  jit[1] = new double[2];

  double xx, yy;
  BasisFunctions *phi = new BasisFunctions(3);
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      GetElementDOF(i, ind);

      Array<double> x(3);
      Array<double> y(3);
      for(int k=0; k<3; k++)
	{
	  x[k] = mesh->GetVertex(ind[k], 0);
	  y[k] = mesh->GetVertex(ind[k], 1);	  
	}

      double t1 = x[1] - x[0];
      double t2 = x[2] - x[0];
      double t3 = y[1] - y[0];
      double t4 = y[2] - y[0];
      double detJ = t1 * t4 - t2 * t3;

      jit[0][0] = y[2] - y[0];
      jit[0][1] = y[0] - y[1];
      jit[1][0] = x[0] - x[2];
      jit[1][1] = x[1] - x[0];

      for(int k=0; k<np; k++)
	{
	  xx = pw[k][0];
	  yy = pw[k][1];
	  double coord[2];
	  coord[0] = x[0] + t1*xx + t2*yy;
	  coord[1] = y[0] + t3*xx + t4*yy;

	  double ux = exactgrad[0]->Eval(coord);
	  double uy = exactgrad[1]->Eval(coord);

	  double femuhx = 0.0;
	  double femuhy = 0.0;
	  for(int j=0; j<10; j++)
	    {
	      femuhx += femsol(ind[j]) * ( jit[0][0] * phi->GradBF2D(j, 0, xx, yy) + jit[0][1] * phi->GradBF2D(j, 1, xx, yy) );
	      femuhy += femsol(ind[j]) * ( jit[1][0] * phi->GradBF2D(j, 0, xx, yy) + jit[1][1] * phi->GradBF2D(j, 1, xx, yy) );
	    }

	  double ppuhx = 0.0;
	  double ppuhy = 0.0;
	  for(int j=0; j<10; j++)
	    {
	      ppuhx += ppsol[j][i] * ( jit[0][0] * phi->GradBF2D(j, 0, xx, yy) + jit[0][1] * phi->GradBF2D(j, 1, xx, yy) );
	      ppuhy += ppsol[j][i] * ( jit[1][0] * phi->GradBF2D(j, 0, xx, yy) + jit[1][1] * phi->GradBF2D(j, 1, xx, yy) );
	    }

	  double temx = ux - femuhx / detJ;
	  double temy = uy - femuhy / detJ;
	  errs[0] += pw[k][2] * (temx*temx + temy*temy) * detJ;

	  temx = ux - ppuhx / detJ;
	  temy = uy - ppuhy / detJ;
	  errs[1] += pw[k][2] * (temx*temx + temy*temy) * detJ;
	  
	  temx = (femuhx - ppuhx) / detJ;
	  temy = (femuhy - ppuhy) / detJ;
	  errs[2] += pw[k][2] * (temx*temx + temy*temy) * detJ;
	}
    }
  
  errs[0] = sqrt(errs[0]);
  errs[1] = sqrt(errs[1]);
  errs[2] = sqrt(errs[2]);
}
