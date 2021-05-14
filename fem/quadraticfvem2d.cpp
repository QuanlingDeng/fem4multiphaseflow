#include "fem_header.h"
#include "../general/general_header.h"

QuadraticFVEM2D::QuadraticFVEM2D(Mesh *_mesh, Data *_data) :
FEM(_mesh, _data)
{
  strcpy(Type, "QuadraticFVEM2D");
  NumGlobalDOF = mesh->GetNV() + mesh->GetNEdges();
  NumLocalDOF = 6;
  rhs.SetSize(NumGlobalDOF);
  mat = new SparseMatrix(NumGlobalDOF, 40);
}

//============================================================================== 

void QuadraticFVEM2D::GetElementDOF(int i, Array<int> &ind)
{
  GetElementVertices(i, ind);

  Array<int> edges;
  mesh->GetElementEdges(i, edges);
  for (int k = 0; k < edges.Size(); k++)
    ind.Append (mesh->GetNV() + edges[k]);
}

//============================================================================== 

void QuadraticFVEM2D::GetBdrElementDOF(int i, Array<int> &ind)
{
  GetBdrElementVertices(i, ind);
  ind.Append(mesh->GetNV() +  mesh->GetBdrElementEdgeIndex(i));
}

//==============================================================================

void QuadraticFVEM2D::ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  Array<double *> coord(mesh->GetDim());
  for (int j=0; j<coord.Size(); j++)
    coord[j] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[1][2] - coord[1][0]) -
    (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);

  Array<double *> normals(3);
  Array<double> elengths;
  Array<double *> ecoord(4);
  Array<double> areas;
  for(int k=0; k<ecoord.Size(); k++) ecoord[k] = new double[2];
  for(int k=0; k<normals.Size(); k++) normals[k] = new double[2];
  
  DualMesh *dualmesh = new DualMesh(mesh);
  dualmesh->GetElemDualInfo(i, ecoord, normals, elengths, areas);  

  Array<double *> qnormals(3);
  Array<double> qelengths;
  Array<double> qareas;
  Array<double *> qecoord(13);
  for(int k=0; k<qecoord.Size(); k++) qecoord[k] = new double[2];
  for(int k=0; k<qnormals.Size(); k++) qnormals[k] = new double[2];

  DualMesh *qdualmesh = new DualMesh(mesh, 2, 0);
  qdualmesh->GetElemDualInfo(i, qecoord, qnormals, qelengths, qareas);  
  
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
  for(int k=0; k<3; k++) tem[k] = new double[6];
  for(int k=0; k<3; k++)
    for(int j=0; j<6; j++)
      tem[k][j] = 0.0;
  
  double xa, xb, ya, yb, tx, ty, ttx, tty;

  BasisFunctions *phi = new BasisFunctions(2);
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
	  tem[k][0] += -pw[l][1] * kk->Eval(bcoord) * phi->GradBF2D(0, 0, ttx, tty) * ( dp0x*normals[k][0] + dp0y*normals[k][1] ) * 0.5 * elengths[k];
	  tem[k][1] += pw[l][1] * kk->Eval(bcoord) *  phi->GradBF2D(1, 0, ttx, tty) * ( dp1x*normals[k][0] + dp1y*normals[k][1] ) * 0.5 * elengths[k];
	  tem[k][2] += pw[l][1] * kk->Eval(bcoord) *  phi->GradBF2D(2, 1, ttx, tty) * ( dp2x*normals[k][0] + dp2y*normals[k][1] ) * 0.5 * elengths[k];
	  
	  tem[k][3] += -pw[l][1] * kk->Eval(bcoord) * ( phi->GradBF2D(3, 0, ttx, tty)  * (dp1x*normals[k][0] + dp1y*normals[k][1]) + phi->GradBF2D(4, 0, ttx, tty)  * (dp2x*normals[k][0] + dp2y*normals[k][1]) ) * 0.5 * elengths[k];
	  
	  tem[k][4] += pw[l][1] * kk->Eval(bcoord) * ( phi->GradBF2D(4, 1, ttx, tty)  * (dp1x*normals[k][0] + dp1y*normals[k][1]) + phi->GradBF2D(4, 0, ttx, tty)  * (dp2x*normals[k][0] + dp2y*normals[k][1]) ) * 0.5 * elengths[k];
	  
	  tem[k][5] += -4.0 * pw[l][1] * kk->Eval(bcoord) * ( phi->GradBF2D(4, 1, ttx, tty)  * (dp1x*normals[k][0] + dp1y*normals[k][1]) + phi->GradBF2D(5, 1, ttx, tty)  * (dp2x*normals[k][0] + dp2y*normals[k][1]) ) * 0.5 * elengths[k];
	}
    }

  Array<double *> qtem(12);
  for(int k=0; k<12; k++) qtem[k] = new double[6];
  for(int k=0; k<12; k++)
    for(int j=0; j<6; j++)
      qtem[k][j] = 0.0;
  
  for(int k=0; k<12; k++)
    {
      int rm = k;
      if(k<3)
	{
	  xa = qecoord[k][0];
	  xb = qecoord[3][0];
	  ya = qecoord[k][1];
	  yb = qecoord[3][1];
	}
      else if(k<6)
	{
	  //rm = (k+1) % 3;
	  xa = qecoord[k+1][0];
	  xb = qecoord[7][0];
	  ya = qecoord[k+1][1];
	  yb = qecoord[7][1];
	}
      else if(k<9)
	{
	  //rm = (k-1) % 3;
	  xa = qecoord[k+2][0];
	  xb = qecoord[11][0];
	  ya = qecoord[k+2][1];
	  yb = qecoord[11][1];
	}
      else
	{
	  //rm = k % 3;
	  if(k==9)
	    {
	      xa = qecoord[9][0];		  
	      ya = qecoord[9][1];
	    }
	  else if(k==10)
	    {
	      xa = qecoord[1][0];		  
	      ya = qecoord[1][1];
	    }
	  else
	    {
	      xa = qecoord[5][0];		  
	      ya = qecoord[5][1];
	    }
	  xb = qecoord[12][0];
	  yb = qecoord[12][1];
	}
      
      for(int l=0; l<n; l++)
	{
	  tx = 0.5 * (xb - xa) * pw[l][0] + 0.5 * (xa + xb);
	  ty = 0.5 * (yb - ya) * pw[l][0] + 0.5 * (ya + yb);
	  
	  ttx = dp1x * (tx - coord[0][0]) + dp1y * (ty - coord[1][0]);
	  tty = dp2x * (tx - coord[0][0]) + dp2y * (ty - coord[1][0]);
	  
	  double bcoord[2];
	  bcoord[0] = tx;
	  bcoord[1] = ty;
	  qtem[k][0] += -pw[l][1] * kk->Eval(bcoord) * phi->GradBF2D(0, 0, ttx, tty) *
	    ( dp0x*qnormals[rm][0] + dp0y*qnormals[rm][1] ) * 0.5 * qelengths[rm];
	  
	  qtem[k][1] += pw[l][1] * kk->Eval(bcoord) * phi->GradBF2D(1, 0, ttx, tty) *
	    ( dp1x*qnormals[rm][0] + dp1y*qnormals[rm][1] ) * 0.5 * qelengths[rm];
	  
	  qtem[k][2] += pw[l][1] * kk->Eval(bcoord) * phi->GradBF2D(2, 1, ttx, tty) *
	    ( dp2x*qnormals[rm][0] + dp2y*qnormals[rm][1] ) * 0.5 * qelengths[rm];
	  
	  qtem[k][3] += -pw[l][1] * kk->Eval(bcoord) * ( phi->GradBF2D(3, 0, ttx, tty) * (dp1x*qnormals[rm][0] + dp1y*qnormals[rm][1]) + phi->GradBF2D(4, 0, ttx, tty) * (dp2x*qnormals[rm][0] + dp2y*qnormals[rm][1]) ) * 0.5 * qelengths[rm];
	  
	  qtem[k][4] += 4.0 * pw[l][1] * kk->Eval(bcoord) * ( phi->GradBF2D(4, 1, ttx, tty)  * (dp1x*qnormals[rm][0] + dp1y*qnormals[rm][1]) + phi->GradBF2D(4, 0, ttx, tty) * (dp2x*qnormals[rm][0] + dp2y*qnormals[rm][1]) ) * 0.5 * qelengths[rm];
	  
	  qtem[k][5] += -4.0 * pw[l][1] * kk->Eval(bcoord) * ( phi->GradBF2D(4, 1, ttx, tty) * (dp1x*qnormals[rm][0] + dp1y*qnormals[rm][1]) + phi->GradBF2D(5, 1, ttx, tty) * (dp2x*qnormals[rm][0] + dp2y*qnormals[rm][1]) ) * 0.5 * qelengths[rm];
	}
    }

  /*
  locmat[0][0] = tem[2][0] - tem[0][0];
  locmat[0][1] = tem[2][1] - tem[0][1];
  locmat[0][2] = tem[2][2] - tem[0][2];
  locmat[0][3] = tem[2][3] - tem[0][3];
  locmat[0][4] = tem[2][4] - tem[0][4];
  locmat[0][5] = tem[2][5] - tem[0][5];
      
  locmat[1][0] = tem[0][0] - tem[1][0];
  locmat[1][1] = tem[0][1] - tem[1][1];
  locmat[1][2] = tem[0][2] - tem[1][2];
  locmat[1][3] = tem[0][3] - tem[1][3];
  locmat[1][4] = tem[0][4] - tem[1][4];
  locmat[1][5] = tem[0][5] - tem[1][5];
  
  locmat[2][0] = tem[1][0] - tem[2][0];
  locmat[2][1] = tem[1][1] - tem[2][1];
  locmat[2][2] = tem[1][2] - tem[2][2]; 
  locmat[2][3] = tem[1][3] - tem[2][3];
  locmat[2][4] = tem[1][4] - tem[2][4];
  locmat[2][5] = tem[1][5] - tem[2][5];
  */

  locmat[0][0] = qtem[2][0] - qtem[0][0];
  locmat[0][1] = qtem[2][1] - qtem[0][1];
  locmat[0][2] = qtem[2][2] - qtem[0][2];
  locmat[0][3] = qtem[2][3] - qtem[0][3];
  locmat[0][4] = qtem[2][4] - qtem[0][4];
  locmat[0][5] = qtem[2][5] - qtem[0][5];
      
  locmat[1][0] = qtem[5][0] - qtem[3][0];
  locmat[1][1] = qtem[5][1] - qtem[3][1];
  locmat[1][2] = qtem[5][2] - qtem[3][2];
  locmat[1][3] = qtem[5][3] - qtem[3][3];
  locmat[1][4] = qtem[5][4] - qtem[3][4];
  locmat[1][5] = qtem[5][5] - qtem[3][5];
  
  locmat[2][0] = qtem[8][0] - qtem[6][0];
  locmat[2][1] = qtem[8][1] - qtem[6][1];
  locmat[2][2] = qtem[8][2] - qtem[6][2]; 
  locmat[2][3] = qtem[8][3] - qtem[6][3];
  locmat[2][4] = qtem[8][4] - qtem[6][4];
  locmat[2][5] = qtem[8][5] - qtem[6][5];

  locmat[3][0] = qtem[0][0] - qtem[1][0] + qtem[4][0] - qtem[5][0] + qtem[11][0] - qtem[10][0];
  locmat[3][1] = qtem[0][1] - qtem[1][1] + qtem[4][1] - qtem[5][1] + qtem[11][1] - qtem[10][1];
  locmat[3][2] = qtem[0][2] - qtem[1][2] + qtem[4][2] - qtem[5][2] + qtem[11][2] - qtem[10][2];
  locmat[3][3] = qtem[0][3] - qtem[1][3] + qtem[4][3] - qtem[5][3] + qtem[11][3] - qtem[10][3];
  locmat[3][4] = qtem[0][4] - qtem[1][4] + qtem[4][4] - qtem[5][4] + qtem[11][4] - qtem[10][4];
  locmat[3][5] = qtem[0][5] - qtem[1][5] + qtem[4][5] - qtem[5][5] + qtem[11][5] - qtem[10][5];
  
  locmat[4][0] = qtem[3][0] - qtem[4][0] + qtem[7][0] - qtem[8][0] + qtem[9][0] - qtem[11][0];
  locmat[4][1] = qtem[3][1] - qtem[4][1] + qtem[7][1] - qtem[8][1] + qtem[9][1] - qtem[11][1];
  locmat[4][2] = qtem[3][2] - qtem[4][2] + qtem[7][2] - qtem[8][2] + qtem[9][2] - qtem[11][2];
  locmat[4][3] = qtem[3][3] - qtem[4][3] + qtem[7][3] - qtem[8][3] + qtem[9][3] - qtem[11][3];
  locmat[4][4] = qtem[3][4] - qtem[4][4] + qtem[7][4] - qtem[8][4] + qtem[9][4] - qtem[11][4];
  locmat[4][5] = qtem[3][5] - qtem[4][5] + qtem[7][5] - qtem[8][5] + qtem[9][5] - qtem[11][5];
  
  locmat[5][0] = qtem[6][0] - qtem[7][0] + qtem[1][0] - qtem[2][0] + qtem[10][0] - qtem[9][0];
  locmat[5][1] = qtem[6][1] - qtem[7][1] + qtem[1][1] - qtem[2][1] + qtem[10][1] - qtem[9][1];
  locmat[5][2] = qtem[6][2] - qtem[7][2] + qtem[1][2] - qtem[2][2] + qtem[10][2] - qtem[9][2];
  locmat[5][3] = qtem[6][3] - qtem[7][3] + qtem[1][3] - qtem[2][3] + qtem[10][3] - qtem[9][3];
  locmat[5][4] = qtem[6][4] - qtem[7][4] + qtem[1][4] - qtem[2][4] + qtem[10][4] - qtem[9][4];
  locmat[5][5] = qtem[6][5] - qtem[7][5] + qtem[1][5] - qtem[2][5] + qtem[10][5] - qtem[9][5];
          
  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
  for (int k=0; k<tem.Size(); k++) { delete []tem[k]; }
  for (int k=0; k<qtem.Size(); k++) { delete []qtem[k]; }
  for (int k=0; k<pw.Size(); k++) { delete []pw[k]; }

  for (int k=0; k<ecoord.Size(); k++) { delete []ecoord[k]; }
  for (int k=0; k<normals.Size(); k++) { delete []normals[k]; }
  for (int k=0; k<qecoord.Size(); k++) { delete []qecoord[k]; }
  for (int k=0; k<qnormals.Size(); k++) { delete []qnormals[k]; }
}

//==============================================================================

void QuadraticFVEM2D::ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  ;
}

//==============================================================================

void QuadraticFVEM2D::ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
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

  Array<double *> qnormals(3);
  Array<double> qelengths;
  Array<double> qareas;
  Array<double *> qecoord(13);
  for(int k=0; k<qecoord.Size(); k++) qecoord[k] = new double[2];
  for(int k=0; k<qnormals.Size(); k++) qnormals[k] = new double[2];
  DualMesh *qdualmesh = new DualMesh(mesh, 2, 0);
  qdualmesh->GetElemDualInfo(i, qecoord, qnormals, qelengths, qareas);  
  
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
  locrhs[3] = 0.0;
  locrhs[4] = 0.0;
  locrhs[5] = 0.0;

  Function *ff = data->GetForceFunction();  
  /*
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
	  bcoord[0] = xx;
	  bcoord[1] = yy;
	  locrhs[l/2] += tpw[k][2] * ff->Eval(bcoord) * detJ;
	} 
    }
  */

  for(int l=0; l<6; l++)
    {
      if(l==0)
	{
	  x[0] = coord[0][0];
	  y[0] = coord[1][0];
	  x[1] = qecoord[0][0];
	  y[1] = qecoord[0][1];
	  x[2] = qecoord[2][0];
	  y[2] = qecoord[2][1];
	}
      else if(l==1)
	{
	  x[0] = qecoord[3][0];
	  y[0] = qecoord[3][1];
	  x[1] = qecoord[2][0];
	  y[1] = qecoord[2][1];
	  x[2] = qecoord[0][0];
	  y[2] = qecoord[0][1];
	}
      else if(l==2)
	{
	  x[0] = coord[0][1];
	  y[0] = coord[1][1];
	  x[1] = qecoord[4][0];
	  y[1] = qecoord[4][1];
	  x[2] = qecoord[6][0];
	  y[2] = qecoord[6][1];
	}
      else if(l==3)
	{
	  x[0] = qecoord[7][0];
	  y[0] = qecoord[7][1];
	  x[1] = qecoord[6][0];
	  y[1] = qecoord[6][1];
	  x[2] = qecoord[4][0];
	  y[2] = qecoord[4][1];
	}
      else if(l==4)
	{
	  x[0] = coord[0][2];
	  y[0] = coord[1][2];
	  x[1] = qecoord[8][0];
	  y[1] = qecoord[8][1];
	  x[2] = qecoord[10][0];
	  y[2] = qecoord[10][1];
	}
      else
	{
	  x[0] = qecoord[11][0];
	  y[0] = qecoord[11][1];
	  x[1] = qecoord[10][0];
	  y[1] = qecoord[10][1];
	  x[2] = qecoord[8][0];
	  y[2] = qecoord[8][1];
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
  
  for(int l=0; l<18; l++)
    {
      if(l==0)
	{
	  x[0] = 0.5 * (coord[0][0] + coord[0][1]);
	  y[0] = 0.5 * (coord[1][0] + coord[1][1]);
	  x[1] = qecoord[1][0];
	  y[1] = qecoord[1][1];
	  x[2] = qecoord[0][0];
	  y[2] = qecoord[0][1];
	}
      else if(l==1)
	{
	  x[0] = qecoord[3][0];
	  y[0] = qecoord[3][1];
	  x[1] = qecoord[0][0];
	  y[1] = qecoord[0][1];
	  x[2] = qecoord[1][0];
	  y[2] = qecoord[1][1];
	}
      else if(l==2)
	{
	  x[0] = 0.5 * (coord[0][0] + coord[0][1]);
	  y[0] = 0.5 * (coord[1][0] + coord[1][1]);
	  x[1] = qecoord[6][0];
	  y[1] = qecoord[6][1];
	  x[2] = qecoord[5][0];
	  y[2] = qecoord[5][1];
	}
      else if(l==3)
	{
	  x[0] = qecoord[7][0];
	  y[0] = qecoord[7][1];
	  x[1] = qecoord[5][0];
	  y[1] = qecoord[5][1];
	  x[2] = qecoord[6][0];
	  y[2] = qecoord[6][1];
	}
      else if(l==4)
	{
	  x[0] = 0.5 * (coord[0][0] + coord[0][1]);
	  y[0] = 0.5 * (coord[1][0] + coord[1][1]);
	  x[1] = qecoord[5][0];
	  y[1] = qecoord[5][1];
	  x[2] = qecoord[1][0];
	  y[2] = qecoord[1][1];
	}
      else if(l==5)
	{
	  x[0] = qecoord[12][0];
	  y[0] = qecoord[12][1];
	  x[1] = qecoord[1][0];
	  y[1] = qecoord[1][1];
	  x[2] = qecoord[5][0];
	  y[2] = qecoord[5][1];
	}
      else if(l==6)
	{
	  x[0] = 0.5 * (coord[0][2] + coord[0][1]);
	  y[0] = 0.5 * (coord[1][2] + coord[1][1]);
	  x[1] = qecoord[5][0];
	  y[1] = qecoord[5][1];
	  x[2] = qecoord[4][0];
	  y[2] = qecoord[4][1];
	}
      else if(l==7)
	{
	  x[0] = qecoord[7][0];
	  y[0] = qecoord[7][1];
	  x[1] = qecoord[4][0];
	  y[1] = qecoord[4][1];
	  x[2] = qecoord[5][0];
	  y[2] = qecoord[5][1];
	}
      else if(l==8)
	{
	  x[0] = 0.5 * (coord[0][2] + coord[0][1]);
	  y[0] = 0.5 * (coord[1][2] + coord[1][1]);
	  x[1] = qecoord[10][0];
	  y[1] = qecoord[10][1];
	  x[2] = qecoord[9][0];
	  y[2] = qecoord[9][1];
	}
      else if(l==9)
	{
	  x[0] = qecoord[11][0];
	  y[0] = qecoord[11][1];
	  x[1] = qecoord[9][0];
	  y[1] = qecoord[9][1];
	  x[2] = qecoord[10][0];
	  y[2] = qecoord[10][1];
	}
      else if(l==10)
	{
	  x[0] = 0.5 * (coord[0][2] + coord[0][1]);
	  y[0] = 0.5 * (coord[1][2] + coord[1][1]);
	  x[1] = qecoord[9][0];
	  y[1] = qecoord[9][1];
	  x[2] = qecoord[5][0];
	  y[2] = qecoord[5][1];
	}
      else if(l==11)
	{
	  x[0] = qecoord[12][0];
	  y[0] = qecoord[12][1];
	  x[1] = qecoord[5][0];
	  y[1] = qecoord[5][1];
	  x[2] = qecoord[9][0];
	  y[2] = qecoord[9][1];
	}
      else if(l==12)
	{
	  x[0] = 0.5 * (coord[0][2] + coord[0][0]);
	  y[0] = 0.5 * (coord[1][2] + coord[1][0]);
	  x[1] = qecoord[2][0];
	  y[1] = qecoord[2][1];
	  x[2] = qecoord[1][0];
	  y[2] = qecoord[1][1];
	}
      else if(l==13)
	{
	  x[0] = qecoord[3][0];
	  y[0] = qecoord[3][1];
	  x[1] = qecoord[1][0];
	  y[1] = qecoord[1][1];
	  x[2] = qecoord[2][0];
	  y[2] = qecoord[2][1];
	}
      else if(l==14)
	{
	  x[0] = 0.5 * (coord[0][2] + coord[0][0]);
	  y[0] = 0.5 * (coord[1][2] + coord[1][0]);
	  x[1] = qecoord[9][0];
	  y[1] = qecoord[9][1];
	  x[2] = qecoord[8][0];
	  y[2] = qecoord[8][1];
	}
      else if(l==15)
	{
	  x[0] = qecoord[11][0];
	  y[0] = qecoord[11][1];
	  x[1] = qecoord[8][0];
	  y[1] = qecoord[8][1];
	  x[2] = qecoord[9][0];
	  y[2] = qecoord[9][1];
	}
      else if(l==16)
	{
	  x[0] = 0.5 * (coord[0][2] + coord[0][0]);
	  y[0] = 0.5 * (coord[1][2] + coord[1][0]);
	  x[1] = qecoord[1][0];
	  y[1] = qecoord[1][1];
	  x[2] = qecoord[9][0];
	  y[2] = qecoord[9][1];
	}
      else
	{
	  x[0] = qecoord[12][0];
	  y[0] = qecoord[12][1];
	  x[1] = qecoord[9][0];
	  y[1] = qecoord[9][1];
	  x[2] = qecoord[1][0];
	  y[2] = qecoord[1][1];
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
	  locrhs[3 + l/6] += tpw[k][2] * ff->Eval(bcoord) * detJ;
	} 
    }
  
  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
  for (int k=0; k<tpw.Size(); k++) { delete []tpw[k]; }

  for (int k=0; k<ecoord.Size(); k++) { delete []ecoord[k]; }
  for (int k=0; k<normals.Size(); k++) { delete []normals[k]; }
  for (int k=0; k<qecoord.Size(); k++) { delete []qecoord[k]; }
  for (int k=0; k<qnormals.Size(); k++) { delete []qnormals[k]; }
}

//==============================================================================

void QuadraticFVEM2D::SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
			       const Array<double *> &Dval)
{
  Element *el;  
  for(int j=0; j<bdry_is_dirichlet.Size(); j++)
    {
      for (int i=0; i<GetNBE(); i++)
	{
	  el = GetBdrElement(i);
	  int attr = el->GetAttribute();
	  int nv = el->GetNVertices() + 1;
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

void QuadraticFVEM2D::ProjectAFunction(double (*func)(Array<double> &), Vector &pval)
{// not right
  Array<double> coord(2);
  for(int i=0; i<mesh->GetNV(); i++)
    {
      coord[0] = mesh->GetVertex(i, 0);
      coord[1] = mesh->GetVertex(i, 1);      
      pval(i) = func(coord);
    }
}
