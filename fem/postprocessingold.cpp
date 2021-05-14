#include "fem_header.h"
#include "../general/general_header.h"

Postprocessing::Postprocessing(DualMesh *_dualmesh, FEM *_fem, Data *_data)
{
  dualmesh = _dualmesh;
  mesh = dualmesh->GetMesh();
  fem = _fem;
  data = _data;
  if(mesh->GetMType()==Element::TRIANGLE)
    {
      if( dualmesh->GetDualMeshOrder()==1 )
	NumLocalDOF = 3;
      else if( dualmesh->GetDualMeshOrder()==2 )
	NumLocalDOF = 6;
      else
	NumLocalDOF = 10;
    }
  else
    {
      if( dualmesh->GetDualMeshOrder()==1 )
	NumLocalDOF = 4;
      else if( dualmesh->GetDualMeshOrder()==2 )
	NumLocalDOF = 9;
      else
	NumLocalDOF = 16;
    }
}

//=========================================================

void Postprocessing::ComputeConservativeFlux(Vector &sol, Array<double*> &flux, Array<double *> &ppsol,
					     Array<double *> &ppflux, Array<double> &ppbdrflux)
{
  int n = NumLocalDOF;
  int NE = mesh->GetNE();

  for(int j=0; j<ppbdrflux.Size(); j++) ppbdrflux[j] = 0.0;
  
  for(int j=0; j<n; j++)
    for(int i=0; i<NE; i++)
      ppsol[j][i] = 0.0;
	  
  for(int j=0; j<flux.Size(); j++)
    {
      for(int i=0; i<NE; i++)
	{
	  flux[j][i] = 0.0;
	  ppflux[j][i] = 0.0;
	}
    }
  
  if(dualmesh->GetDualMeshOrder()==1)
    {
      if(mesh->GetMType()==Element::TRIANGLE)
	ComputeLinearConservativeTriangularFlux(sol, flux, ppsol, ppflux, ppbdrflux);
      else
	ComputeLinearConservativeRectangularFlux(sol, flux, ppsol, ppflux, ppbdrflux);
    }
  else if(dualmesh->GetDualMeshOrder()==2)
    {
      if(mesh->GetMType()==Element::TRIANGLE)
	ComputeQuadraticConservativeTriangularFlux(sol, flux, ppsol, ppflux, ppbdrflux);
      else
	ComputeQuadraticConservativeRectangularFlux(sol, flux, ppsol, ppflux, ppbdrflux);
    }
  else
    {
      if(mesh->GetMType()==Element::TRIANGLE)
	ComputeCubicConservativeTriangularFlux(sol, flux, ppsol, ppflux, ppbdrflux);
      else
	ComputeCubicConservativeRectangularFlux(sol, flux, ppsol, ppflux, ppbdrflux);
    }
}

//===========================================================

void Postprocessing::ComputeLinearConservativeTriangularFlux(Vector &sol, Array<double*> &flux, Array<double*> &ppsol,
							     Array<double*> &ppflux, Array<double> &ppbdrflux)
{
  
  mesh->GetEdgeToElementTable();
    
  Array<double *> normals(3);
  Array<double> elengths(3);
  Array<double *> ecoord(4);
  Array<double> areas(3);  
  for(int k=0; k<ecoord.Size(); k++) ecoord[k] = new double[2];
  for(int k=0; k<normals.Size(); k++) normals[k] = new double[2];

  int n = 4;
  Array<double *> pw(n);
  for(int i=0; i<n; i++) pw[i] = new double[2];
  GetStandLineQuadPW(n, pw);

  int np = ConvertN2PN(n+1);
  Array<double *> tpw(np);
  for(int i=0; i<np; i++) tpw[i] = new double[3];
  GetStandTriQuadPW(n+1, tpw);

  Array<double *> sflux(mesh->GetNE());
  for(int i=0; i<sflux.Size(); i++) sflux[i] = new double[6];
  for(int j=0; j<sflux.Size(); j++)
    for(int i=0; i<6; i++)
      sflux[j][i] = 0.0;

  Array<double *> tflux(mesh->GetNE());
  for(int i=0; i<tflux.Size(); i++) tflux[i] = new double[6];
  for(int j=0; j<tflux.Size(); j++)
    for(int i=0; i<6; i++)
      tflux[j][i] = 0.0;

  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      mesh->GetElementVertices(i, ind);

      Array<double *> coord(mesh->GetDim());
      for (int k=0; k<coord.Size(); k++) coord[k] = new double[3];
      
      mesh->GetElementVerticesCoord(i, coord);
      double det = (coord[0][1] - coord[0][0]) * (coord[1][2] - coord[1][0])
	- (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);
      
      double dp0x, dp0y, dp1x, dp1y, dp2x, dp2y;
      dp0x = (coord[1][2]-coord[1][1])/det;
      dp0y = (coord[0][1]-coord[0][2])/det;
      dp1x = (coord[1][2]-coord[1][0])/det;
      dp1y = (coord[0][0]-coord[0][2])/det;
      dp2x = (coord[1][0]-coord[1][1])/det;
      dp2y = (coord[0][1]-coord[0][0])/det;
      
      dualmesh->GetElemDualInfo(i, ecoord, normals, elengths, areas); 

      Vector locsol(3);
      for(int k=0; k<3; k++) locsol(k) = sol(ind[k]);

      Array<double *> btem(6);
      for(int k=0; k<6; k++) btem[k] = new double[3];
      for(int k=0; k<6; k++)
	for(int j=0; j<3; j++)
	  btem[k][j] = 0.0;
      
      double xa, xb, ya, yb, tx, ty, ttx, tty;
      Function *kk = data->GetEllipticFunction();
      
      for(int k=0; k<6; k++)
	{
	  double len = 0.0;
	  double nn[2];
	  if(k==0)
	    {
	      xa = coord[0][0];
	      xb = ecoord[0][0];
	      ya = coord[1][0];
	      yb = ecoord[0][1];
	    }
	  else if(k==1)
	    {
	      xa = ecoord[0][0];
	      xb = coord[0][1];
	      ya = ecoord[0][1];
	      yb = coord[1][1];
	    }
	  else if(k==2)
	    {
	      xa = coord[0][1];
	      xb = ecoord[1][0];
	      ya = coord[1][1];
	      yb = ecoord[1][1];
	    }
	  else if(k==3)
	    {
	      xa = ecoord[1][0];
	      xb = coord[0][2];
	      ya = ecoord[1][1];
	      yb = coord[1][2];
	    }
	  else if(k==4)
	    {
	      xa = coord[0][2];
	      xb = ecoord[2][0];
	      ya = coord[1][2];
	      yb = ecoord[2][1];
	    }
	  else
	    {
	      xa = ecoord[2][0];
	      xb = coord[0][0];
	      ya = ecoord[2][1];
	      yb = coord[1][0];
	    }

	  len = sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya) );
	  nn[0] = (yb - ya)/len;
	  nn[1] = (xa - xb)/len;

	  for(int l=0; l<n; l++)
	    {
	      tx = 0.5 * (xb - xa) * pw[l][0] + 0.5 * (xa + xb);
	      ty = 0.5 * (yb - ya) * pw[l][0] + 0.5 * (ya + yb);

	      double bcoord[2];
	      bcoord[0] = tx;
	      bcoord[1] = ty;
	      btem[k][0] += -pw[l][1] * kk->Eval(bcoord) * ( dp0x*nn[0] + dp0y*nn[1] ) * 0.5 * len;
	      btem[k][1] += pw[l][1] * kk->Eval(bcoord) * ( dp1x*nn[0] + dp1y*nn[1] ) * 0.5 * len;
	      btem[k][2] += pw[l][1] * kk->Eval(bcoord) * ( dp2x*nn[0] + dp2y*nn[1] ) * 0.5 * len; 
	    }
	  for(int l=0; l<3; l++) sflux[i][k] += btem[k][l] * locsol(l);	  
	}

      for(int k=0; k<6; k++)
	for(int j=0; j<3; j++)
	  btem[k][j] = 0.0;
      
      for(int k=0; k<6; k++)
	{
	  double len = 0.0;
	  double nn[2];
	  if(k<2)
	    {
	      xa = coord[0][0];
	      xb = coord[0][1];
	      ya = coord[1][0];
	      yb = coord[1][1];
	    }
	  else if(k<4)
	    {
	      xa = coord[0][1];
	      xb = coord[0][2];
	      ya = coord[1][1];
	      yb = coord[1][2];
	    }
	  else
	    {
	      xa = coord[0][2];
	      xb = coord[0][0];
	      ya = coord[1][2];
	      yb = coord[1][0];
	    }

	  len = sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya) );
	  nn[0] = (yb - ya)/len;
	  nn[1] = (xa - xb)/len;

	  for(int l=0; l<n; l++)
	    {
	      tx = 0.5 * (xb - xa) * pw[l][0] + 0.5 * (xa + xb);
	      ty = 0.5 * (yb - ya) * pw[l][0] + 0.5 * (ya + yb);
	      ttx = dp1x * (tx - coord[0][0]) + dp1y * (ty - coord[1][0]);
	      tty = dp2x * (tx - coord[0][0]) + dp2y * (ty - coord[1][0]);

	      double bcoord[2];
	      bcoord[0] = tx;
	      bcoord[1] = ty;
	      double temprod0 = pw[l][1] * kk->Eval(bcoord) * (dp0x*nn[0] + dp0y*nn[1]) * 0.5 * len;
	      double temprod1 = pw[l][1] * kk->Eval(bcoord) * (dp1x*nn[0] + dp1y*nn[1]) * 0.5 * len;
	      double temprod2 = pw[l][1] * kk->Eval(bcoord) * (dp2x*nn[0] + dp2y*nn[1]) * 0.5 * len;
	      
	      if(k==0)
		{
		  btem[k][0] += -(1.0 - ttx - tty) * temprod0;
		  btem[k][1] += (1.0 - ttx - tty) * temprod1;
		  btem[k][2] += (1.0 - ttx - tty) * temprod2;
		}
	      else if(k==1)
		{
		  btem[k][0] += -ttx * temprod0;
		  btem[k][1] += ttx * temprod1;
		  btem[k][2] += ttx * temprod2;
		}
	      else if(k==2)
		{
		  btem[k][0] += -ttx * temprod0;
		  btem[k][1] += ttx * temprod1;
		  btem[k][2] += ttx * temprod2;
		}
	      else if(k==3)
		{
		  btem[k][0] += -tty * temprod0;
		  btem[k][1] += tty * temprod1;
		  btem[k][2] += tty * temprod2;
		}
	      else if(k==4)
		{
		  btem[k][0] += -tty * temprod0;
		  btem[k][1] += tty * temprod1;
		  btem[k][2] += tty * temprod2;
		}
	      else
		{
		  btem[k][0] += -(1.0 - ttx - tty) * temprod0;
		  btem[k][1] += (1.0 - ttx - tty) * temprod1;
		  btem[k][2] += (1.0 - ttx - tty) * temprod2;
		}
	    }
	  for(int l=0; l<3; l++) tflux[i][k] += btem[k][l] * locsol(l);	  
	}
      for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
      for (int k=0; k<btem.Size(); k++) { delete []btem[k]; }
    }

  Vector cverr(mesh->GetNV());
  cverr = 0.0;

  Vector cverruh(mesh->GetNV());
  cverruh = 0.0;

  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      mesh->GetElementVertices(i, ind);
      
      Array<double *> coord(mesh->GetDim());
      for (int k=0; k<coord.Size(); k++) coord[k] = new double[3];
      
      mesh->GetElementVerticesCoord(i, coord);
      double det = (coord[0][1] - coord[0][0]) * (coord[1][2] 
	    - coord[1][0]) - (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);

      double dp0x, dp0y, dp1x, dp1y, dp2x, dp2y;
      dp0x = (coord[1][2]-coord[1][1])/det;
      dp0y = (coord[0][1]-coord[0][2])/det;
      dp1x = (coord[1][2]-coord[1][0])/det;
      dp1y = (coord[0][0]-coord[0][2])/det;
      dp2x = (coord[1][0]-coord[1][1])/det;
      dp2y = (coord[0][1]-coord[0][0])/det;

      //------------------------------------
      dualmesh->GetElemDualInfo(i, ecoord, normals, elengths, areas); 

      Array<double *> tem(3);
      for(int k=0; k<3; k++) tem[k] = new double[3];
      for(int k=0; k<3; k++)
	for(int j=0; j<3; j++)
	  tem[k][j] = 0.0;

      double xa, xb, ya, yb, tx, ty;
      Function *kk = data->GetEllipticFunction();
      Function *ff = data->GetForceFunction();
      
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
	      
	      double bcoord[2];
	      bcoord[0] = tx;
	      bcoord[1] = ty;
	      tem[k][0] += -pw[l][1] * kk->Eval(bcoord) * ( dp0x*normals[k][0] + dp0y*normals[k][1] ) * 0.5 * elengths[k];
	      tem[k][1] += pw[l][1] * kk->Eval(bcoord) * ( dp1x*normals[k][0] + dp1y*normals[k][1] ) * 0.5 * elengths[k];
	      tem[k][2] += pw[l][1] * kk->Eval(bcoord) * ( dp2x*normals[k][0] + dp2y*normals[k][1] ) * 0.5 * elengths[k];
	    }
	}

      SparseMatrix *locmat = new SparseMatrix(3,3);
      locmat->Elem(0,0) = tem[2][0] - tem[0][0];
      locmat->Elem(0,1) = tem[2][1] - tem[0][1];
      locmat->Elem(0,2) = tem[2][2] - tem[0][2];
      
      locmat->Elem(1,0) = tem[0][0] - tem[1][0];
      locmat->Elem(1,1) = tem[0][1] - tem[1][1];
      locmat->Elem(1,2) = tem[0][2] - tem[1][2];
      
      locmat->Elem(2,0) = tem[1][0] - tem[2][0];
      locmat->Elem(2,1) = tem[1][1] - tem[2][1];
      locmat->Elem(2,2) = tem[1][2] - tem[2][2];
      
      locmat->Finalize();

      Vector Q(3);
      Q = 0.0;
      Array<double *> femat(3);
      Vector locsol(3);
      for(int k=0; k<3; k++)
	{
	  femat[k] = new double[3];
	  locsol(k) = sol(ind[k]);
	}
      
      fem->GetEllipticLocalSystem(i, ind, femat);
 
      for(int k=0; k<3; k++)
	for(int j=0; j<3; j++)
	  Q(k) += femat[k][j] * locsol(j);

      Array<double> F(3);
      fem->GetForceLocalSystem(i, ind, F);

      Vector f(3);
      f = 0.0;

      Array<double> x(3);
      Array<double> y(3);
      double t1, t2, t3, t4, detJ;

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
	      f(l/2) += tpw[k][2] * ff->Eval(bcoord) * detJ;
	    } 
	}

      Vector RHS(3);
      for(int k=0; k<3; k++) RHS(k) = Q(k) - F[k] + f(k);
      
      Array<int> edges;
      mesh->GetElementEdges(i, edges);
      if(i%2==0)
	{
	  Array<int> eind;
	  int k = edges[0];
	  mesh->GetEdgeElements(k, eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[1]][3] - tflux[i][0] + tflux[eind[1]][3] );
		  RHS(1) += 0.5 * ( sflux[i][1] - sflux[eind[1]][2] - tflux[i][1] + tflux[eind[1]][2] );
		}
	      else
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[0]][3] - tflux[i][0] + tflux[eind[0]][3] );
		  RHS(1) += 0.5 * ( sflux[i][1] - sflux[eind[0]][2] - tflux[i][1] + tflux[eind[0]][2] );
		}
	    }
	  else
	    {
	      RHS(0) += -tflux[i][0] + sflux[i][0];
	      RHS(1) += -tflux[i][1] + sflux[i][1];
	    }

	  k = edges[1];
	  mesh->GetEdgeElements(k, eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(1) += 0.5 * ( sflux[i][2] - sflux[eind[1]][5] - tflux[i][2] + tflux[eind[1]][5] );
		  RHS(2) += 0.5 * ( sflux[i][3] - sflux[eind[1]][4] - tflux[i][3] + tflux[eind[1]][4] );
		}
	      else
		{
		  RHS(1) += 0.5 * ( sflux[i][2] - sflux[eind[0]][5] - tflux[i][2] + tflux[eind[0]][5] );
		  RHS(2) += 0.5 * ( sflux[i][3] - sflux[eind[0]][4] - tflux[i][3] + tflux[eind[0]][4] );
		}
	    }
	  else
	    {
	      RHS(1) += -tflux[i][2] + sflux[i][2];
	      RHS(2) += -tflux[i][3] + sflux[i][3];
	    }

	 k = edges[2];
	 mesh->GetEdgeElements(k, eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(0) += 0.5 * ( sflux[i][5] - sflux[eind[1]][0] - tflux[i][5] + tflux[eind[1]][0] );
		  RHS(2) += 0.5 * ( sflux[i][4] - sflux[eind[1]][1] - tflux[i][4] + tflux[eind[1]][1] );
		}
	      else
		{
		  RHS(0) += 0.5 * ( sflux[i][5] - sflux[eind[0]][0] - tflux[i][5] + tflux[eind[0]][0] );
		  RHS(2) += 0.5 * ( sflux[i][4] - sflux[eind[0]][1] - tflux[i][4] + tflux[eind[0]][1] );
		}
	    }
	  else
	    {
	      RHS(0) += -tflux[i][5] + sflux[i][5];
	      RHS(2) += -tflux[i][4] + sflux[i][4];
	    }
	}
      else
	{
	  Array<int> eind;
	  int k = edges[0];
	  mesh->GetEdgeElements(k, eind);
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[1]][5] - tflux[i][0] + tflux[eind[1]][5] );
		  RHS(1) += 0.5 * ( sflux[i][1] - sflux[eind[1]][4] - tflux[i][1] + tflux[eind[1]][4] );
		}
	      else
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[0]][5] - tflux[i][0] + tflux[eind[0]][5] );
		  RHS(1) += 0.5 * ( sflux[i][1] - sflux[eind[0]][4] - tflux[i][1] + tflux[eind[0]][4] );
		}
	    }
	  else
	    {
	      RHS(0) += -tflux[i][0] + sflux[i][0];
	      RHS(1) += -tflux[i][1] + sflux[i][1];
	    }

	  k = edges[1];
	  mesh->GetEdgeElements(k, eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(1) += 0.5 * ( sflux[i][2] - sflux[eind[1]][1] - tflux[i][2] + tflux[eind[1]][1] );
		  RHS(2) += 0.5 * ( sflux[i][3] - sflux[eind[1]][0] - tflux[i][3] + tflux[eind[1]][0] );
		}
	      else
		{
		  RHS(1) += 0.5 * ( sflux[i][2] - sflux[eind[0]][1] - tflux[i][2] + tflux[eind[0]][1] );
		  RHS(2) += 0.5 * ( sflux[i][3] - sflux[eind[0]][0] - tflux[i][3] + tflux[eind[0]][0] );
		}
	    }
	  else
	    {
	      RHS(1) += -tflux[i][2] + sflux[i][2];
	      RHS(2) += -tflux[i][3] + sflux[i][3];
	    }

	 k = edges[2];
	 mesh->GetEdgeElements(k, eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(0) += 0.5 * ( sflux[i][5] - sflux[eind[1]][2] - tflux[i][5] + tflux[eind[1]][2] );
		  RHS(2) += 0.5 * ( sflux[i][4] - sflux[eind[1]][3] - tflux[i][4] + tflux[eind[1]][3] );
		}
	      else
		{
		  RHS(0) += 0.5 * ( sflux[i][5] - sflux[eind[0]][2] - tflux[i][5] + tflux[eind[0]][2] );
		  RHS(2) += 0.5 * ( sflux[i][4] - sflux[eind[0]][3] - tflux[i][4] + tflux[eind[0]][3] );
		}
	    }
	  else
	    {
	      RHS(0) += -tflux[i][5] + sflux[i][5];
	      RHS(2) += -tflux[i][4] + sflux[i][4];
	    }
	}
            
      Vector s(3);
      s = 0.0;
      Vector rhs(3);
      rhs = RHS;
            cout<<RHS(0)+RHS(1)+RHS(2)<<endl;
      RectangularMatrix localconsmat(3);
      
      for(int j=0; j<3; j++)
	for(int k=0; k<3; k++)
	  localconsmat(j,k) = locmat->Elem(j,k);

      for(int j=1; j<3; j++) localconsmat(0, j) = 0.0;
      localconsmat(0,0) = 1.0;
      RHS(0) = locsol(0); 

      DenseMatrixInverse invmat(localconsmat);
      invmat.Mult(RHS, s);

      Vector sss(3);
      sss = 0.0;
      locmat->Mult(s, sss);
      for(int k=0; k<3; k++) cverr(ind[k]) += sss(k) - f(k); 

      locmat->Mult(locsol, sss);
      for(int k=0; k<3; k++) cverruh(ind[k]) += sss(k) - f(k); 

      for(int k=0; k<3; k++) { ppsol[k][i] = s(k); }
      for(int k=0; k<3; k++)
	{
	  for(int j=0; j<3; j++)
	    {
	      ppflux[k][i] += -tem[k][j] * s(j);
	      flux[k][i] += -tem[k][j] * locsol(j);
	    }
	}
      
      for(int k=0; k<femat.Size(); k++) delete []femat[k];
      for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
      for (int k=0; k<tem.Size(); k++) { delete []tem[k]; }
      delete locmat;
    }      
 
  for(int i=0; i<mesh->GetNBE(); i++)
    {
      Array<int> bdrind(2);
      mesh->GetBdrElementVertices(i, bdrind);
      if( fabs(cverr(bdrind[0])) > 1.0e-15 )
	{
	  ppbdrflux[bdrind[0]] = -cverr(bdrind[0]);
	  cverr(bdrind[0]) = 0.0;
	  cverruh(bdrind[0]) = 0.0;	  
	}
      if( fabs(cverr(bdrind[1])) > 1.0e-15 )
	{
	  ppbdrflux[bdrind[1]] = -cverr(bdrind[1]);
	  cverr(bdrind[1]) = 0.0;
	  cverruh(bdrind[1]) = 0.0;
	}
    }

  ofstream fileout("llce.out");
  for(int i=0; i<dualmesh->GetDualMeshNumDOF(); i++)
    {
      fileout<<i<<"\t"<<cverruh(i)<< "\t" << cverr(i) << endl;
      cout<<setprecision(6)<<i<<"\t"<<cverruh(i)<< "\t" <<cverr(i)<<endl;
    }
  fileout.close(); 
  
  for(int i=0; i<tflux.Size(); i++) delete []tflux[i];
  for(int i=0; i<sflux.Size(); i++) delete []sflux[i];
  for(int i=0; i<pw.Size(); i++) delete []pw[i];
  for(int i=0; i<tpw.Size(); i++) delete []tpw[i];
  for (int k=0; k<ecoord.Size(); k++) { delete []ecoord[k]; }
  for (int k=0; k<normals.Size(); k++) { delete []normals[k]; }
}

//===========================================================

void Postprocessing::ComputeQuadraticConservativeTriangularFlux(Vector &sol, Array<double*> &flux, Array<double*> &ppsol,
								Array<double*> &ppflux, Array<double> &ppbdrflux)
{
  mesh->GetEdgeToElementTable();
    
  Array<double *> nmls(3);
  Array<double> lens(3);
  Array<double *> cod(13);
  Array<double> areas(1);  
  for(int k=0; k<cod.Size(); k++) cod[k] = new double[2];
  for(int k=0; k<nmls.Size(); k++) nmls[k] = new double[2];

  int n = 4;
  Array<double *> pw(n);
  for(int i=0; i<n; i++) pw[i] = new double[2];
  GetStandLineQuadPW(n, pw);

  int np = ConvertN2PN(5);
  Array<double *> tpw(np);
  for(int i=0; i<np; i++) tpw[i] = new double[3];
  GetStandTriQuadPW(5, tpw);

  Array<double *> sflux(mesh->GetNE());
  for(int i=0; i<sflux.Size(); i++) sflux[i] = new double[9];
  for(int j=0; j<sflux.Size(); j++)
    for(int i=0; i<9; i++)
      sflux[j][i] = 0.0;

  Array<double *> tflux(mesh->GetNE());
  for(int i=0; i<tflux.Size(); i++) tflux[i] = new double[9];
  for(int j=0; j<tflux.Size(); j++)
    for(int i=0; i<9; i++)
      tflux[j][i] = 0.0;

  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      fem->getElementDOF(i, ind);
      Array<int> edges;
      mesh->GetElementEdges(i, edges);

      Array<double *> coord(mesh->GetDim());
      for (int k=0; k<coord.Size(); k++) coord[k] = new double[3];
      
      mesh->GetElementVerticesCoord(i, coord);

      double t1 = coord[0][1] - coord[0][0];
      double t2 = coord[0][2] - coord[0][0];
      double t3 = coord[1][1] - coord[1][0];
      double t4 = coord[1][2] - coord[1][0];
      double detJ = t1 * t4 - t2 * t3;
      
      // Jacobian Inverse Transpose * detJ
      Array<double *> jit(2);
      jit[0] = new double[2];
      jit[1] = new double[2];
      jit[0][0] = coord[1][2] - coord[1][0];
      jit[0][1] = coord[1][0] - coord[1][1];
      jit[1][0] = coord[0][0] - coord[0][2];
      jit[1][1] = coord[0][1] - coord[0][0];

      for(int k=0; k<2; k++)
	for(int j=0; j<2; j++)
	  jit[k][j] /= detJ;

      dualmesh->GetElemDualInfo(i, cod, nmls, lens, areas); 

      Vector locsol(6);
      for(int k=0; k<6; k++) locsol(k) = sol(ind[k]);

      Array<double *> btem(9);
      for(int k=0; k<9; k++) btem[k] = new double[6];
      for(int k=0; k<9; k++)
	for(int j=0; j<6; j++)
	  btem[k][j] = 0.0;
      
      double xa, xb, ya, yb, tx, ty, ttx, tty;

      BasisFunctions *phi = new BasisFunctions(2);
      Function *kk = data->GetEllipticFunction();
      for(int k=0; k<9; k++)
	{
	  double len = 0.0;
	  double nn[2];
	  if(k==1)
	    {
	      xa = cod[0][0];
	      xb = cod[6][0];
	      ya = cod[0][1];
	      yb = cod[6][1];
	    }
	  else if(k==4)
	    {
	      xa = cod[4][0];
	      xb = cod[10][0];
	      ya = cod[4][1];
	      yb = cod[10][1];
	    }
	  else if(k==7)
	    {
	      xa = cod[8][0];
	      xb = cod[2][0];
	      ya = cod[8][1];
	      yb = cod[2][1];
	    }
	  else if(k==0)
	    {
	      xa = coord[0][0];
	      xb = cod[0][0];
	      ya = coord[1][0];
	      yb = cod[0][1];
	    }
	  else if(k==2)
	    {
	      xa = cod[6][0];
	      xb = coord[0][1];
	      ya = cod[6][1];
	      yb = coord[1][1];
	    }
	  else if(k==3)
	    {
	      xa = coord[0][1];
	      xb = cod[4][0];
	      ya = coord[1][1];
	      yb = cod[4][1];
	    }
	  else if(k==5)
	    {
	      xa = cod[10][0];
	      xb = coord[0][2];
	      ya = cod[10][1];
	      yb = coord[1][2];
	    }
	  else if(k==6)
	    {
	      xa = coord[0][2];
	      xb = cod[8][0];
	      ya = coord[1][2];
	      yb = cod[8][1];
	    }
	  else
	    {
	      xa = cod[2][0];
	      xb = coord[0][0];
	      ya = cod[2][1];
	      yb = coord[1][0];
	    }

	  len = sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya) );
	  nn[0] = (yb - ya)/len;
	  nn[1] = (xa - xb)/len;

	  for(int l=0; l<n; l++)
	    {
	      tx = 0.5 * (xb - xa) * pw[l][0] + 0.5 * (xa + xb);
	      ty = 0.5 * (yb - ya) * pw[l][0] + 0.5 * (ya + yb);
	      ttx = jit[0][0] * (tx - coord[0][0]) + jit[1][0] * (ty - coord[1][0]);
	      tty = jit[0][1] * (tx - coord[0][0]) + jit[1][1] * (ty - coord[1][0]);

	      double bcoord[2];
	      bcoord[0] = tx;
	      bcoord[1] = ty;
	      for(int j=0; j<6; j++)
		{
		  double drx = jit[0][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[0][1] * phi->GradBF2D(j, 1, ttx, tty);
		  double dry = jit[1][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[1][1] * phi->GradBF2D(j, 1, ttx, tty);
		  btem[k][j] += pw[l][1] * kk->Eval(bcoord) * (drx*nn[0] + dry*nn[1]) * 0.5 * len; 
		}
	      //sflux[i][k] += pw[l][1] * kk-Eval(bcoord) * ( uux(tx, ty)*nn[0] + uuy(tx, ty)*nn[1] ) * 0.5 * len;
	    }
	  for(int l=0; l<6; l++) sflux[i][k] += btem[k][l] * locsol(l);	  
	}

      for(int k=0; k<9; k++)
	for(int j=0; j<6; j++)
	  btem[k][j] = 0.0;
      
      for(int k=0; k<9; k++)
	{
	  double len = 0.0;
	  double nn[2];
	  if(k<3)
	    {
	      xa = coord[0][0];
	      xb = coord[0][1];
	      ya = coord[1][0];
	      yb = coord[1][1];
	    }
	  else if(k<6)
	    {
	      xa = coord[0][1];
	      xb = coord[0][2];
	      ya = coord[1][1];
	      yb = coord[1][2];
	    }
	  else
	    {
	      xa = coord[0][2];
	      xb = coord[0][0];
	      ya = coord[1][2];
	      yb = coord[1][0];
	    }

	  len = sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya) );
	  nn[0] = (yb - ya)/len;
	  nn[1] = (xa - xb)/len;

	  for(int l=0; l<n; l++)
	    {
	      tx = 0.5 * (xb - xa) * pw[l][0] + 0.5 * (xa + xb);
	      ty = 0.5 * (yb - ya) * pw[l][0] + 0.5 * (ya + yb);
	      ttx = jit[0][0] * (tx - coord[0][0]) + jit[1][0] * (ty - coord[1][0]);
	      tty = jit[0][1] * (tx - coord[0][0]) + jit[1][1] * (ty - coord[1][0]);
	      
	      double bcoord[2];
	      bcoord[0] = tx;
	      bcoord[1] = ty;
	      for(int j=0; j<6; j++)
		{
		  double drx = jit[0][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[0][1] * phi->GradBF2D(j, 1, ttx, tty);
		  double dry = jit[1][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[1][1] * phi->GradBF2D(j, 1, ttx, tty);
		  double temprod = pw[l][1] * kk->Eval(bcoord) * (drx*nn[0] + dry*nn[1]) * 0.5 * len;

		  if(k==0)
		    btem[k][j] += phi->BF2D(0, ttx, tty) * temprod;
		  else if(k==1)
		    btem[k][j] += phi->BF2D(3, ttx, tty) * temprod;
		  else if(k==2)
		    btem[k][j] += phi->BF2D(1, ttx, tty) * temprod;
		  else if(k==3)
		    btem[k][j] += phi->BF2D(1, ttx, tty) * temprod;
		  else if(k==4)
		    btem[k][j] += phi->BF2D(4, ttx, tty) * temprod;
		  else if(k==5)
		    btem[k][j] += phi->BF2D(2, ttx, tty) * temprod;
		  else if(k==6)
		    btem[k][j] += phi->BF2D(2, ttx, tty) * temprod;
		  else if(k==7)
		    btem[k][j] += phi->BF2D(5, ttx, tty) * temprod;
		  else
		    btem[k][j] += phi->BF2D(0, ttx, tty) * temprod;
		}
	      //tflux[i][k] += pw[l][1] * kk->Eval(bcoord) * ( uux(tx, ty)*nn[0] + uuy(tx, ty)*nn[1] ) * 0.5 * len;
	    }
	  for(int l=0; l<6; l++) tflux[i][k] += btem[k][l] * locsol(l);	  
	}

      for (int k=0; k<jit.Size(); k++) { delete []jit[k]; }	    
      for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
      for (int k=0; k<btem.Size(); k++) { delete []btem[k]; }
    }

  Vector cverr(dualmesh->GetDualMeshNumDOF());
  cverr = 0.0;
  
  Vector cverruh(dualmesh->GetDualMeshNumDOF());
  cverruh = 0.0;
  
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      fem->getElementDOF(i, ind);
      Array<int> edges;
      mesh->GetElementEdges(i, edges);
      
      Array<double *> coord(mesh->GetDim());
      for (int k=0; k<coord.Size(); k++) coord[k] = new double[3];
      
      mesh->GetElementVerticesCoord(i, coord);

      double t1 = coord[0][1] - coord[0][0];
      double t2 = coord[0][2] - coord[0][0];
      double t3 = coord[1][1] - coord[1][0];
      double t4 = coord[1][2] - coord[1][0];
      double detJ = t1 * t4 - t2 * t3;
      
      // Jacobian Inverse Transpose * detJ
      Array<double *> jit(2);
      jit[0] = new double[2];
      jit[1] = new double[2];
      jit[0][0] = coord[1][2] - coord[1][0];
      jit[0][1] = coord[1][0] - coord[1][1];
      jit[1][0] = coord[0][0] - coord[0][2];
      jit[1][1] = coord[0][1] - coord[0][0];

      for(int k=0; k<2; k++)
	for(int j=0; j<2; j++)
	  jit[k][j] /= detJ;
      
      dualmesh->GetElemDualInfo(i, cod, nmls, lens, areas);
      
      double xa, xb, ya, yb, tx, ty, ttx, tty;

      BasisFunctions *phi = new BasisFunctions(2);
      Function *kk = data->GetEllipticFunction();
      Function *ff = data->GetForceFunction();

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
	      xa = cod[k][0];
	      xb = cod[3][0];
	      ya = cod[k][1];
	      yb = cod[3][1];
	    }
	  else if(k<6)
	    {
	      rm = (k+1) % 3;
	      xa = cod[k+1][0];
	      xb = cod[7][0];
	      ya = cod[k+1][1];
	      yb = cod[7][1];
	    }
	  else if(k<9)
	    {
	      rm = (k-1) % 3;
	      xa = cod[k+2][0];
	      xb = cod[11][0];
	      ya = cod[k+2][1];
	      yb = cod[11][1];
	    }
	  else
	    {
	      rm = k % 3;
	      if(k==9)
		{
		  xa = cod[9][0];		  
		  ya = cod[9][1];
		}
	      else if(k==10)
		{
		  xa = cod[1][0];		  
		  ya = cod[1][1];
		}
	      else
		{
		  xa = cod[5][0];		  
		  ya = cod[5][1];
		}
	      xb = cod[12][0];
	      yb = cod[12][1];
	    }
	  
	  for(int l=0; l<n; l++)
	    {
	      tx = 0.5 * (xb - xa) * pw[l][0] + 0.5 * (xa + xb);
	      ty = 0.5 * (yb - ya) * pw[l][0] + 0.5 * (ya + yb);
	      ttx = jit[0][0] * (tx - coord[0][0]) + jit[1][0] * (ty - coord[1][0]);
	      tty = jit[0][1] * (tx - coord[0][0]) + jit[1][1] * (ty - coord[1][0]);

	      double bcoord[2];
	      bcoord[0] = tx;
	      bcoord[1] = ty;
	      for(int j=0; j<6; j++)
		{
		  double drx = jit[0][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[0][1] * phi->GradBF2D(j, 1, ttx, tty);
		  double dry = jit[1][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[1][1] * phi->GradBF2D(j, 1, ttx, tty);
		  qtem[k][j] += pw[l][1] * kk->Eval(bcoord) * (drx*nmls[rm][0] + dry*nmls[rm][1]) * 0.5 * lens[rm];
		}
	    }
	} 

      SparseMatrix *locmat = new SparseMatrix(6,6);

      for(int j=0; j<6; j++)
	locmat->Elem(0,j) = qtem[2][j] - qtem[0][j];
      for(int j=0; j<6; j++)
	locmat->Elem(1,j) = qtem[5][j] - qtem[3][j];
      for(int j=0; j<6; j++)
	locmat->Elem(2,j) = qtem[8][j] - qtem[6][j];

      for(int j=0; j<6; j++)
	locmat->Elem(3,j) = qtem[0][j] - qtem[1][j] + qtem[4][j] - qtem[5][j] + qtem[11][j] - qtem[10][j];
      for(int j=0; j<6; j++)
	locmat->Elem(4,j) = qtem[3][j] - qtem[4][j] + qtem[7][j] - qtem[8][j] + qtem[9][j] - qtem[11][j];
      for(int j=0; j<6; j++)
	locmat->Elem(5,j) = qtem[6][j] - qtem[7][j] + qtem[1][j] - qtem[2][j] + qtem[10][j] - qtem[9][j];
            
      locmat->Finalize();
      //locmat->Print();

      /*
      for(int j=0; j<6; j++)
	{
	  for(int k=0; k<6; k++)
	    cout<<setprecision(16)<<locmat->Elem(j, k)<<"\t";
	  cout<<endl;
	}
      cout<<endl<<endl;
      */
	
      Vector Q(6);
      Q = 0.0;
      Array<double *> femat(6);
      Vector locsol(6);
      for(int k=0; k<6; k++)
	{
	  femat[k] = new double[6];
	  locsol(k) = sol(ind[k]);
	}

      fem->GetEllipticLocalSystem(i, ind, femat);
 
      for(int k=0; k<6; k++)
	for(int j=0; j<6; j++)
	  Q(k) += femat[k][j] * locsol(j);
	  
      Array<double> F(6);
      fem->GetForceLocalSystem(i, ind, F);

      Vector f(6);
      f = 0.0;

      Array<double> x(3);
      Array<double> y(3);
      double gpx, gpy, xx, yy;

      for(int l=0; l<6; l++)
	{
	  if(l==0)
	    {
	      x[0] = coord[0][0];
	      y[0] = coord[1][0];
	      x[1] = cod[0][0];
	      y[1] = cod[0][1];
	      x[2] = cod[2][0];
	      y[2] = cod[2][1];
	    }
	  else if(l==1)
	    {
	      x[0] = cod[3][0];
	      y[0] = cod[3][1];
	      x[1] = cod[2][0];
	      y[1] = cod[2][1];
	      x[2] = cod[0][0];
	      y[2] = cod[0][1];
	    }
	  else if(l==2)
	    {
	      x[0] = coord[0][1];
	      y[0] = coord[1][1];
	      x[1] = cod[4][0];
	      y[1] = cod[4][1];
	      x[2] = cod[6][0];
	      y[2] = cod[6][1];
	    }
	  else if(l==3)
	    {
	      x[0] = cod[7][0];
	      y[0] = cod[7][1];
	      x[1] = cod[6][0];
	      y[1] = cod[6][1];
	      x[2] = cod[4][0];
	      y[2] = cod[4][1];
	    }
	  else if(l==4)
	    {
	      x[0] = coord[0][2];
	      y[0] = coord[1][2];
	      x[1] = cod[8][0];
	      y[1] = cod[8][1];
	      x[2] = cod[10][0];
	      y[2] = cod[10][1];
	    }
	  else
	    {
	      x[0] = cod[11][0];
	      y[0] = cod[11][1];
	      x[1] = cod[10][0];
	      y[1] = cod[10][1];
	      x[2] = cod[8][0];
	      y[2] = cod[8][1];
	    }
	  
	  t1 = x[1] - x[0];
	  t2 = x[2] - x[0];
	  t3 = y[1] - y[0];
	  t4 = y[2] - y[0];
	  detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<np; k++)
	    {
	      gpx = tpw[k][0];
	      gpy = tpw[k][1];
	      xx = x[0] + t1*gpx + t2*gpy;
	      yy = y[0] + t3*gpx + t4*gpy;
	      
	      double bcoord[2];
	      bcoord[0] = xx;
	      bcoord[1] = yy;
	      f(l/2) += tpw[k][2] * ff->Eval(bcoord) * detJ;
	    } 
	}

      for(int l=0; l<18; l++)
	{
	  if(l==0)
	    {
	      x[0] = 0.5 * (coord[0][0] + coord[0][1]);
	      y[0] = 0.5 * (coord[1][0] + coord[1][1]);
	      x[1] = cod[1][0];
	      y[1] = cod[1][1];
	      x[2] = cod[0][0];
	      y[2] = cod[0][1];
	    }
	  else if(l==1)
	    {
	      x[0] = cod[3][0];
	      y[0] = cod[3][1];
	      x[1] = cod[0][0];
	      y[1] = cod[0][1];
	      x[2] = cod[1][0];
	      y[2] = cod[1][1];
	    }
	  else if(l==2)
	    {
	      x[0] = 0.5 * (coord[0][0] + coord[0][1]);
	      y[0] = 0.5 * (coord[1][0] + coord[1][1]);
	      x[1] = cod[6][0];
	      y[1] = cod[6][1];
	      x[2] = cod[5][0];
	      y[2] = cod[5][1];
	    }
	  else if(l==3)
	    {
	      x[0] = cod[7][0];
	      y[0] = cod[7][1];
	      x[1] = cod[5][0];
	      y[1] = cod[5][1];
	      x[2] = cod[6][0];
	      y[2] = cod[6][1];
	    }
	  else if(l==4)
	    {
	      x[0] = 0.5 * (coord[0][0] + coord[0][1]);
	      y[0] = 0.5 * (coord[1][0] + coord[1][1]);
	      x[1] = cod[5][0];
	      y[1] = cod[5][1];
	      x[2] = cod[1][0];
	      y[2] = cod[1][1];
	    }
	  else if(l==5)
	    {
	      x[0] = cod[12][0];
	      y[0] = cod[12][1];
	      x[1] = cod[1][0];
	      y[1] = cod[1][1];
	      x[2] = cod[5][0];
	      y[2] = cod[5][1];
	    }
	  else if(l==6)
	    {
	      x[0] = 0.5 * (coord[0][2] + coord[0][1]);
	      y[0] = 0.5 * (coord[1][2] + coord[1][1]);
	      x[1] = cod[5][0];
	      y[1] = cod[5][1];
	      x[2] = cod[4][0];
	      y[2] = cod[4][1];
	    }
	  else if(l==7)
	    {
	      x[0] = cod[7][0];
	      y[0] = cod[7][1];
	      x[1] = cod[4][0];
	      y[1] = cod[4][1];
	      x[2] = cod[5][0];
	      y[2] = cod[5][1];
	    }
	  else if(l==8)
	    {
	      x[0] = 0.5 * (coord[0][2] + coord[0][1]);
	      y[0] = 0.5 * (coord[1][2] + coord[1][1]);
	      x[1] = cod[10][0];
	      y[1] = cod[10][1];
	      x[2] = cod[9][0];
	      y[2] = cod[9][1];
	    }
	  else if(l==9)
	    {
	      x[0] = cod[11][0];
	      y[0] = cod[11][1];
	      x[1] = cod[9][0];
	      y[1] = cod[9][1];
	      x[2] = cod[10][0];
	      y[2] = cod[10][1];
	    }
	  else if(l==10)
	    {
	      x[0] = 0.5 * (coord[0][2] + coord[0][1]);
	      y[0] = 0.5 * (coord[1][2] + coord[1][1]);
	      x[1] = cod[9][0];
	      y[1] = cod[9][1];
	      x[2] = cod[5][0];
	      y[2] = cod[5][1];
	    }
	  else if(l==11)
	    {
	      x[0] = cod[12][0];
	      y[0] = cod[12][1];
	      x[1] = cod[5][0];
	      y[1] = cod[5][1];
	      x[2] = cod[9][0];
	      y[2] = cod[9][1];
	    }
	  else if(l==12)
	    {
	      x[0] = 0.5 * (coord[0][2] + coord[0][0]);
	      y[0] = 0.5 * (coord[1][2] + coord[1][0]);
	      x[1] = cod[2][0];
	      y[1] = cod[2][1];
	      x[2] = cod[1][0];
	      y[2] = cod[1][1];
	    }
	  else if(l==13)
	    {
	      x[0] = cod[3][0];
	      y[0] = cod[3][1];
	      x[1] = cod[1][0];
	      y[1] = cod[1][1];
	      x[2] = cod[2][0];
	      y[2] = cod[2][1];

	    }
	  else if(l==14)
	    {
	      x[0] = 0.5 * (coord[0][2] + coord[0][0]);
	      y[0] = 0.5 * (coord[1][2] + coord[1][0]);
	      x[1] = cod[9][0];
	      y[1] = cod[9][1];
	      x[2] = cod[8][0];
	      y[2] = cod[8][1];
	    }
	  else if(l==15)
	    {
	      x[0] = cod[11][0];
	      y[0] = cod[11][1];
	      x[1] = cod[8][0];
	      y[1] = cod[8][1];
	      x[2] = cod[9][0];
	      y[2] = cod[9][1];
	    }
	  else if(l==16)
	    {
	      x[0] = 0.5 * (coord[0][2] + coord[0][0]);
	      y[0] = 0.5 * (coord[1][2] + coord[1][0]);
	      x[1] = cod[1][0];
	      y[1] = cod[1][1];
	      x[2] = cod[9][0];
	      y[2] = cod[9][1];
	    }
	  else
	    {
	      x[0] = cod[12][0];
	      y[0] = cod[12][1];
	      x[1] = cod[9][0];
	      y[1] = cod[9][1];
	      x[2] = cod[1][0];
	      y[2] = cod[1][1];
	    }

	  t1 = x[1] - x[0];
	  t2 = x[2] - x[0];
	  t3 = y[1] - y[0];
	  t4 = y[2] - y[0];
	  detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<np; k++)
	    {
	      gpx = tpw[k][0];
	      gpy = tpw[k][1];
	      xx = x[0] + t1*gpx + t2*gpy;
	      yy = y[0] + t3*gpx + t4*gpy;
	      
	      double bcoord[2];
	      bcoord[0] = xx;
	      bcoord[1] = yy;
	      f(3+l/6) += tpw[k][2] * ff->Eval(bcoord) * detJ;
	    } 
	}

      Vector RHS(6);
      for(int k=0; k<6; k++) RHS(k) = Q(k) - F[k] + f(k);

      Array<int> eind;
      if(i%2==0)
	{
	  mesh->GetEdgeElements(edges[0], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[1]][5] - tflux[i][0] + tflux[eind[1]][5] );
		  RHS(3) += 0.5 * ( sflux[i][1] - sflux[eind[1]][4] - tflux[i][1] + tflux[eind[1]][4] ); 
		  RHS(1) += 0.5 * ( sflux[i][2] - sflux[eind[1]][3] - tflux[i][2] + tflux[eind[1]][3] ); 
		}
	      else
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[0]][5] - tflux[i][0] + tflux[eind[0]][5] );
		  RHS(3) += 0.5 * ( sflux[i][1] - sflux[eind[0]][4] - tflux[i][1] + tflux[eind[0]][4] ); 
		  RHS(1) += 0.5 * ( sflux[i][2] - sflux[eind[0]][3] - tflux[i][2] + tflux[eind[0]][3] ); 
		}
	    }
	  else
	    {
	      RHS(0) += -tflux[i][0] + sflux[i][0];
	      RHS(3) += -tflux[i][1] + sflux[i][1];
	      RHS(1) += -tflux[i][2] + sflux[i][2];
	    }

	  mesh->GetEdgeElements(edges[1], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(1) += 0.5 * ( sflux[i][3] - sflux[eind[1]][8] - tflux[i][3] + tflux[eind[1]][8] );
		  RHS(4) += 0.5 * ( sflux[i][4] - sflux[eind[1]][7] - tflux[i][4] + tflux[eind[1]][7] ); 
		  RHS(2) += 0.5 * ( sflux[i][5] - sflux[eind[1]][6] - tflux[i][5] + tflux[eind[1]][6] ); 
		}
	      else
		{
		  RHS(1) += 0.5 * ( sflux[i][3] - sflux[eind[0]][8] - tflux[i][3] + tflux[eind[0]][8] );
		  RHS(4) += 0.5 * ( sflux[i][4] - sflux[eind[0]][7] - tflux[i][4] + tflux[eind[0]][7] ); 
		  RHS(2) += 0.5 * ( sflux[i][5] - sflux[eind[0]][6] - tflux[i][5] + tflux[eind[0]][6] ); 
		}
	    }
	  else
	    {
	      RHS(1) += -tflux[i][3] + sflux[i][3];
	      RHS(4) += -tflux[i][4] + sflux[i][4];
	      RHS(2) += -tflux[i][4] + sflux[i][5];
	    }

	  mesh->GetEdgeElements(edges[2], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(2) += 0.5 * ( sflux[i][6] - sflux[eind[1]][2] - tflux[i][6] + tflux[eind[1]][2] );
		  RHS(5) += 0.5 * ( sflux[i][7] - sflux[eind[1]][1] - tflux[i][7] + tflux[eind[1]][1] ); 
		  RHS(0) += 0.5 * ( sflux[i][8] - sflux[eind[1]][0] - tflux[i][8] + tflux[eind[1]][0] ); 
		}
	      else
		{
		  RHS(2) += 0.5 * ( sflux[i][6] - sflux[eind[0]][2] - tflux[i][6] + tflux[eind[0]][2] );
		  RHS(5) += 0.5 * ( sflux[i][7] - sflux[eind[0]][1] - tflux[i][7] + tflux[eind[0]][1] ); 
		  RHS(0) += 0.5 * ( sflux[i][8] - sflux[eind[0]][0] - tflux[i][8] + tflux[eind[0]][0] ); 
		}
	    }
	  else
	    {
	      RHS(2) += -tflux[i][6] + sflux[i][6];
	      RHS(5) += -tflux[i][7] + sflux[i][7];
	      RHS(0) += -tflux[i][8] + sflux[i][8];
	    }
	}
      else
	{
	  mesh->GetEdgeElements(edges[0], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[1]][8] - tflux[i][0] + tflux[eind[1]][8] );
		  RHS(3) += 0.5 * ( sflux[i][1] - sflux[eind[1]][7] - tflux[i][1] + tflux[eind[1]][7] ); 
		  RHS(1) += 0.5 * ( sflux[i][2] - sflux[eind[1]][6] - tflux[i][2] + tflux[eind[1]][6] ); 
		}
	      else
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[0]][8] - tflux[i][0] + tflux[eind[0]][8] );
		  RHS(3) += 0.5 * ( sflux[i][1] - sflux[eind[0]][7] - tflux[i][1] + tflux[eind[0]][7] ); 
		  RHS(1) += 0.5 * ( sflux[i][2] - sflux[eind[0]][6] - tflux[i][2] + tflux[eind[0]][6] ); 
		}
	    }
	  else
	    {
	      RHS(0) += -tflux[i][0] + sflux[i][0];
	      RHS(3) += -tflux[i][1] + sflux[i][1];
	      RHS(1) += -tflux[i][2] + sflux[i][2];
	    }

	  mesh->GetEdgeElements(edges[1], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(1) += 0.5 * ( sflux[i][3] - sflux[eind[1]][2] - tflux[i][3] + tflux[eind[1]][2] );
		  RHS(4) += 0.5 * ( sflux[i][4] - sflux[eind[1]][1] - tflux[i][4] + tflux[eind[1]][1] ); 
		  RHS(2) += 0.5 * ( sflux[i][5] - sflux[eind[1]][0] - tflux[i][5] + tflux[eind[1]][0] ); 
		}
	      else
		{
		  RHS(1) += 0.5 * ( sflux[i][3] - sflux[eind[0]][2] - tflux[i][3] + tflux[eind[0]][2] );
		  RHS(4) += 0.5 * ( sflux[i][4] - sflux[eind[0]][1] - tflux[i][4] + tflux[eind[0]][1] ); 
		  RHS(2) += 0.5 * ( sflux[i][5] - sflux[eind[0]][0] - tflux[i][5] + tflux[eind[0]][0] ); 
		}
	    }
	  else
	    {
	      RHS(1) += -tflux[i][3] + sflux[i][3];
	      RHS(4) += -tflux[i][4] + sflux[i][4];
	      RHS(2) += -tflux[i][5] + sflux[i][5];
	    }

	  mesh->GetEdgeElements(edges[2], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(2) += 0.5 * ( sflux[i][6] - sflux[eind[1]][5] - tflux[i][6] + tflux[eind[1]][5] );
		  RHS(5) += 0.5 * ( sflux[i][7] - sflux[eind[1]][4] - tflux[i][7] + tflux[eind[1]][4] ); 
		  RHS(0) += 0.5 * ( sflux[i][8] - sflux[eind[1]][3] - tflux[i][8] + tflux[eind[1]][3] ); 
		}
	      else
		{
		  RHS(2) += 0.5 * ( sflux[i][6] - sflux[eind[0]][5] - tflux[i][6] + tflux[eind[0]][5] );
		  RHS(5) += 0.5 * ( sflux[i][7] - sflux[eind[0]][4] - tflux[i][7] + tflux[eind[0]][4] ); 
		  RHS(0) += 0.5 * ( sflux[i][8] - sflux[eind[0]][3] - tflux[i][8] + tflux[eind[0]][3] ); 
		}
	    }
	  else
	    {
	      RHS(2) += -tflux[i][6] + sflux[i][6];
	      RHS(5) += -tflux[i][7] + sflux[i][7];
	      RHS(0) += -tflux[i][8] + sflux[i][8];
	    }
	}

      Vector s(6);
      s = 0.0;
      Vector rhs(6); rhs = RHS;
      RectangularMatrix localconsmat(6);      
      for(int j=0; j<6; j++)
	{
	  for(int k=0; k<6; k++)
	    localconsmat(j,k) = locmat->Elem(j,k);
	}
      
      for(int j=0; j<6; j++) localconsmat(3, j) = 0.0;
      localconsmat(3,3) = 1.0;
      rhs(3) = locsol(3);

      DenseMatrixInverse invmat(localconsmat);
      invmat.Mult(rhs, s);
      //s.Print();
      
      Vector sss(6);
      sss = 0.0;
      locmat->Mult(s, sss);
      for(int k=0; k<6; k++)
	cverr(ind[k]) += sss(k) - f(k);
      locmat->Mult(locsol, sss);
      for(int k=0; k<6; k++)
	cverruh(ind[k]) += sss(k) - f(k);
      //cverr(ind[k]) += Q(k) - F[k];
	
      for(int k=0; k<6; k++) { ppsol[k][i] = s(k); }
      for(int k=0; k<12; k++)
	{
	  for(int j=0; j<6; j++)
	    {
	      ppflux[k][i] += -qtem[k][j] * s(j);
	      flux[k][i] += -qtem[k][j] * locsol(j);
	    }
	}
      
      for(int k=0; k<femat.Size(); k++) delete []femat[k];
      for(int k=0; k<jit.Size(); k++) delete []jit[k];
      for(int k=0; k<coord.Size(); k++) { delete []coord[k]; }
      for(int k=0; k<qtem.Size(); k++) { delete []qtem[k]; }
      delete locmat;
    }      

  for(int i=0; i<mesh->GetNBE(); i++)
    {
      Array<int> bdrind(2);
      mesh->GetBdrElementVertices(i, bdrind);
      bdrind.Append(mesh->GetNV() +  mesh->GetBdrElementEdgeIndex(i));

      if( fabs(cverr(bdrind[0])) > 1.0e-16 )
	{
	  ppbdrflux[bdrind[0]] = -cverr(bdrind[0]);
	  cverr(bdrind[0]) = 0.0;
	  cverruh(bdrind[0]) = 0.0;
	}
      if( fabs(cverr(bdrind[1])) > 1.0e-16 )
	{
	  ppbdrflux[bdrind[1]] = -cverr(bdrind[1]);
	  cverr(bdrind[1]) = 0.0;
	  cverruh(bdrind[1]) = 0.0;
	}
      
      ppbdrflux[bdrind[2]] = -cverr(bdrind[2]);
      cverr(bdrind[2]) = 0.0;
      cverruh(bdrind[2]) = 0.0;
    }

  ofstream fileout("qlce.out");
  for(int i=0; i<dualmesh->GetDualMeshNumDOF(); i++)
    {
      fileout<<i<<"\t"<<cverruh(i)<< "\t" << cverr(i) << endl;
      cout<<setprecision(3)<<i<<"\t"<<cverruh(i)<< "\t"<<cverr(i)<<endl;
    }
  fileout.close(); 

  for(int i=0; i<tflux.Size(); i++) delete []tflux[i];
  for(int i=0; i<sflux.Size(); i++) delete []sflux[i];
  for(int i=0; i<pw.Size(); i++) delete []pw[i];
  for(int i=0; i<tpw.Size(); i++) delete []tpw[i];
  for (int k=0; k<cod.Size(); k++) { delete []cod[k]; }
  for (int k=0; k<nmls.Size(); k++) { delete []nmls[k]; }
}

//===========================================================

void Postprocessing::ComputeCubicConservativeTriangularFlux(Vector &sol, Array<double*> &flux, Array<double*> &ppsol,
							    Array<double*> &ppflux, Array<double> &ppbdrflux)
{
  mesh->GetEdgeToElementTable();
    
  int nl = NumLocalDOF;
  Array<double *> nmls(18);
  Array<double> lens(18);
  Array<double> areas(10);
  Array<double *> cod(19);
  for(int k=0; k<cod.Size(); k++) cod[k] = new double[2];
  for(int k=0; k<nmls.Size(); k++) nmls[k] = new double[2];

  int n = 5;
  Array<double *> pw(n);
  for(int i=0; i<n; i++) pw[i] = new double[2];
  GetStandLineQuadPW(n, pw);

  int np = ConvertN2PN(6);
  Array<double *> tpw(np);
  for(int i=0; i<np; i++) tpw[i] = new double[3];
  GetStandTriQuadPW(6, tpw);

  Array<double *> sflux(mesh->GetNE());
  for(int i=0; i<sflux.Size(); i++) sflux[i] = new double[12];
  for(int j=0; j<sflux.Size(); j++)
    for(int i=0; i<12; i++)
      sflux[j][i] = 0.0;

  Array<double *> tflux(mesh->GetNE());
  for(int i=0; i<tflux.Size(); i++) tflux[i] = new double[12];
  for(int j=0; j<tflux.Size(); j++)
    for(int i=0; i<12; i++)
      tflux[j][i] = 0.0;

  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      fem->getElementDOF(i, ind);
      Array<int> edges;
      mesh->GetElementEdges(i, edges);

      Array<double *> coord(mesh->GetDim());
      for (int k=0; k<coord.Size(); k++) coord[k] = new double[3];
      
      mesh->GetElementVerticesCoord(i, coord);

      double t1 = coord[0][1] - coord[0][0];
      double t2 = coord[0][2] - coord[0][0];
      double t3 = coord[1][1] - coord[1][0];
      double t4 = coord[1][2] - coord[1][0];
      double detJ = t1 * t4 - t2 * t3;
      
      // Jacobian Inverse Transpose * detJ
      Array<double *> jit(2);
      jit[0] = new double[2];
      jit[1] = new double[2];
      jit[0][0] = coord[1][2] - coord[1][0];
      jit[0][1] = coord[1][0] - coord[1][1];
      jit[1][0] = coord[0][0] - coord[0][2];
      jit[1][1] = coord[0][1] - coord[0][0];

      for(int k=0; k<2; k++)
	for(int j=0; j<2; j++)
	  jit[k][j] /= detJ;
      
      dualmesh->GetElemDualInfo(i, cod, nmls, lens, areas);

      Vector locsol(nl);
      for(int k=0; k<nl; k++) locsol(k) = sol(ind[k]);

      Array<double *> btem(12);
      for(int k=0; k<12; k++) btem[k] = new double[nl];
      for(int k=0; k<12; k++)
	for(int j=0; j<nl; j++)
	  btem[k][j] = 0.0;
      
      double xa, xb, ya, yb, tx, ty, ttx, tty;

      BasisFunctions *phi = new BasisFunctions(3); 
      Function *kk = data->GetEllipticFunction();
      for(int k=0; k<12; k++)
	{
	  double len = 0.0;
	  double nn[2];
	  if(k==0)
	    {
	      xa = coord[0][0];
	      xb = cod[0][0];
	      ya = coord[1][0];
	      yb = cod[0][1];
	    }
	  else if(k==1)
	    {
	      xa = cod[0][0];
	      xb = cod[1][0];
	      ya = cod[0][1];
	      yb = cod[1][1];
	    }
	  else if(k==2)
	    {
	      xa = cod[1][0];
	      xb = cod[2][0];
	      ya = cod[1][1];
	      yb = cod[2][1];
	    }
	  else if(k==3)
	    {
	      xa = cod[2][0];
	      xb = coord[0][1];
	      ya = cod[2][1];
	      yb = coord[1][1];
	    }
	  else if(k==4)
	    {
	      xa = coord[0][1];
	      xb = cod[3][0];
	      ya = coord[1][1];
	      yb = cod[3][1];
	    }
	  else if(k==5)
	    {
	      xa = cod[3][0];
	      xb = cod[4][0];
	      ya = cod[3][1];
	      yb = cod[4][1];
	    }
	  else if(k==6)
	    {
	      xa = cod[4][0];
	      xb = cod[5][0];
	      ya = cod[4][1];
	      yb = cod[5][1];
	    }
	  else if(k==7)
	    {
	      xa = cod[5][0];
	      xb = coord[0][2];
	      ya = cod[5][1];
	      yb = coord[1][2];
	    }
	  else if(k==8)
	    {
	      xa = coord[0][2];
	      xb = cod[6][0];
	      ya = coord[1][2];
	      yb = cod[6][1];
	    }
	  else if(k==9)
	    {
	      xa = cod[6][0];
	      xb = cod[7][0];
	      ya = cod[6][1];
	      yb = cod[7][1];
	    }
	  else if(k==10)
	    {
	      xa = cod[7][0];
	      xb = cod[8][0];
	      ya = cod[7][1];
	      yb = cod[8][1];
	    }
	  else
	    {
	      xa = cod[8][0];
	      xb = coord[0][0];
	      ya = cod[8][1];
	      yb = coord[1][0];
	    }

	  len = sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya) );
	  nn[0] = (yb - ya)/len;
	  nn[1] = (xa - xb)/len;
	  
	  for(int l=0; l<n; l++)
	    {
	      tx = 0.5 * (xb - xa) * pw[l][0] + 0.5 * (xa + xb);
	      ty = 0.5 * (yb - ya) * pw[l][0] + 0.5 * (ya + yb);
	      ttx = jit[0][0] * (tx - coord[0][0]) + jit[1][0] * (ty - coord[1][0]);
	      tty = jit[0][1] * (tx - coord[0][0]) + jit[1][1] * (ty - coord[1][0]);

	      double bcoord[2];
	      bcoord[0] = tx;
	      bcoord[1] = ty;
	      for(int j=0; j<nl; j++)
		{
		  double drx = jit[0][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[0][1] * phi->GradBF2D(j, 1, ttx, tty);
		  double dry = jit[1][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[1][1] * phi->GradBF2D(j, 1, ttx, tty);
		  btem[k][j] += pw[l][1] * kk->Eval(bcoord) * (drx*nn[0] + dry*nn[1]) * 0.5 * len;
		}
	    }
	  for(int l=0; l<nl; l++) sflux[i][k] += btem[k][l] * locsol(l);	  
	}

      for(int k=0; k<12; k++)
	for(int j=0; j<nl; j++)
	  btem[k][j] = 0.0;
      
      for(int k=0; k<12; k++)
	{
	  double len = 0.0;
	  double nn[2];
	  if(k<4)
	    {
	      xa = coord[0][0];
	      xb = coord[0][1];
	      ya = coord[1][0];
	      yb = coord[1][1];
	    }
	  else if(k<8)
	    {
	      xa = coord[0][1];
	      xb = coord[0][2];
	      ya = coord[1][1];
	      yb = coord[1][2];
	    }
	  else
	    {
	      xa = coord[0][2];
	      xb = coord[0][0];
	      ya = coord[1][2];
	      yb = coord[1][0];
	    }

	  len = sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya) );
	  nn[0] = (yb - ya)/len;
	  nn[1] = (xa - xb)/len;
	  
	  for(int l=0; l<n; l++)
	    {
	      tx = 0.5 * (xb - xa) * pw[l][0] + 0.5 * (xa + xb);
	      ty = 0.5 * (yb - ya) * pw[l][0] + 0.5 * (ya + yb);
	      ttx = jit[0][0] * (tx - coord[0][0]) + jit[1][0] * (ty - coord[1][0]);
	      tty = jit[0][1] * (tx - coord[0][0]) + jit[1][1] * (ty - coord[1][0]);

	      double bcoord[2];
	      bcoord[0] = tx;
	      bcoord[1] = ty;
	      for(int j=0; j<nl; j++)
		{
		  double drx = jit[0][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[0][1] * phi->GradBF2D(j, 1, ttx, tty);
		  double dry = jit[1][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[1][1] * phi->GradBF2D(j, 1, ttx, tty);
		  double temprod = pw[l][1] * kk->Eval(bcoord) * (drx*nn[0] + dry*nn[1]) * 0.5 * len;
		  if(k==0)
		    btem[k][j] += phi->BF2D(0, ttx, tty) * temprod;
		  else if(k==1)
		    btem[k][j] += phi->BF2D(3, ttx, tty) * temprod;
		  else if(k==2)
		    btem[k][j] += phi->BF2D(4, ttx, tty) * temprod;
		  else if(k==3)
		    btem[k][j] += phi->BF2D(1, ttx, tty) * temprod;
		  else if(k==4)
		    btem[k][j] += phi->BF2D(1, ttx, tty) * temprod;
		  else if(k==5)
		    btem[k][j] += phi->BF2D(5, ttx, tty) * temprod;
		  else if(k==6)
		    btem[k][j] += phi->BF2D(6, ttx, tty) * temprod;
		  else if(k==7)
		    btem[k][j] += phi->BF2D(2, ttx, tty) * temprod;
		  else if(k==8)
		    btem[k][j] += phi->BF2D(2, ttx, tty) * temprod;
		  else if(k==9)
		    btem[k][j] += phi->BF2D(7, ttx, tty) * temprod;
		  else if(k==10)
		    btem[k][j] += phi->BF2D(8, ttx, tty) * temprod;
		  else if(k==11)
		    btem[k][j] += phi->BF2D(0, ttx, tty) * temprod;
		}
 	    }
	  for(int l=0; l<nl; l++) tflux[i][k] += btem[k][l] * locsol(l);	  
	}

      for (int k=0; k<jit.Size(); k++) { delete []jit[k]; }
      for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
      for (int k=0; k<btem.Size(); k++) { delete []btem[k]; }
    }

  Vector cverr( dualmesh->GetDualMeshNumDOF() );
  cverr = 0.0;
  
  Vector cverruh( dualmesh->GetDualMeshNumDOF() );
  cverruh = 0.0;
  
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      fem->getElementDOF(i, ind);
      Array<int> edges;
      mesh->GetElementEdges(i, edges);

      Array<double *> coord(mesh->GetDim());
      for (int k=0; k<coord.Size(); k++) coord[k] = new double[3];
      
      mesh->GetElementVerticesCoord(i, coord);

      double t1 = coord[0][1] - coord[0][0];
      double t2 = coord[0][2] - coord[0][0];
      double t3 = coord[1][1] - coord[1][0];
      double t4 = coord[1][2] - coord[1][0];
      double detJ = t1 * t4 - t2 * t3;
      
      // Jacobian Inverse Transpose * detJ
      Array<double *> jit(2);
      jit[0] = new double[2];
      jit[1] = new double[2];
      jit[0][0] = coord[1][2] - coord[1][0];
      jit[0][1] = coord[1][0] - coord[1][1];
      jit[1][0] = coord[0][0] - coord[0][2];
      jit[1][1] = coord[0][1] - coord[0][0];

      for(int k=0; k<2; k++)
	for(int j=0; j<2; j++)
	  jit[k][j] /= detJ;
      
      dualmesh->GetElemDualInfo(i, cod, nmls, lens, areas);

      double xa, xb, ya, yb, tx, ty, ttx, tty;
      
      BasisFunctions *phi = new BasisFunctions(3);
      Function *kk = data->GetEllipticFunction();
      Function *ff = data->GetForceFunction();

      Array<double *> ctem(18);
      for(int k=0; k<18; k++) ctem[k] = new double[nl];
      for(int k=0; k<18; k++)
	for(int j=0; j<nl; j++)
	  ctem[k][j] = 0.0;

      for(int k=0; k<18; k++)
	{
	  if(k==0)
	    {
	      xa = cod[9][0];
	      xb = cod[0][0];
	      ya = cod[9][1];
	      yb = cod[0][1];
	    }
	  else if(k==1)
	    {
	      xa = cod[9][0];
	      xb = cod[8][0];
	      ya = cod[9][1];
	      yb = cod[8][1];
	    }
	  else if(k==9)
	    {
	      xa = cod[9][0];
	      xb = cod[12][0];
	      ya = cod[9][1];
	      yb = cod[12][1];
	    }
	  else if(k==2)
	    {
	      xa = cod[10][0];
	      xb = cod[3][0];
	      ya = cod[10][1];
	      yb = cod[3][1];
	    }
	  else if(k==3)
	    {
	      xa = cod[10][0];
	      xb = cod[2][0];
	      ya = cod[10][1];
	      yb = cod[2][1];
	    }
	  else if(k==10)
	    {
	      xa = cod[10][0];
	      xb = cod[14][0];
	      ya = cod[10][1];
	      yb = cod[14][1];
	    }
	  else if(k==4)
	    {
	      xa = cod[11][0];
	      xb = cod[6][0];
	      ya = cod[11][1];
	      yb = cod[6][1];
	    }
	  else if(k==5)
	    {
	      xa = cod[11][0];
	      xb = cod[5][0];
	      ya = cod[11][1];
	      yb = cod[5][1];
	    }
	  else if(k==11)
	    {
	      xa = cod[11][0];
	      xb = cod[16][0];
	      ya = cod[11][1];
	      yb = cod[16][1];
	    }
	  else if(k==6)
	    {
	      xa = cod[13][0];
	      xb = cod[1][0];
	      ya = cod[13][1];
	      yb = cod[1][1];
	    }
	  else if(k==12)
	    {
	      xa = cod[13][0];
	      xb = cod[12][0];
	      ya = cod[13][1];
	      yb = cod[12][1];
	    }
	  else if(k==13)
	    {
	      xa = cod[13][0];
	      xb = cod[14][0];
	      ya = cod[13][1];
	      yb = cod[14][1];
	    }
	  else if(k==7)
	    {
	      xa = cod[15][0];
	      xb = cod[4][0];
	      ya = cod[15][1];
	      yb = cod[4][1];
	    }
	  else if(k==14)
	    {
	      xa = cod[15][0];
	      xb = cod[14][0];
	      ya = cod[15][1];
	      yb = cod[14][1];
	    }
	  else if(k==15)
	    {
	      xa = cod[15][0];
	      xb = cod[16][0];
	      ya = cod[15][1];
	      yb = cod[16][1];
	    }
	  else if(k==8)
	    {
	      xa = cod[17][0];
	      xb = cod[7][0];
	      ya = cod[17][1];
	      yb = cod[7][1];
	    }
	  else if(k==16)
	    {
	      xa = cod[17][0];
	      xb = cod[16][0];
	      ya = cod[17][1];
	      yb = cod[16][1];
	    }
	  else
	    {
	      xa = cod[17][0];
	      xb = cod[12][0];
	      ya = cod[17][1];
	      yb = cod[12][1];
	    }
	  
	  for(int l=0; l<n; l++)
	    {
	      tx = 0.5 * (xb - xa) * pw[l][0] + 0.5 * (xa + xb);
	      ty = 0.5 * (yb - ya) * pw[l][0] + 0.5 * (ya + yb);
	      ttx = jit[0][0] * (tx - coord[0][0]) + jit[1][0] * (ty - coord[1][0]);
	      tty = jit[0][1] * (tx - coord[0][0]) + jit[1][1] * (ty - coord[1][0]);
	      
	      double bcoord[2];
	      bcoord[0] = tx;
	      bcoord[1] = ty;
	      for(int j=0; j<nl; j++)
		{
		  double drx = jit[0][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[0][1] * phi->GradBF2D(j, 1, ttx, tty);
		  double dry = jit[1][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[1][1] * phi->GradBF2D(j, 1, ttx, tty);
		  ctem[k][j] += pw[l][1] * kk->Eval(bcoord) * (drx*nmls[k][0] + dry*nmls[k][1]) * 0.5 * lens[k];
		}
	    }
	}
      
      SparseMatrix *locmat = new SparseMatrix(nl,nl);

      for(int j=0; j<nl; j++)
	locmat->Elem(0,j) = ctem[1][j] - ctem[0][j];
      for(int j=0; j<nl; j++)
	locmat->Elem(1,j) = ctem[3][j] - ctem[2][j];
      for(int j=0; j<nl; j++)
	locmat->Elem(2,j) = ctem[5][j] - ctem[4][j];
      
      for(int j=0; j<nl; j++)
	locmat->Elem(3,j) = ctem[0][j] - ctem[9][j] + ctem[12][j] - ctem[6][j];
      for(int j=0; j<nl; j++)
	locmat->Elem(4,j) = ctem[6][j] - ctem[13][j] + ctem[10][j] - ctem[3][j];

      for(int j=0; j<nl; j++)
	locmat->Elem(5,j) = ctem[2][j] - ctem[10][j] + ctem[14][j] - ctem[7][j];
      for(int j=0; j<nl; j++)
	locmat->Elem(6,j) = ctem[7][j] - ctem[15][j] + ctem[11][j] - ctem[5][j];

      for(int j=0; j<nl; j++)
	locmat->Elem(7,j) = ctem[4][j] - ctem[11][j] + ctem[16][j] - ctem[8][j];
      for(int j=0; j<nl; j++)
	locmat->Elem(8,j) = ctem[8][j] - ctem[17][j] + ctem[9][j] - ctem[1][j];

      for(int j=0; j<nl; j++)
	locmat->Elem(9,j) = ctem[17][j] - ctem[12][j] + ctem[13][j] - ctem[14][j] + ctem[15][j] - ctem[16][j];

      locmat->Finalize();
      //locmat->Print();

      Vector Q(nl);
      Q = 0.0;
      Array<double *> femat(nl);
      Vector locsol(nl);
      for(int k=0; k<nl; k++)
	{
	  femat[k] = new double[nl];
	  locsol(k) = sol(ind[k]);
	}

      fem->GetEllipticLocalSystem(i, ind, femat);
 
      for(int k=0; k<nl; k++)
	for(int j=0; j<nl; j++)
	  Q(k) += femat[k][j] * locsol(j);
	  
      Array<double> F(nl);
      fem->GetForceLocalSystem(i, ind, F);

      Vector f(nl);
      f = 0.0;

      Array<double> x(3);
      Array<double> y(3);
      double gpx, gpy, xx, yy;

      for(int l=0; l<6; l++)
	{
	  if(l==0)
	    {
	      x[0] = coord[0][0];
	      y[0] = coord[1][0];
	      x[1] = cod[0][0];
	      y[1] = cod[0][1];
	      x[2] = cod[8][0];
	      y[2] = cod[8][1];
	    }
	  else if(l==1)
	    {
	      x[0] = cod[0][0];
	      y[0] = cod[0][1];
	      x[2] = cod[8][0];
	      y[2] = cod[8][1];
	      x[1] = cod[9][0];
	      y[1] = cod[9][1];
	    }
	  else if(l==2)
	    {
	      x[1] = coord[0][1];
	      y[1] = coord[1][1];
	      x[0] = cod[2][0];
	      y[0] = cod[2][1];
	      x[2] = cod[3][0];
	      y[2] = cod[3][1];
	    }
	  else if(l==3)
	    {
	      x[1] = cod[2][0];
	      y[1] = cod[2][1];
	      x[2] = cod[3][0];
	      y[2] = cod[3][1];
	      x[0] = cod[10][0];
	      y[0] = cod[10][1];
	    }
	  else if(l==4)
	    {
	      x[0] = coord[0][2];
	      y[0] = coord[1][2];
	      x[1] = cod[6][0];
	      y[1] = cod[6][1];
	      x[2] = cod[5][0];
	      y[2] = cod[5][1];
	    }
	  else
	    {
	      x[0] = cod[6][0];
	      y[0] = cod[6][1];
	      x[2] = cod[5][0];
	      y[2] = cod[5][1];
	      x[1] = cod[11][0];
	      y[1] = cod[11][1];
	    }
	  
	  t1 = x[1] - x[0];
	  t2 = x[2] - x[0];
	  t3 = y[1] - y[0];
	  t4 = y[2] - y[0];
	  detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<np; k++)
	    {
	      gpx = tpw[k][0];
	      gpy = tpw[k][1];
	      xx = x[0] + t1*gpx + t2*gpy;
	      yy = y[0] + t3*gpx + t4*gpy;
	      
	      double bcoord[2];
	      bcoord[0] = xx;
	      bcoord[1] = yy;
	      f(l/2) += tpw[k][2] * ff->Eval(bcoord) * detJ;
	    } 
	}

      for(int l=0; l<24; l++)
	{
	  if(l==0)
	    {
	      x[1] = coord[0][0]*2.0/3.0 + coord[0][1]/3.0;
	      y[1] = coord[1][0]*2.0/3.0 + coord[1][1]/3.0;
	      x[0] = cod[0][0];
	      y[0] = cod[0][1];
	      x[2] = cod[12][0];
	      y[2] = cod[12][1];
	    }
	  else if(l==1)
	    {
	      x[0] = cod[0][0];
	      y[0] = cod[0][1];
	      x[1] = cod[12][0];
	      y[1] = cod[12][1];
	      x[2] = cod[9][0];
	      y[2] = cod[9][1];
	    }
	  else if(l==2)
	    {
	      x[1] = cod[1][0];
	      y[1] = cod[1][1];
	      x[2] = cod[12][0];
	      y[2] = cod[12][1];
	      x[0] = coord[0][0]*2.0/3.0 + coord[0][1]/3.0;
	      y[0] = coord[1][0]*2.0/3.0 + coord[1][1]/3.0;
	    }
	  else if(l==3)
	    {
	      x[0] = cod[1][0];
	      y[0] = cod[1][1];
	      x[2] = cod[12][0];
	      y[2] = cod[12][1];
	      x[1] = cod[13][0];
	      y[1] = cod[13][1];
	    }
	  else if(l==4)
	    {
	      x[1] = cod[2][0];
	      y[1] = cod[2][1];
	      x[2] = cod[14][0];
	      y[2] = cod[14][1];
	      x[0] = coord[0][1]*2.0/3.0 + coord[0][0]/3.0;
	      y[0] = coord[1][1]*2.0/3.0 + coord[1][0]/3.0;
	    }
	  else if(l==5)
	    {
	      x[0] = cod[2][0];
	      y[0] = cod[2][1];
	      x[2] = cod[14][0];
	      y[2] = cod[14][1];
	      x[1] = cod[10][0];
	      y[1] = cod[10][1];
	    }
	  else if(l==6)
	    {
	      x[0] = cod[1][0];
	      y[0] = cod[1][1];
	      x[2] = cod[14][0];
	      y[2] = cod[14][1];
	      x[1] = coord[0][1]*2.0/3.0 + coord[0][0]/3.0;
	      y[1] = coord[1][1]*2.0/3.0 + coord[1][0]/3.0;
	    }
	  else if(l==7)
	    {
	      x[1] = cod[1][0];
	      y[1] = cod[1][1];
	      x[2] = cod[14][0];
	      y[2] = cod[14][1];
	      x[0] = cod[13][0];
	      y[0] = cod[13][1];
	    }
	  else if(l==8)
	    {
	      x[1] = cod[3][0];
	      y[1] = cod[3][1];
	      x[0] = cod[14][0];
	      y[0] = cod[14][1];
	      x[2] = coord[0][1]*2.0/3.0 + coord[0][2]/3.0;
	      y[2] = coord[1][1]*2.0/3.0 + coord[1][2]/3.0;
	    }
	  else if(l==9)
	    {
	      x[1] = cod[3][0];
	      y[1] = cod[3][1];
	      x[2] = cod[14][0];
	      y[2] = cod[14][1];
	      x[0] = cod[10][0];
	      y[0] = cod[10][1];
	    }
	  else if(l==10)
	    {
	      x[1] = cod[4][0];
	      y[1] = cod[4][1];
	      x[2] = cod[14][0];
	      y[2] = cod[14][1];
	      x[0] = coord[0][1]*2.0/3.0 + coord[0][2]/3.0;
	      y[0] = coord[1][1]*2.0/3.0 + coord[1][2]/3.0;
	    }
	  else if(l==11)
	    {
	      x[0] = cod[4][0];
	      y[0] = cod[4][1];
	      x[2] = cod[14][0];
	      y[2] = cod[14][1];
	      x[1] = cod[15][0];
	      y[1] = cod[15][1];
	    }
	  else if(l==12)
	    {
	      x[1] = cod[5][0];
	      y[1] = cod[5][1];
	      x[2] = cod[16][0];
	      y[2] = cod[16][1];
	      x[0] = coord[0][2]*2.0/3.0 + coord[0][1]/3.0;
	      y[0] = coord[1][2]*2.0/3.0 + coord[1][1]/3.0;
	    }
	  else if(l==13)
	    {
	      x[0] = cod[5][0];
	      y[0] = cod[5][1];
	      x[2] = cod[16][0];
	      y[2] = cod[16][1];
	      x[1] = cod[11][0];
	      y[1] = cod[11][1];
	    }
	  else if(l==14)
	    {
	      x[0] = cod[4][0];
	      y[0] = cod[4][1];
	      x[2] = cod[16][0];
	      y[2] = cod[16][1];
	      x[1] = coord[0][2]*2.0/3.0 + coord[0][1]/3.0;
	      y[1] = coord[1][2]*2.0/3.0 + coord[1][1]/3.0;
	    }
	  else if(l==15)
	    {
	      x[1] = cod[4][0];
	      y[1] = cod[4][1];
	      x[2] = cod[16][0];
	      y[2] = cod[16][1];
	      x[0] = cod[15][0];
	      y[0] = cod[15][1];
	    }
	  else if(l==16)
	    {
	      x[0] = cod[6][0];
	      y[0] = cod[6][1];
	      x[2] = cod[16][0];
	      y[2] = cod[16][1];
	      x[1] = coord[0][2]*2.0/3.0 + coord[0][0]/3.0;
	      y[1] = coord[1][2]*2.0/3.0 + coord[1][0]/3.0;
	    }
	  else if(l==17)
	    {
	      x[1] = cod[6][0];
	      y[1] = cod[6][1];
	      x[2] = cod[16][0];
	      y[2] = cod[16][1];
	      x[0] = cod[11][0];
	      y[0] = cod[11][1];
	    }
	  else if(l==18)
	    {
	      x[1] = cod[7][0];
	      y[1] = cod[7][1];
	      x[2] = cod[16][0];
	      y[2] = cod[16][1];
	      x[0] = coord[0][2]*2.0/3.0 + coord[0][0]/3.0;
	      y[0] = coord[1][2]*2.0/3.0 + coord[1][0]/3.0;
	    }
	  else if(l==19)
	    {
	      x[0] = cod[7][0];
	      y[0] = cod[7][1];
	      x[2] = cod[16][0];
	      y[2] = cod[16][1];
	      x[1] = cod[17][0];
	      y[1] = cod[17][1];
	    }
	  else if(l==20)
	    {
	      x[1] = cod[8][0];
	      y[1] = cod[8][1];
	      x[2] = cod[12][0];
	      y[2] = cod[12][1];
	      x[0] = coord[0][0]*2.0/3.0 + coord[0][2]/3.0;
	      y[0] = coord[1][0]*2.0/3.0 + coord[1][2]/3.0;
	    }
	  else if(l==21)
	    {
	      x[0] = cod[8][0];
	      y[0] = cod[8][1];
	      x[2] = cod[12][0];
	      y[2] = cod[12][1];
	      x[1] = cod[9][0];
	      y[1] = cod[9][1];
	    }
	  else if(l==23)
	    {
	      x[0] = cod[7][0];
	      y[0] = cod[7][1];
	      x[2] = cod[12][0];
	      y[2] = cod[12][1];
	      x[1] = coord[0][0]*2.0/3.0 + coord[0][2]/3.0;
	      y[1] = coord[1][0]*2.0/3.0 + coord[1][2]/3.0;
	    }
	  else
	    {
	      x[1] = cod[7][0];
	      y[1] = cod[7][1];
	      x[2] = cod[12][0];
	      y[2] = cod[12][1];
	      x[0] = cod[17][0];
	      y[0] = cod[17][1];
	    }

	  t1 = x[1] - x[0];
	  t2 = x[2] - x[0];
	  t3 = y[1] - y[0];
	  t4 = y[2] - y[0];
	  detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<np; k++)
	    {
	      gpx = tpw[k][0];
	      gpy = tpw[k][1];
	      xx = x[0] + t1*gpx + t2*gpy;
	      yy = y[0] + t3*gpx + t4*gpy;

	      double bcoord[2];
	      bcoord[0] = xx;
	      bcoord[1] = yy;
	      f(3 + l/4) += tpw[k][2] * ff->Eval(bcoord) * detJ;
	    } 
	}

      for(int l=0; l<6; l++)
	{
	  x[0] = cod[18][0];
	  y[0] = cod[18][1];
	  x[1] = cod[12+l][0];
	  y[1] = cod[12+l][1];
	  x[2] = cod[12+(l+1)%6][0];
	  y[2] = cod[12+(l+1)%6][1];
	  
	  t1 = x[1] - x[0];
	  t2 = x[2] - x[0];
	  t3 = y[1] - y[0];
	  t4 = y[2] - y[0];
	  detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<np; k++)
	    {
	      gpx = tpw[k][0];
	      gpy = tpw[k][1];
	      xx = x[0] + t1*gpx + t2*gpy;
	      yy = y[0] + t3*gpx + t4*gpy;

	      double bcoord[2];
	      bcoord[0] = xx;
	      bcoord[1] = yy;
	      f(9) += tpw[k][2] * ff->Eval(bcoord) * detJ;
	    } 
	}

      Vector RHS(nl);
      for(int k=0; k<nl; k++) RHS(k) = Q(k) - F[k] + f(k);

      if(i%2==0)
	{
	  Array<int> eind;
	  mesh->GetEdgeElements(edges[0], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[1]][7] - tflux[i][0] + tflux[eind[1]][7] );
		  RHS(3) += 0.5 * ( sflux[i][1] - sflux[eind[1]][6] - tflux[i][1] + tflux[eind[1]][6] ); 
		  RHS(4) += 0.5 * ( sflux[i][2] - sflux[eind[1]][5] - tflux[i][2] + tflux[eind[1]][5] );
		  RHS(1) += 0.5 * ( sflux[i][3] - sflux[eind[1]][4] - tflux[i][3] + tflux[eind[1]][4] ); 
		}
	      else
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[0]][7] - tflux[i][0] + tflux[eind[0]][7] );
		  RHS(3) += 0.5 * ( sflux[i][1] - sflux[eind[0]][6] - tflux[i][1] + tflux[eind[0]][6] ); 
		  RHS(4) += 0.5 * ( sflux[i][2] - sflux[eind[0]][5] - tflux[i][2] + tflux[eind[0]][5] );
		  RHS(1) += 0.5 * ( sflux[i][3] - sflux[eind[0]][4] - tflux[i][3] + tflux[eind[0]][4] ); 
		}
	    }
	  else
	    {
	      RHS(0) += -tflux[i][0] + sflux[i][0];
	      RHS(3) += -tflux[i][1] + sflux[i][1];
	      RHS(4) += -tflux[i][2] + sflux[i][2];
	      RHS(1) += -tflux[i][3] + sflux[i][3];
	    }

	  mesh->GetEdgeElements(edges[1], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(1) += 0.5 * ( sflux[i][4] - sflux[eind[1]][11] - tflux[i][4] + tflux[eind[1]][11] );
		  RHS(5) += 0.5 * ( sflux[i][5] - sflux[eind[1]][10] - tflux[i][5] + tflux[eind[1]][10] ); 
		  RHS(6) += 0.5 * ( sflux[i][6] - sflux[eind[1]][9] - tflux[i][6] + tflux[eind[1]][9] );
		  RHS(2) += 0.5 * ( sflux[i][7] - sflux[eind[1]][8] - tflux[i][7] + tflux[eind[1]][8] ); 
		}
	      else
		{
		  RHS(1) += 0.5 * ( sflux[i][4] - sflux[eind[0]][11] - tflux[i][4] + tflux[eind[0]][11] );
		  RHS(5) += 0.5 * ( sflux[i][5] - sflux[eind[0]][10] - tflux[i][5] + tflux[eind[0]][10] ); 
		  RHS(6) += 0.5 * ( sflux[i][6] - sflux[eind[0]][9] - tflux[i][6] + tflux[eind[0]][9] );
		  RHS(2) += 0.5 * ( sflux[i][7] - sflux[eind[0]][8] - tflux[i][7] + tflux[eind[0]][8] ); 
		}
	    }
	  else
	    {
	      RHS(1) += -tflux[i][4] + sflux[i][4];
	      RHS(5) += -tflux[i][5] + sflux[i][5];
	      RHS(6) += -tflux[i][6] + sflux[i][6];
	      RHS(2) += -tflux[i][7] + sflux[i][7];
	    }

	  mesh->GetEdgeElements(edges[2], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(2) += 0.5 * ( sflux[i][8] - sflux[eind[1]][3] - tflux[i][8] + tflux[eind[1]][3] );
		  RHS(7) += 0.5 * ( sflux[i][9] - sflux[eind[1]][2] - tflux[i][9] + tflux[eind[1]][2] ); 
		  RHS(8) += 0.5 * ( sflux[i][10] - sflux[eind[1]][1] - tflux[i][10] + tflux[eind[1]][1] );
		  RHS(0) += 0.5 * ( sflux[i][11] - sflux[eind[1]][0] - tflux[i][11] + tflux[eind[1]][0] ); 
		}
	      else
		{
		  RHS(2) += 0.5 * ( sflux[i][8] - sflux[eind[0]][3] - tflux[i][8] + tflux[eind[0]][3] );
		  RHS(7) += 0.5 * ( sflux[i][9] - sflux[eind[0]][2] - tflux[i][9] + tflux[eind[0]][2] ); 
		  RHS(8) += 0.5 * ( sflux[i][10] - sflux[eind[0]][1] - tflux[i][10] + tflux[eind[0]][1] );
		  RHS(0) += 0.5 * ( sflux[i][11] - sflux[eind[0]][0] - tflux[i][11] + tflux[eind[0]][0] ); 
		}
	    }
	  else
	    {
	      RHS(2) += -tflux[i][8] + sflux[i][8];
	      RHS(7) += -tflux[i][9] + sflux[i][9];
	      RHS(8) += -tflux[i][10] + sflux[i][10];
	      RHS(0) += -tflux[i][11] + sflux[i][11];
	    }
	}
      else
	{
	  Array<int> eind;
	  mesh->GetEdgeElements(edges[0], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[1]][11] - tflux[i][0] + tflux[eind[1]][11] );
		  RHS(3) += 0.5 * ( sflux[i][1] - sflux[eind[1]][10] - tflux[i][1] + tflux[eind[1]][10] ); 
		  RHS(4) += 0.5 * ( sflux[i][2] - sflux[eind[1]][9] - tflux[i][2] + tflux[eind[1]][9] );
		  RHS(1) += 0.5 * ( sflux[i][3] - sflux[eind[1]][8] - tflux[i][3] + tflux[eind[1]][8] ); 
		}
	      else
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[0]][11] - tflux[i][0] + tflux[eind[0]][11] );
		  RHS(3) += 0.5 * ( sflux[i][1] - sflux[eind[0]][10] - tflux[i][1] + tflux[eind[0]][10] ); 
		  RHS(4) += 0.5 * ( sflux[i][2] - sflux[eind[0]][9] - tflux[i][2] + tflux[eind[0]][9] );
		  RHS(1) += 0.5 * ( sflux[i][3] - sflux[eind[0]][8] - tflux[i][3] + tflux[eind[0]][8] ); 
		}
	    }
	  else
	    {
	      RHS(0) += -tflux[i][0] + sflux[i][0];
	      RHS(3) += -tflux[i][1] + sflux[i][1];
	      RHS(4) += -tflux[i][2] + sflux[i][2];
	      RHS(1) += -tflux[i][3] + sflux[i][3];
	    }

	  mesh->GetEdgeElements(edges[1], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(1) += 0.5 * ( sflux[i][4] - sflux[eind[1]][3] - tflux[i][4] + tflux[eind[1]][3] );
		  RHS(5) += 0.5 * ( sflux[i][5] - sflux[eind[1]][2] - tflux[i][5] + tflux[eind[1]][2] ); 
		  RHS(6) += 0.5 * ( sflux[i][6] - sflux[eind[1]][1] - tflux[i][6] + tflux[eind[1]][1] );
		  RHS(2) += 0.5 * ( sflux[i][7] - sflux[eind[1]][0] - tflux[i][7] + tflux[eind[1]][0] ); 
		}
	      else
		{
		  RHS(1) += 0.5 * ( sflux[i][4] - sflux[eind[0]][3] - tflux[i][4] + tflux[eind[0]][3] );
		  RHS(5) += 0.5 * ( sflux[i][5] - sflux[eind[0]][2] - tflux[i][5] + tflux[eind[0]][2] ); 
		  RHS(6) += 0.5 * ( sflux[i][6] - sflux[eind[0]][1] - tflux[i][6] + tflux[eind[0]][1] );
		  RHS(2) += 0.5 * ( sflux[i][7] - sflux[eind[0]][0] - tflux[i][7] + tflux[eind[0]][0] ); 
		}
	    }
	  else
	    {
	      RHS(1) += -tflux[i][4] + sflux[i][4];
	      RHS(5) += -tflux[i][5] + sflux[i][5];
	      RHS(6) += -tflux[i][6] + sflux[i][6];
	      RHS(2) += -tflux[i][7] + sflux[i][7];
	    }

	  mesh->GetEdgeElements(edges[2], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(2) += 0.5 * ( sflux[i][8] - sflux[eind[1]][7] - tflux[i][8] + tflux[eind[1]][7] );
		  RHS(7) += 0.5 * ( sflux[i][9] - sflux[eind[1]][6] - tflux[i][9] + tflux[eind[1]][6] ); 
		  RHS(8) += 0.5 * ( sflux[i][10] - sflux[eind[1]][5] - tflux[i][10] + tflux[eind[1]][5] );
		  RHS(0) += 0.5 * ( sflux[i][11] - sflux[eind[1]][4] - tflux[i][11] + tflux[eind[1]][4] ); 
		}
	      else
		{
		  RHS(2) += 0.5 * ( sflux[i][8] - sflux[eind[0]][7] - tflux[i][8] + tflux[eind[0]][7] );
		  RHS(7) += 0.5 * ( sflux[i][9] - sflux[eind[0]][6] - tflux[i][9] + tflux[eind[0]][6] ); 
		  RHS(8) += 0.5 * ( sflux[i][10] - sflux[eind[0]][5] - tflux[i][10] + tflux[eind[0]][5] );
		  RHS(0) += 0.5 * ( sflux[i][11] - sflux[eind[0]][4] - tflux[i][11] + tflux[eind[0]][4] ); 
		}
	    }
	  else
	    {
	      RHS(2) += -tflux[i][8] + sflux[i][8];
	      RHS(7) += -tflux[i][9] + sflux[i][9];
	      RHS(8) += -tflux[i][10] + sflux[i][10];
	      RHS(0) += -tflux[i][11] + sflux[i][11];
	    }
	}

      Vector s(nl);
      s = 0.0;
      
      RectangularMatrix localconsmat(nl);
      
      for(int j=0; j<nl; j++)
	for(int k=0; k<nl; k++)
	  localconsmat(j,k) = locmat->Elem(j,k);

      for(int j=0; j<nl; j++) localconsmat(9, j) = 0.0;
      localconsmat(9,9) = 1.0;
      RHS(9) = locsol(9); 

      DenseMatrixInverse invmat(localconsmat);
      invmat.Mult(RHS, s);
      //s.Print();

      Vector sss(10);
      sss = 0.0;
      locmat->Mult(s, sss);
      for(int k=0; k<10; k++) cverr(ind[k]) += sss(k) - f(k);
      locmat->Mult(locsol, sss);
      for(int k=0; k<10; k++) cverruh(ind[k]) += sss(k) - f(k);
      
      for(int k=0; k<nl; k++)
	{
	  ppsol[k][i] = s(k);
	}

      for(int k=0; k<18; k++)
	{
	  for(int j=0; j<nl; j++)
	    {
	      ppflux[k][i] += -ctem[k][j] * s(j);
	      flux[k][i] += -ctem[k][j] * locsol(j);
	    }
	}

      for(int k=0; k<femat.Size(); k++) delete []femat[k];
      for(int k=0; k<jit.Size(); k++) delete []jit[k];
      for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
      for (int k=0; k<ctem.Size(); k++) { delete []ctem[k]; }
      delete locmat;
    }

  for(int i=0; i<mesh->GetNBE(); i++)
    {
      Array<int> bdrind(2);
      mesh->GetBdrElementVertices(i, bdrind);

      if(bdrind[0] < bdrind[1])
	{
	  bdrind.Append( mesh->GetNV() + 2 * mesh->GetBdrElementEdgeIndex(i) );
	  bdrind.Append( mesh->GetNV() + 2 * mesh->GetBdrElementEdgeIndex(i) + 1);
	}
      else
	{
	  bdrind.Append( mesh->GetNV() + 2 * mesh->GetBdrElementEdgeIndex(i) + 1);
	  bdrind.Append( mesh->GetNV() + 2 * mesh->GetBdrElementEdgeIndex(i) );
	}

      if( fabs(cverr(bdrind[0])) > 1.0e-15 )
	{
	  ppbdrflux[bdrind[0]] = -cverr(bdrind[0]);
	  cverr(bdrind[0]) = 0.0;
	  cverruh(bdrind[0]) = 0.0;
	}
      if( fabs(cverr(bdrind[1])) > 1.0e-15 )
	{
	  ppbdrflux[bdrind[1]] = -cverr(bdrind[1]);
	  cverr(bdrind[1]) = 0.0;
	  cverruh(bdrind[1]) = 0.0;
	}
      
      ppbdrflux[bdrind[2]] = -cverr(bdrind[2]);
      cverr(bdrind[2]) = 0.0;
      cverruh(bdrind[2]) = 0.0;
      ppbdrflux[bdrind[3]] = -cverr(bdrind[3]);
      cverr(bdrind[3]) = 0.0;
      cverruh(bdrind[3]) = 0.0;
    }

  ofstream fileout("clce.out");
  for(int i=0; i<dualmesh->GetDualMeshNumDOF(); i++)
    {
      fileout<<i<<"\t"<<cverruh(i) << "\t" << cverr(i) <<endl;
      cout<<i<<"\t"<<cverruh(i) << "\t"<<cverr(i)<<endl;
    }
  fileout.close(); 
  
  for(int i=0; i<tflux.Size(); i++) delete []tflux[i];
  for(int i=0; i<sflux.Size(); i++) delete []sflux[i];
  for(int i=0; i<pw.Size(); i++) delete []pw[i];
  for(int i=0; i<tpw.Size(); i++) delete []tpw[i];
  for (int k=0; k<cod.Size(); k++) { delete []cod[k]; }
  for (int k=0; k<nmls.Size(); k++) { delete []nmls[k]; }
}

//===========================================================

Postprocessing::~Postprocessing()
{
  ;
}
