#include "fem_header.h"

QuadraticFEM2D::QuadraticFEM2D(Mesh *_mesh, Data *_data) :
FEM(_mesh, _data)
{
  strcpy(Type, "QuadraticFEM2D");
  NumGlobalDOF = mesh->GetNV() + mesh->GetNEdges();
  NumLocalDOF = 6;
  rhs.SetSize(NumGlobalDOF);
  mat = new SparseMatrix(NumGlobalDOF, 40);
}

//============================================================================== 

void QuadraticFEM2D::GetElementDOF(int i, Array<int> &ind)
{
  GetElementVertices(i, ind);

  Array<int> edges;
  mesh->GetElementEdges(i, edges);
  for (int k = 0; k < edges.Size(); k++)
    ind.Append (mesh->GetNV() + edges[k]);
}

//============================================================================== 

void QuadraticFEM2D::GetBdrElementDOF(int i, Array<int> &ind)
{
  GetBdrElementVertices(i, ind);
  ind.Append(mesh->GetNV() +  mesh->GetBdrElementEdgeIndex(i));
  //for(int k=0; k<3; k++) cout<<ind[k] <<endl; cout<<endl;
}

//==============================================================================

void QuadraticFEM2D::ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  
  Array<double> k(ind.Size()); // k for diffusion coefficient
  for (int j=0; j<ind.Size(); j++)
    k[j] = data->GetNodalEllipticCoeff(ind[j]);

  Array<double *> c(mesh->GetDim()); // c for coordinate
  for (int j=0; j<c.Size(); j++)
    c[j] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, c);
  double detJ = (c[0][1] - c[0][0]) * (c[1][2] - c[1][0]) - (c[1][1] - c[1][0]) * (c[0][2] - c[0][0]);

  double t1 = (c[1][2] - c[1][1]) * (c[1][2] - c[1][1]) + (c[0][1] - c[0][2])*(c[0][1] - c[0][2]);
  double t2 = (c[1][2] - c[1][1]) * (c[1][2] - c[1][0]) + (c[0][1] - c[0][2])*(c[0][0] - c[0][2]);
  double t3 = (c[1][2] - c[1][1]) * (c[1][0] - c[1][1]) + (c[0][1] - c[0][2])*(c[0][1] - c[0][0]);

  double t4 = (c[1][2] - c[1][0]) * (c[1][2] - c[1][0]) + (c[0][0] - c[0][2])*(c[0][0] - c[0][2]);
  double t5 = (c[1][2] - c[1][0]) * (c[1][0] - c[1][1]) + (c[0][0] - c[0][2])*(c[0][1] - c[0][0]);

  double t6 = (c[1][0] - c[1][1]) * (c[1][0] - c[1][1]) + (c[0][1] - c[0][0])*(c[0][1] - c[0][0]);
  double t7 = (c[1][0] - c[1][1]) * (c[1][2] - c[1][0]) + (c[0][1] - c[0][0])*(c[0][0] - c[0][2]);

  t1 /= detJ;
  t2 /= detJ;
  t3 /= detJ;
  t4 /= detJ;
  t5 /= detJ;
  t6 /= detJ;
  t7 /= detJ;
  
  ///------------------------------------ row 0 -----------------------------------------------------------------
  locmat[0][0] = t1 * ( 2.0*k[0]/15.0 - k[1]/45.0 - k[2]/45.0 + k[3]/6.0 + 7.0*k[4]/90.0 + k[5]/6.0 );  
  locmat[0][1] = t2 * ( k[0]/30.0 + k[1]/30.0 - k[2]/45.0 + k[3]/90.0 + k[4]/18.0 + k[5]/18.0 );  
  locmat[0][2] = t3 * ( k[0]/30.0 + k[2]/30.0 - k[1]/45.0 + k[3]/18.0 + k[4]/18.0 + k[5]/90.0 ); 

  double s1 = 4.0 * t3 * ( -k[0]/120.0 + k[1]/60.0 + k[2]/360.0 - k[3]/45.0 + k[4]/45.0 - k[5]/90.0 );
  double s2 = 4.0 * t2 * ( -k[0]/120.0 + k[1]/360.0 + k[2]/60.0 - k[3]/90.0 + k[4]/45.0 - k[5]/45.0 );

  locmat[0][3] = -s1 - 4.0 * t2 * ( k[0]/24.0 - k[2]/90.0 + k[1]/360.0 + k[4]/30.0 + k[3]*2.0/45.0 + k[5]/18.0 );
  locmat[0][4] = s1 + s2;
  locmat[0][5] = -s2 - 4.0 * t3 * ( k[0]/24.0 - k[1]/90.0 + k[2]/360.0 + k[4]/30.0 + k[5]*2.0/45.0 + k[3]/18.0 );

  ///------------------------------------ row 1 -----------------------------------------------------------------
  locmat[1][0] = locmat[0][1];  
  locmat[1][1] = t4 * ( 2.0*k[1]/15.0 - k[0]/45.0 - k[2]/45.0 + 7.0*k[5]/90.0 + k[3]/6.0 + k[4]/6.0 );  
  locmat[1][2] = t5 * ( -k[2]/30.0 - k[1]/30.0 + k[0]/45.0 - k[3]/18.0 - k[5]/18.0 - k[4]/90.0 );

  s1 = 4.0 * t4 * ( -k[0]/360.0 + k[1]/120.0 - k[2]/60.0 + k[4]/45.0 - k[5]/45.0 + k[3]/90.0 );
  s2 = 4.0 * t5 * ( -k[0]/72.0 + k[1]/20.0 - k[2]/72.0 + k[3]/15.0 + k[5]/90.0 + k[4]/15.0 );
 
  locmat[1][3] = -s2 - 4.0 * t4 * ( k[1]/24.0 - k[2]/90.0 + k[0]/360.0 + k[5]/30.0 + k[3]*2.0/45.0 + k[4]/18.0 );
  locmat[1][4] = s1 + s2;
  locmat[1][5] = -s1 - 4.0 * t5 * ( k[0]/72.0 - k[2]/72.0 - k[3]/90.0 + k[4]/90.0 ); 

  ///------------------------------------ row 2 -----------------------------------------------------------------
  locmat[2][0] = locmat[0][2];  
  locmat[2][1] = locmat[1][2];  
  locmat[2][2] = t6 * ( 2.0*k[2]/15.0 - k[0]/45.0 - k[1]/45.0 + 7.0*k[3]/90.0 + k[4]/6.0 + k[5]/6.0 );  

  s1 = 4.0 * t7 * ( -k[0]/72.0 + k[2]/20.0 - k[1]/72.0 + k[5]/15.0 + k[3]/90.0 + k[4]/15.0 );  
  s2 = 4.0 * t6 * ( -k[0]/360.0 + k[2]/120.0 - k[1]/60.0 + k[4]/45.0 - k[3]/45.0 + k[5]/90.0 );

  locmat[2][3] = -s2 - 4.0 * t7 * ( k[0]/72.0 - k[1]/72.0 + k[4]/90.0 - k[5]/90.0 );  
  locmat[2][4] = s1 + s2;  
  locmat[2][5] = -s1 - 4.0 * t6 * ( k[0]/360.0 - k[1]/90.0 + k[2]/24.0 + k[4]/18.0 + k[5]*2.0/45.0 + k[3]/30.0 );

  ///------------------------------------ row 3 -----------------------------------------------------------------
  locmat[3][0] = locmat[0][3];  
  locmat[3][1] = locmat[1][3];  
  locmat[3][2] = locmat[2][3];

  s1 = 16.0 * t4 * ( -k[0]/180.0 - k[1]/180.0 + k[2]/60.0 + k[5]/30.0 + k[4]/30.0 + k[3]/90.0 );
  s2 = 16.0 * t6 * ( -k[0]/180.0 - k[2]/180.0 + k[1]/60.0 + k[3]/30.0 + k[4]/30.0 + k[5]/90.0 );
  
  locmat[3][3] = s2 + 16.0 * t4 * ( k[0]/90.0 - k[2]/180.0 + k[1]/90.0 + k[3]/45.0 + k[4]/45.0 + k[5]/45.0 )
    + 32.0 * t5 * ( -k[0]/180.0 - k[2]/360.0 + k[1]/60.0 + k[4]/45.0 + k[3]/90.0 );  

  locmat[3][4] = -s2 - 16.0 * t4 * ( -k[0]/360.0 + k[1]/360.0 + k[4]/90.0 - k[5]/90.0 )
    - 16.0 * t5 * ( -k[0]/120.0 - k[2]/360.0 + k[1]/60.0 + k[4]*2.0/45.0 + k[3]/45.0 + k[5]/90.0 );  
  locmat[3][5] = 16.0 * t4 * ( -k[0]/360.0 + k[1]/360.0 + k[4]/90.0 - k[5]/90.0 )
    + 16.0 * t6 * ( -k[0]/360.0 + k[2]/360.0 - k[3]/90.0 + k[4]/90.0 )
    + 16.0 * t5 * ( k[0]/90.0 - k[1]/360.0 - k[2]/360.0 + k[3]/45.0 + k[5]/45.0 + k[4]/30.0 );  

  ///------------------------------------ row 4 -----------------------------------------------------------------
  locmat[4][0] = locmat[0][4];  
  locmat[4][1] = locmat[1][4];  
  locmat[4][2] = locmat[2][4];  
  locmat[4][3] = locmat[3][4];  
  locmat[4][4] = s1 + s2 + 32.0 * t5 * ( -k[0]/360.0 + k[4]/45.0 + k[3]/90.0 + k[5]/90.0 );

  locmat[4][5] = -s1 - 16.0 * t6 * ( -k[0]/360.0 + k[2]/360.0 - k[3]/90.0 + k[4]/90.0 )
    - 16.0 * t5 * ( -k[0]/120.0 - k[1]/360.0 + k[2]/60.0 + k[3]/90.0 + 2.0*k[4]/45.0 + k[5]/45.0 );  

  ///------------------------------------ row 5 -----------------------------------------------------------------
  locmat[5][0] = locmat[0][5];  
  locmat[5][1] = locmat[1][5];  
  locmat[5][2] = locmat[2][5];  
  locmat[5][3] = locmat[3][5];  
  locmat[5][4] = locmat[4][5];  
  locmat[5][5] = s1 + 16.0 * t6 * ( k[0]/90.0 - k[1]/180.0 + k[2]/90.0 + k[3]/45.0 + k[4]/45.0 + k[5]/45.0 )
    + 32.0 * t5 * ( -k[0]/180.0 - k[1]/360.0 + k[2]/60.0 + k[4]/45.0 + k[5]/90.0 );
    
  for (int k=0; k<c.Size(); k++) { delete []c[k]; }
  
  /*  
  for(int k=0; k<locmat.Size(); k++)
    for(int j=0; j<locmat.Size(); j++)
      locmat[k][j] = 0.0;
	  
  int np = ConvertN2PN(5);
  Array<double *> tpw(np);
  for(int j=0; j<np; j++) tpw[j] = new double[3];
  GetStandTriQuadPW(5, tpw);

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

  BasisFunctions *phi = new BasisFunctions(2);
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
  */
}

//==============================================================================

void QuadraticFEM2D::ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  cout << "Not Tested: Test before you use. " << endl;
  Array<double> r(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    r[j] = data->GetNodalReactionCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[1][2] - coord[1][0]) - (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);

  locmat[0][0] = detJ * ( 9.0*r[0] - r[1] - r[2] + 2.0*r[3] + 6.0*r[4] + 6.0*r[5] ) / 1260.0;  
  locmat[0][1] = detJ * ( -2.0*r[0] - 2.0*r[1] + r[2] - 4.0*r[5] ) / 2520.0;
  locmat[0][2] = detJ * ( -2.0*r[0] + r[1] - 2.0*r[2] - 4.0*r[4] ) / 2520.0;  
  locmat[0][3] = detJ * ( r[0] - 4.0*r[3] - 2.0*r[4] - 2.0*r[5] ) / 630.0;  
  locmat[0][4] = detJ * ( 3.0*r[0] - r[2] - 2.0*r[3] ) / 630.0;  
  locmat[0][5] = detJ * ( 3.0*r[0] - r[1] - 2.0*r[3] ) / 630.0;  

  locmat[1][0] = locmat[0][1];  
  locmat[1][1] = detJ * ( -r[0] + 9.0*r[1] - r[2] + 6.0*r[3] + 2.0*r[4] + 6.0*r[5] ) / 1260.0;  
  locmat[1][2] = detJ * ( r[0] - 2.0*r[1] - 2.0*r[2] - 4.0*r[3] ) / 2520.0;;  
  locmat[1][3] = detJ * (3.0*r[1] - r[2] - 2.0*r[4] ) / 630.0;  
  locmat[1][4] = detJ * ( r[1] - 2.0*r[3] - 4.0*r[4] - 2.0*r[5] ) / 630.0;  
  locmat[1][5] = detJ * (3.0*r[1] - r[0] - 2.0*r[4] ) / 630.0;  

  locmat[2][0] = locmat[0][2];  
  locmat[2][1] = locmat[1][2];  
  locmat[2][2] = detJ * ( -r[0] - r[1] + 9.0*r[2] + 6.0*r[3] + 6.0*r[4] + 2.0*r[5] ) / 1260.0;  
  locmat[2][3] = detJ * ( -r[1] + 3.0*r[2] - 2.0*r[5] ) / 630.0;  
  locmat[2][4] = detJ * ( -r[0] + 3.0*r[2] - 2.0*r[5] ) / 630.0;    
  locmat[2][5] = detJ * ( r[2] - 2.0*r[3] - 2.0*r[4] - 4.0*r[5] ) / 630.0;  

  locmat[3][0] = locmat[0][3];  
  locmat[3][1] = locmat[1][3];  
  locmat[3][2] = locmat[2][3];  
  locmat[3][3] = detJ * ( -2.0*r[0] + 9.0*r[3] + 6.0*r[4] + 6.0*r[5] ) / 315.0;  
  locmat[3][4] = detJ* ( -2.0*r[0] - 2.0*r[1] + 12.0*r[3] + 12.0*r[4] + 8.0*r[5] ) / 630.0;  
  locmat[3][5] = detJ* ( -2.0*r[0] - 2.0*r[2] + 12.0*r[3] + 8.0*r[4] + 12.0*r[5] ) / 630.0;  

  locmat[4][0] = locmat[0][4];  
  locmat[4][1] = locmat[1][4];  
  locmat[4][2] = locmat[2][4];  
  locmat[4][3] = locmat[3][4];  
  locmat[4][4] = detJ * ( -2.0*r[1] + 6.0*r[3] + 9.0*r[4] + 6.0*r[5] ) / 315.0;  
  locmat[4][5] = detJ* ( -2.0*r[1] - 2.0*r[2] + 8.0*r[3] + 12.0*r[4] + 12.0*r[5] ) / 630.0;  

  locmat[5][0] = locmat[0][5];  
  locmat[5][1] = locmat[1][5];  
  locmat[5][2] = locmat[2][5];  
  locmat[5][3] = locmat[3][5];  
  locmat[5][4] = locmat[4][5];  
  locmat[5][5] = detJ * ( -2.0*r[2] + 6.0*r[3] + 6.0*r[4] + 9.0*r[5] ) / 315.0; 

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void QuadraticFEM2D::ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
{
  
  Array<double> f(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    f[j] = data->GetNodalForceCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[1][2] - coord[1][0])
    - (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);

  locrhs[0] = detJ * ( 6.0*f[0] - f[1] - f[2] - 4.0*f[4] ) / 360.0;
  locrhs[1] = detJ * ( -f[0] + 6.0*f[1] - f[2] - 4.0*f[5] ) / 360.0;
  locrhs[2] = detJ * ( -f[0] - f[1] + 6.0*f[2] - 4.0*f[3] ) / 360.0;
  locrhs[3] = detJ * ( -f[2] + 8.0*f[3] + 4.0*f[4] + 4.0*f[5] ) / 90.0;
  locrhs[4] = detJ * ( -f[0] + 4.0*f[3] + 8.0*f[4] + 4.0*f[5] ) / 90.0;
  locrhs[5] = detJ * ( -f[1] + 4.0*f[3] + 4.0*f[4] + 8.0*f[5] ) / 90.0;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
  
  /*
  for(int j=0; j<locrhs.Size(); j++) locrhs[j] = 0.0;

  int np = ConvertN2PN(5);
  Array<double *> tpw(np);
  for(int j=0; j<np; j++) tpw[j] = new double[3];
  GetStandTriQuadPW(5, tpw);

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

  BasisFunctions *phi = new BasisFunctions(2);
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
  */
}

//==============================================================================

void QuadraticFEM2D::ComputeNeumannLocalSystem(int b, int i, Array<int> &ind, Array<double> &locrhs)
{
  Array<double> locdat(ind.Size());

  for (int j=0; j<ind.Size(); j++)
    {
      locdat[j] = data->GetNeumannBdryVal(b, ind[j]);
    }

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetBdrElement(i)->GetNVertices()];

  mesh->GetBdrElementVerticesCoord(i, coord);

  double h = (coord[0][1] - coord[0][0]) * (coord[0][1] - coord[0][0]) 
    + (coord[1][1] - coord[1][0]) * (coord[1][1] - coord[1][0]);

  h = sqrt(h);
  locrhs[0] = -h * ( 4.0*locdat[0] - locdat[1] + 2.0*locdat[2] ) / 30.0;
  locrhs[1] = -h * ( -locdat[0] + 4.0*locdat[1] + 2.0*locdat[2] ) / 30.0;
  locrhs[2] = -h * ( locdat[0] + 8.0*locdat[1] + locdat[2] ) / 15.0;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

void QuadraticFEM2D::ComputeNeumannLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalNeumannCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetBdrElement(i)->GetNVertices()];
  
  mesh->GetBdrElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[0][1] - coord[0][0]) 
                + (coord[1][1] - coord[1][0]) * (coord[1][1] - coord[1][0]);

  detJ = sqrt(detJ);

  locrhs[0] = -detJ * ( 4.0*locdat[0] - locdat[1] + 2.0*locdat[2] ) / 30.0;
  locrhs[1] = -detJ * ( -locdat[0] + 4.0*locdat[1] + 2.0*locdat[2] ) / 30.0;
  locrhs[2] = -detJ * ( locdat[0] + 8.0*locdat[1] + locdat[2] ) / 15.0;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void QuadraticFEM2D::ComputeRobinLocalSystem(int i, int b, Array<int> &ind, Array<double *> &locmat, Array<double> &locrhs)
{
  cout << "Not Implemented" << endl;
}

void QuadraticFEM2D::ComputeRobinLocalSystem(int i, Array<int> &ind, Array<double *> &locmat, Array<double> &locrhs)
{
  cout << "Not Implemented" << endl;
}

//==============================================================================

void QuadraticFEM2D::SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
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
  // mat->Print();
}

//==============================================================================

double QuadraticFEM2D::ComputeL2Error(Function *exactsol, Vector &femsol, double time)
{
  int N = 5; // accuracy of Gaussian Quadrature
  int np = ConvertN2PN(N);
  Array<double *> pw(np);
  for(int i=0; i<np; i++) pw[i] = new double[3];
  GetStandTriQuadPW(N, pw);  
  
  BasisFunctions *phi = new BasisFunctions(2);
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
	  for(int j=0; j<6; j++) uh += femsol(ind[j]) * phi->BF2D(j, xx, yy);
	  uh -= exactsol->Eval(coord);
	  err += pw[k][2] * uh * uh * detJ;
	}
    }  
  return sqrt(err);
}

//==============================================================================

void QuadraticFEM2D::ComputePPL2Error(Array<double> &errs, Function *exactsol, Vector &femsol, Array<double *> &ppsol, double time)
{
  errs[0] = 0.0;
  errs[1] = 0.0;
  errs[2] = 0.0;
  
  int N = 5; // accuracy of Gaussian Quadrature
  int np = ConvertN2PN(N);
  Array<double *> pw(np);
  for(int i=0; i<np; i++) pw[i] = new double[3];
  GetStandTriQuadPW(N, pw);  
  
  BasisFunctions *phi = new BasisFunctions(2);
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
	  for(int j=0; j<6; j++) femuh += femsol(ind[j]) * phi->BF2D(j, xx, yy);
	  double ppuh = 0.0;
	  for(int j=0; j<6; j++) ppuh += ppsol[j][i] * phi->BF2D(j, xx, yy);
	  
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

double QuadraticFEM2D::ComputeH1Error(Array<Function *> &exactgrad, Vector &femsol)
{
  int N = 5; // accuracy of Gaussian Quadrature
  int np = ConvertN2PN(N);
  Array<double *> pw(np);
  for(int i=0; i<np; i++) pw[i] = new double[3];
  GetStandTriQuadPW(N, pw);  
  
  Array<double *> jit(2);
  jit[0] = new double[2];
  jit[1] = new double[2];

  double xx, yy, err = 0.0;
  BasisFunctions *phi = new BasisFunctions(2);
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
	  for(int j=0; j<6; j++)
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

void QuadraticFEM2D::ComputePPH1Error(Array<double> &errs, Array<Function *> &exactgrad, Vector &femsol, Array<double *> &ppsol)
{
  errs[0] = 0.0;
  errs[1] = 0.0;
  errs[2] = 0.0;
  
  int N = 5; // accuracy of Gaussian Quadrature
  int np = ConvertN2PN(N);
  Array<double *> pw(np);
  for(int i=0; i<np; i++) pw[i] = new double[3];
  GetStandTriQuadPW(N, pw);  
  
  Array<double *> jit(2);
  jit[0] = new double[2];
  jit[1] = new double[2];

  double xx, yy;
  BasisFunctions *phi = new BasisFunctions(2);
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
	  for(int j=0; j<6; j++)
	    {
	      femuhx += femsol(ind[j]) * ( jit[0][0] * phi->GradBF2D(j, 0, xx, yy) + jit[0][1] * phi->GradBF2D(j, 1, xx, yy) );
	      femuhy += femsol(ind[j]) * ( jit[1][0] * phi->GradBF2D(j, 0, xx, yy) + jit[1][1] * phi->GradBF2D(j, 1, xx, yy) );
	    }

	  double ppuhx = 0.0;
	  double ppuhy = 0.0;
	  for(int j=0; j<6; j++)
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
