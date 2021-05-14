#include "mpi_fem_header.h"

//==============================================================================

QuadraticFEM2D_p::QuadraticFEM2D_p(Mesh_p *_mesh, Data *_data) :
FEM_p(_mesh, _data)
{
    strcpy(Type, "QuadraticFEM2D_p");
    NumGlobalDOF = NumVertices;
    NumLocalDOF = 6;
    rhs.SetSize(NumGlobalDOF);
    mat = new SparseMatrix(NumGlobalDOF, 10);
}

//============================================================================== 

void QuadraticFEM2D_p::GetElementDOF(int i, Array<int> &ind)
{
  GetElementVertices(i, ind);
}

//============================================================================== 

void QuadraticFEM2D_p::GetBdrElementDOF(int i, Array<int> &ind)
{
  GetBdrElementVertices(i, ind);
}

//==============================================================================

void QuadraticFEM2D_p::ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
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
  
  ///------------------------------------ row 1 -----------------------------------------------------------------
  locmat[0][0] = t1 * ( 2.0*k[0]/15.0 - k[1]/45.0 - k[2]/45.0 + 7.0*k[3]/90.0 + k[4]/6.0 + k[5]/6.0 );  
  locmat[0][1] = t2 * ( k[0]/30.0 + k[1]/30.0 - k[2]/45.0 + k[3]/18.0 + k[4]/18.0 + k[5]/90.0 );  
  locmat[0][2] = t3 * ( k[0]/30.0 + k[2]/30.0 - k[1]/45.0 + k[3]/18.0 + k[5]/18.0 + k[4]/90.0 ); 

  double s1 = 4.0 * t2 * ( -k[0]/120.0 + k[1]/360.0 + k[2]/60.0 + k[3]/45.0 - k[4]/45.0 - k[5]/90.0 );
  double s2 = 4.0 * t3 * ( -k[0]/120.0 + k[1]/60.0 + k[2]/360.0 + k[3]/45.0 - k[4]/90.0 - k[5]/45.0 );

  locmat[0][3] = s1 + s2;  
  locmat[0][4] = -s1 - 4.0 * t3 * ( k[0]/24.0 - k[1]/90.0 + k[2]/360.0 + k[3]/30.0 + k[4]*2.0/45.0 + k[5]/18.0 );  
  locmat[0][5] = -4.0 * t2 * ( k[0]/24.0 - k[2]/90.0 + k[1]/360.0 + k[3]/30.0 + k[5]*2.0/45.0 + k[4]/18.0 ) - s2;  

  ///------------------------------------ row 2 -----------------------------------------------------------------
  locmat[1][0] = locmat[0][1];  
  locmat[1][1] = t4 * ( 2.0*k[1]/15.0 - k[0]/45.0 - k[2]/45.0 + 7.0*k[4]/90.0 + k[3]/6.0 + k[5]/6.0 );  
  locmat[1][2] = t5 * ( -k[2]/30.0 - k[1]/30.0 + k[0]/45.0 - k[5]/18.0 - k[4]/18.0 - k[3]/90.0 );

  s1 = 4.0 * t4 * ( -k[0]/360.0 + k[1]/120.0 - k[2]/60.0 + k[3]/45.0 - k[4]/45.0 + k[5]/90.0 );
  s2 = 4.0 * t5 * ( -k[0]/72.0 + k[1]/20.0 - k[2]/72.0 + k[3]/15.0 + k[4]/90.0 + k[5]/15.0 );
 
  locmat[1][3] = s1 + s2;  
  locmat[1][4] = -s1 - 4.0 * t5 * ( k[0]/72.0 - k[2]/72.0 + k[3]/90.0 - k[5]/90.0 );  
  locmat[1][5] = -4.0 * t4 * ( k[1]/24.0 - k[2]/90.0 + k[0]/360.0 + k[4]/30.0 + k[5]*2.0/45.0 + k[3]/18.0 ) - s2;  

  ///------------------------------------ row 3 -----------------------------------------------------------------
  locmat[2][0] = locmat[0][2];  
  locmat[2][1] = locmat[1][2];  
  locmat[2][2] = t6 * ( 2.0*k[2]/15.0 - k[0]/45.0 - k[1]/45.0 + 7.0*k[5]/90.0 + k[3]/6.0 + k[3]/6.0 );  

  s1 = 4.0 * t7 * ( -k[0]/72.0 + k[2]/20.0 - k[1]/72.0 + k[3]/15.0 + k[5]/90.0 + k[4]/15.0 );  
  s2 = 4.0 * t6 * ( -k[0]/360.0 + k[2]/120.0 - k[1]/60.0 + k[3]/45.0 - k[5]/45.0 + k[4]/90.0 );

  locmat[2][3] = s1 + s2;  
  locmat[2][4] = -s1 - 4.0 * t6 * ( k[0]/360.0 - k[1]/90.0 + k[2]/24.0 + k[3]/18.0 + k[4]*2.0/45.0 + k[5]/30.0 ) ;  
  locmat[2][5] = -4.0 * t7 * ( k[0]/72.0 - k[1]/72.0 + k[3]/90.0 - k[4]/90.0 ) - s2;  

  ///------------------------------------ row 4 -----------------------------------------------------------------
  locmat[3][0] = locmat[0][3];  
  locmat[3][1] = locmat[1][3];  
  locmat[3][2] = locmat[2][3];

  s1 = 16.0 * t4 * ( -k[0]/180.0 - k[1]/180.0 + k[2]/60.0 + k[3]/30.0 + k[4]/30.0 + k[5]/90.0 );
  s2 = 16.0 * t6 * ( -k[0]/180.0 - k[2]/180.0 + k[1]/60.0 + k[3]/30.0 + k[5]/30.0 + k[4]/90.0 );
  
  locmat[3][3] = s1 + s2 + 32.0 * t5 * ( -k[0]/360.0 + k[3]/45.0 + k[4]/90.0 + k[5]/90.0 );

  double s3 = -16.0 * t6 * ( -k[0]/360.0 + k[2]/360.0 + k[3]/90.0 - k[5]/90.0 );
  double s4 = -16.0 * t4 * ( -k[0]/360.0 + k[1]/360.0 + k[3]/90.0 - k[4]/90.0 );

  locmat[3][4] = -s1 + s3 - 16.0 * t5 * ( -k[0]/120.0 - k[1]/360.0 + k[2]/60.0 + k[3]*2.0/45.0 + k[4]/45.0 + k[5]/90.0 );  
  locmat[3][5] = -s2 + s4 - 16.0 * t5 * ( -k[0]/120.0 - k[2]/360.0 + k[1]/60.0 + k[3]*2.0/45.0 + k[5]/45.0 + k[4]/90.0 );  

  ///------------------------------------ row 5 -----------------------------------------------------------------
  locmat[4][0] = locmat[0][4];  
  locmat[4][1] = locmat[1][4];  
  locmat[4][2] = locmat[2][4];  
  locmat[4][3] = locmat[3][4];  
  locmat[4][4] = s1 + 16.0 * t6 * ( k[0]/90.0 - k[1]/180.0 + k[2]/90.0 + k[3]/45.0 + k[4]/45.0 + k[5]/45.0 ) + 32.0 * t5 * ( -k[0]/180.0 - k[1]/360.0 + k[2]/60.0 + k[3]/45.0 + k[4]/90.0 );  
  locmat[4][5] = -s3 - s4 + 16.0 * t5 * ( k[0]/90.0 - k[1]/360.0 - k[2]/360.0 + k[3]/30.0 + k[4]/45.0 + k[5]/45.0 );  

  ///------------------------------------ row 6 -----------------------------------------------------------------
  locmat[5][0] = locmat[0][5];  
  locmat[5][1] = locmat[1][5];  
  locmat[5][2] = locmat[2][5];  
  locmat[5][3] = locmat[3][5];  
  locmat[5][4] = locmat[4][5];  
  locmat[5][5] = s2 + 16.0 * t6 * ( k[0]/90.0 - k[2]/180.0 + k[1]/90.0 + k[3]/45.0 + k[4]/45.0 + k[5]/45.0 ) + 32.0 * t5 * ( -k[0]/180.0 - k[2]/360.0 + k[1]/60.0 + k[3]/45.0 + k[5]/90.0 );  
    
  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void QuadraticFEM2D_p::ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
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

void QuadraticFEM2D_p::ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
{
  Array<double> f(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    f[j] = data->GetNodalForceCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[1][2] - coord[1][0]) - (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);

  locrhs[0] = detJ * ( 6.0*f[0] - f[1] - f[2] - 4.0*f[3] ) / 360.0;
  locrhs[1] = detJ * ( -f[0] + 6.0*f[1] - f[2] - 4.0*f[4] ) / 360.0;
  locrhs[2] = detJ * ( -f[0] - f[1] + 6.0*f[2] - 4.0*f[5] ) / 360.0;
  locrhs[3] = detJ * ( -f[0] + 8.0*f[3] + 4.0*f[4] + 4.0*f[5] ) / 90.0;
  locrhs[4] = detJ * ( -f[1] + 4.0*f[3] + 8.0*f[4] + 4.0*f[5] ) / 90.0;
  locrhs[5] = detJ * ( -f[2] + 4.0*f[3] + 4.0*f[4] + 8.0*f[5] ) / 90.0;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

void QuadraticFEM2D_p::ProjectAFunction(double (*func)(Array<double> &), Vector &pval)
{
  Array<double> coord(2);
  for(int i=0; i<mesh->GetNV(); i++)
    {
      coord[0] = mesh->GetVertex(i, 0);
      coord[1] = mesh->GetVertex(i, 1);      
      pval(i) = func(coord);
    }
}
