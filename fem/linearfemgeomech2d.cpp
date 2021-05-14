#include "fem_header.h"

//==============================================================================

LinearFEMGeoMech2D::LinearFEMGeoMech2D(Mesh *_mesh, Data *_data) :
FEM(_mesh, _data)
{
  strcpy(Type, "LinearFEMGeoMech2D");
  NumGlobalDOF = 3*mesh->GetNV();
  NumLocalDOF = 9;
  rhs.SetSize(NumGlobalDOF);
  mat = new SparseMatrix(NumGlobalDOF, 30);
}

//============================================================================= 

void LinearFEMGeoMech2D::GetElementDOF(int i, Array<int> &ind)
{
  Array<int> tempind;
  GetElementVertices(i, tempind);
  ind.SetSize(9);  
  for(int i=0; i<3; i++)
    {
      ind[i] = tempind[i];
      ind[i+3] = tempind[i] + mesh->GetNV();
      ind[i+6] = tempind[i] + 2 * mesh->GetNV();
    }
}

//==============================================================================

void LinearFEMGeoMech2D::GetBdrElementDOF(int i, Array<int> &ind)
{
  Array<int> tempind;
  GetBdrElementVertices(i, tempind);
  ind.SetSize(6);
  for(int i=0; i<2; i++)
    {
      ind[i] = tempind[i];
      ind[i+2] = tempind[i] + mesh->GetNV();
      ind[i+4] = tempind[i] + 2 * mesh->GetNV();
    } 
}

//==============================================================================

void LinearFEMGeoMech2D::ComputeEllipticLocalSystem(int i, Array<int> &ind, 
						  Array<double *> &locmat)
{
 Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];

  mesh->GetElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[1][2] 
		   - coord[1][0]) - (coord[1][1] - coord[1][0]) * 
		    (coord[0][2] - coord[0][0]);
  
  double intval = 1.0/(6.0*detJ);
 
  double musum, lambdasum;
  lambdasum = 0.0;
  musum = 0.0;
  for(int j=0; j<3; j++)
    {
      lambdasum += data->GetNodalEllipticCoeff(ind[j], 0);
      musum += data->GetNodalEllipticCoeff(ind[j], 1);
    }

  double dp0x, dp0y, dp1x, dp1y, dp2x, dp2y;
  dp0x = (coord[1][1]-coord[1][2]);
  dp0y = (coord[0][2]-coord[0][1]);
  dp1x = (coord[1][2]-coord[1][0]);
  dp1y = (coord[0][0]-coord[0][2]);
  dp2x = (coord[1][0]-coord[1][1]);
  dp2y = (coord[0][1]-coord[0][0]);

  locmat[0][0] = (musum + lambdasum)*dp0x*dp0x + 0.5*musum*dp0y*dp0y;
  locmat[0][1] = (musum + lambdasum)*dp1x*dp0x + 0.5*musum*dp1y*dp0y;
  locmat[0][2] = (musum + lambdasum)*dp2x*dp0x + 0.5*musum*dp2y*dp0y;
  locmat[0][3] = 0.5*musum*dp0x*dp0y + lambdasum*dp0y*dp0x;
  locmat[0][4] = 0.5*musum*dp1x*dp0y + lambdasum*dp1y*dp0x;
  locmat[0][5] = 0.5*musum*dp2x*dp0y + lambdasum*dp2y*dp0x;
  
  locmat[1][0] = (musum + lambdasum)*dp0x*dp1x + 0.5*musum*dp0y*dp1y;
  locmat[1][1] = (musum + lambdasum)*dp1x*dp1x + 0.5*musum*dp1y*dp1y;
  locmat[1][2] = (musum + lambdasum)*dp2x*dp1x + 0.5*musum*dp2y*dp1y;
  locmat[1][3] = 0.5*musum*dp0x*dp1y + lambdasum*dp0y*dp1x;
  locmat[1][4] = 0.5*musum*dp1x*dp1y + lambdasum*dp1y*dp1x;
  locmat[1][5] = 0.5*musum*dp2x*dp1y + lambdasum*dp2y*dp1x;
  
  locmat[2][0] = (musum + lambdasum)*dp0x*dp2x + 0.5*musum*dp0y*dp2y;
  locmat[2][1] = (musum + lambdasum)*dp1x*dp2x + 0.5*musum*dp1y*dp2y;
  locmat[2][2] = (musum + lambdasum)*dp2x*dp2x + 0.5*musum*dp2y*dp2y;
  locmat[2][3] = 0.5*musum*dp0x*dp2y + lambdasum*dp0y*dp2x;
  locmat[2][4] = 0.5*musum*dp1x*dp2y + lambdasum*dp1y*dp2x;
  locmat[2][5] = 0.5*musum*dp2x*dp2y + lambdasum*dp2y*dp2x;
  
  locmat[3][0] = 0.5*musum*dp0y*dp0x + lambdasum*dp0x*dp0y;
  locmat[3][1] = 0.5*musum*dp1y*dp0x + lambdasum*dp1x*dp0y;
  locmat[3][2] = 0.5*musum*dp2y*dp0x + lambdasum*dp2x*dp0y;
  locmat[3][3] = (lambdasum + musum)*dp0y*dp0y + 0.5*musum*dp0x*dp0x;
  locmat[3][4] = (lambdasum + musum)*dp1y*dp0y + 0.5*musum*dp1x*dp0x;
  locmat[3][5] = (lambdasum + musum)*dp2y*dp0y + 0.5*musum*dp2x*dp0x;

  locmat[4][0] = 0.5*musum*dp0y*dp1x + lambdasum*dp0x*dp1y;
  locmat[4][1] = 0.5*musum*dp1y*dp1x + lambdasum*dp1x*dp1y;
  locmat[4][2] = 0.5*musum*dp2y*dp1x + lambdasum*dp2x*dp1y;
  locmat[4][3] = (lambdasum + musum)*dp0y*dp1y + 0.5*musum*dp0x*dp1x;
  locmat[4][4] = (lambdasum + musum)*dp1y*dp1y + 0.5*musum*dp1x*dp1x;
  locmat[4][5] = (lambdasum + musum)*dp2y*dp1y + 0.5*musum*dp2x*dp1x;

  locmat[5][0] = 0.5*musum*dp0y*dp2x + lambdasum*dp0x*dp2y;
  locmat[5][1] = 0.5*musum*dp1y*dp2x + lambdasum*dp1x*dp2y;
  locmat[5][2] = 0.5*musum*dp2y*dp2x + lambdasum*dp2x*dp2y;
  locmat[5][3] = (lambdasum + musum)*dp0y*dp2y + 0.5*musum*dp0x*dp2x;
  locmat[5][4] = (lambdasum + musum)*dp1y*dp2y + 0.5*musum*dp1x*dp2x;
  locmat[5][5] = (lambdasum + musum)*dp2y*dp2y + 0.5*musum*dp2x*dp2x;

  for(int j=0; j<6; j++)
    {
      for(int k=0; k<6; k++)
	{
	  locmat[k][j] *= intval;
	}
    }

  Array<double> locdat(3);
  for (int j=0; j<3; j++)
    locdat[j] = data->GetNodalEllipticCoeff(ind[j], 2);
  
  double t = ( locdat[0] + locdat[1] + locdat[2] ) / ( 6.0*detJ );

  locmat[6][6] = t * ( (coord[1][1] - coord[1][2]) * (coord[1][1] 
                 - coord[1][2]) + (coord[0][2] - coord[0][1]) * 
                 (coord[0][2] - coord[0][1]) );
  locmat[6][7] = t * ( (coord[1][1] - coord[1][2]) * (coord[1][2] 
                 - coord[1][0]) + (coord[0][2] - coord[0][1]) * 
                 (coord[0][0] - coord[0][2]) );
  locmat[6][8] = t * ( (coord[1][1] - coord[1][2]) * (coord[1][0] 
                 - coord[1][1]) + (coord[0][2] - coord[0][1]) * 
                 (coord[0][1] - coord[0][0]) );

  locmat[7][6] = locmat[6][7];
  locmat[7][7] = t * ( (coord[1][0] - coord[1][2]) * (coord[1][0] 
                 - coord[1][2]) + (coord[0][2] - coord[0][0]) * 
                 (coord[0][2] - coord[0][0]) );
  locmat[7][8] = t * ( (coord[1][2] - coord[1][0]) * (coord[1][0] 
                 - coord[1][1]) + (coord[0][0] - coord[0][2]) * 
                 (coord[0][1] - coord[0][0]) );

  locmat[8][6] = locmat[6][8];
  locmat[8][7] = locmat[7][8];
  locmat[8][8] =  t * ( (coord[1][0] - coord[1][1]) * (coord[1][0] 
                  - coord[1][1]) + (coord[0][1] - coord[0][0]) * 
                  (coord[0][1] - coord[0][0]) );

  Array<double *> locd(2);
  for (int k=0; k<locd.Size(); k++)
    locd[k] = new double[3];

  for (int j=0; j<3; j++)
    {
      locd[0][j] = 1.0;
      locd[1][j] = 0.0;
    }
  
  locmat[0][6] = ( coord[1][1] - coord[1][2] ) * ( 2.0*locd[0][0] + locd[0][1] + locd[0][2] ) / 24.0 + ( coord[0][2] - coord[0][1] ) * ( 2.0*locd[1][0] + locd[1][1] + locd[1][2] ) / 24.0;
  locmat[0][7] = ( coord[1][2] - coord[1][0] ) * ( 2.0*locd[0][0] + locd[0][1] + locd[0][2] ) / 24.0 + ( coord[0][0] - coord[0][2] ) * ( 2.0*locd[1][0] + locd[1][1] + locd[1][2] ) / 24.0;
  locmat[0][8] = ( coord[1][0] - coord[1][1] ) * ( 2.0*locd[0][0] + locd[0][1] + locd[0][2] ) / 24.0 + ( coord[0][1] - coord[0][0] ) * ( 2.0*locd[1][0] + locd[1][1] + locd[1][2] ) / 24.0;

  locmat[1][6] = ( coord[1][1] - coord[1][2] ) * ( locd[0][0] + 2.0*locd[0][1] + locd[0][2] ) / 24.0 + ( coord[0][2] - coord[0][1] ) * ( locd[1][0] + 2.0*locd[1][1] + locd[1][2] ) / 24.0;
  locmat[1][7] = ( coord[1][2] - coord[1][0] ) * ( locd[0][0] + 2.0*locd[0][1] + locd[0][2] ) / 24.0 + ( coord[0][0] - coord[0][2] ) * ( locd[1][0] + 2.0*locd[1][1] + locd[1][2] ) / 24.0;
  locmat[1][8] = ( coord[1][0] - coord[1][1] ) * ( locd[0][0] + 2.0*locd[0][1] + locd[0][2] ) / 24.0 + ( coord[0][1] - coord[0][0] ) * ( locd[1][0] + 2.0*locd[1][1] + locd[1][2] ) / 24.0;

  locmat[2][6] = ( coord[1][1] - coord[1][2] ) * ( locd[0][0] + locd[0][1] + 2.0*locd[0][2] ) / 24.0 + ( coord[0][2] - coord[0][1] ) * ( locd[1][0] + locd[1][1] + 2.0*locd[1][2] ) / 24.0;
  locmat[2][7] = ( coord[1][2] - coord[1][0] ) * ( locd[0][0] + locd[0][1] + 2.0*locd[0][2] ) / 24.0 + ( coord[0][0] - coord[0][2] ) * ( locd[1][0] + locd[1][1] + 2.0*locd[1][2] ) / 24.0;
  locmat[2][8] = ( coord[1][0] - coord[1][1] ) * ( locd[0][0] + locd[0][1] + 2.0*locd[0][2] ) / 24.0 + ( coord[0][1] - coord[0][0] ) * ( locd[1][0] + locd[1][1] + 2.0*locd[1][2] ) / 24.0;

  for (int j=0; j<3; j++)
    {
      locd[0][j] = 0.0;
      locd[1][j] = 1.0;
    }
  
  locmat[3][6] = ( coord[1][1] - coord[1][2] ) * ( 2.0*locd[0][0] + locd[0][1] + locd[0][2] ) / 24.0 + ( coord[0][2] - coord[0][1] ) * ( 2.0*locd[1][0] + locd[1][1] + locd[1][2] ) / 24.0;
  locmat[3][7] = ( coord[1][2] - coord[1][0] ) * ( 2.0*locd[0][0] + locd[0][1] + locd[0][2] ) / 24.0 + ( coord[0][0] - coord[0][2] ) * ( 2.0*locd[1][0] + locd[1][1] + locd[1][2] ) / 24.0;
  locmat[3][8] = ( coord[1][0] - coord[1][1] ) * ( 2.0*locd[0][0] + locd[0][1] + locd[0][2] ) / 24.0 + ( coord[0][1] - coord[0][0] ) * ( 2.0*locd[1][0] + locd[1][1] + locd[1][2] ) / 24.0;

  locmat[4][6] = ( coord[1][1] - coord[1][2] ) * ( locd[0][0] + 2.0*locd[0][1] + locd[0][2] ) / 24.0 + ( coord[0][2] - coord[0][1] ) * ( locd[1][0] + 2.0*locd[1][1] + locd[1][2] ) / 24.0;
  locmat[4][7] = ( coord[1][2] - coord[1][0] ) * ( locd[0][0] + 2.0*locd[0][1] + locd[0][2] ) / 24.0 + ( coord[0][0] - coord[0][2] ) * ( locd[1][0] + 2.0*locd[1][1] + locd[1][2] ) / 24.0;
  locmat[4][8] = ( coord[1][0] - coord[1][1] ) * ( locd[0][0] + 2.0*locd[0][1] + locd[0][2] ) / 24.0 + ( coord[0][1] - coord[0][0] ) * ( locd[1][0] + 2.0*locd[1][1] + locd[1][2] ) / 24.0;

  locmat[5][6] = ( coord[1][1] - coord[1][2] ) * ( locd[0][0] + locd[0][1] + 2.0*locd[0][2] ) / 24.0 + ( coord[0][2] - coord[0][1] ) * ( locd[1][0] + locd[1][1] + 2.0*locd[1][2] ) / 24.0;
  locmat[5][7] = ( coord[1][2] - coord[1][0] ) * ( locd[0][0] + locd[0][1] + 2.0*locd[0][2] ) / 24.0 + ( coord[0][0] - coord[0][2] ) * ( locd[1][0] + locd[1][1] + 2.0*locd[1][2] ) / 24.0;
  locmat[5][8] = ( coord[1][0] - coord[1][1] ) * ( locd[0][0] + locd[0][1] + 2.0*locd[0][2] ) / 24.0 + ( coord[0][1] - coord[0][0] ) * ( locd[1][0] + locd[1][1] + 2.0*locd[1][2] ) / 24.0;

  //--------------------------------------------------
  double dt = data->GetNodalEllipticCoeff(mesh->GetNV(), 3);
  for (int j=0; j<3; j++)
    {
      locd[0][j] = 1.0/dt;
      locd[1][j] = 0.0;
    }
  
  locmat[6][0] = -( coord[1][1] - coord[1][2] ) * ( 2.0*locd[0][0] + locd[0][1] + locd[0][2] ) / 24.0 - ( coord[0][2] - coord[0][1] ) * ( 2.0*locd[1][0] + locd[1][1] + locd[1][2] ) / 24.0;
  locmat[6][1] = -( coord[1][1] - coord[1][2] ) * ( locd[0][0] + 2.0*locd[0][1] + locd[0][2] ) / 24.0 - ( coord[0][2] - coord[0][1] ) * ( locd[1][0] + 2.0*locd[1][1] + locd[1][2] ) / 24.0;
  locmat[6][2] = -( coord[1][1] - coord[1][2] ) * ( locd[0][0] + locd[0][1] + 2.0*locd[0][2] ) / 24.0 - ( coord[0][2] - coord[0][1] ) * ( locd[1][0] + locd[1][1] + 2.0*locd[1][2] ) / 24.0;

  locmat[7][0] = -( coord[1][2] - coord[1][0] ) * ( 2.0*locd[0][0] + locd[0][1] + locd[0][2] ) / 24.0 - ( coord[0][0] - coord[0][2] ) * ( 2.0*locd[1][0] + locd[1][1] + locd[1][2] ) / 24.0;
  locmat[7][1] = -( coord[1][2] - coord[1][0] ) * ( locd[0][0] + 2.0*locd[0][1] + locd[0][2] ) / 24.0 - ( coord[0][0] - coord[0][2] ) * ( locd[1][0] + 2.0*locd[1][1] + locd[1][2] ) / 24.0;
  locmat[7][2] = -( coord[1][2] - coord[1][0] ) * ( locd[0][0] + locd[0][1] + 2.0*locd[0][2] ) / 24.0 - ( coord[0][0] - coord[0][2] ) * ( locd[1][0] + locd[1][1] + 2.0*locd[1][2] ) / 24.0;

  locmat[8][0] = -( coord[1][0] - coord[1][1] ) * ( 2.0*locd[0][0] + locd[0][1] + locd[0][2] ) / 24.0 - ( coord[0][1] - coord[0][0] ) * ( 2.0*locd[1][0] + locd[1][1] + locd[1][2] ) / 24.0;
  locmat[8][1] = -( coord[1][0] - coord[1][1] ) * ( locd[0][0] + 2.0*locd[0][1] + locd[0][2] ) / 24.0 - ( coord[0][1] - coord[0][0] ) * ( locd[1][0] + 2.0*locd[1][1] + locd[1][2] ) / 24.0;
  locmat[8][2] = -( coord[1][0] - coord[1][1] ) * ( locd[0][0] + locd[0][1] + 2.0*locd[0][2] ) / 24.0 - ( coord[0][1] - coord[0][0] ) * ( locd[1][0] + locd[1][1] + 2.0*locd[1][2] ) / 24.0;

  for (int j=0; j<3; j++)
    {
      locd[0][j] = 0.0;
      locd[1][j] = 1.0/dt;
    }
  
  locmat[6][3] = -( coord[1][1] - coord[1][2] ) * ( 2.0*locd[0][0] + locd[0][1] + locd[0][2] ) / 24.0 - ( coord[0][2] - coord[0][1] ) * ( 2.0*locd[1][0] + locd[1][1] + locd[1][2] ) / 24.0;
  locmat[6][4] = -( coord[1][1] - coord[1][2] ) * ( locd[0][0] + 2.0*locd[0][1] + locd[0][2] ) / 24.0 - ( coord[0][2] - coord[0][1] ) * ( locd[1][0] + 2.0*locd[1][1] + locd[1][2] ) / 24.0;
  locmat[6][5] = -( coord[1][1] - coord[1][2] ) * ( locd[0][0] + locd[0][1] + 2.0*locd[0][2] ) / 24.0 - ( coord[0][2] - coord[0][1] ) * ( locd[1][0] + locd[1][1] + 2.0*locd[1][2] ) / 24.0;

  locmat[7][3] = -( coord[1][2] - coord[1][0] ) * ( 2.0*locd[0][0] + locd[0][1] + locd[0][2] ) / 24.0 - ( coord[0][0] - coord[0][2] ) * ( 2.0*locd[1][0] + locd[1][1] + locd[1][2] ) / 24.0;
  locmat[7][4] = -( coord[1][2] - coord[1][0] ) * ( locd[0][0] + 2.0*locd[0][1] + locd[0][2] ) / 24.0 - ( coord[0][0] - coord[0][2] ) * ( locd[1][0] + 2.0*locd[1][1] + locd[1][2] ) / 24.0;
  locmat[7][5] = -( coord[1][2] - coord[1][0] ) * ( locd[0][0] + locd[0][1] + 2.0*locd[0][2] ) / 24.0 - ( coord[0][0] - coord[0][2] ) * ( locd[1][0] + locd[1][1] + 2.0*locd[1][2] ) / 24.0;

  locmat[8][3] = -( coord[1][0] - coord[1][1] ) * ( 2.0*locd[0][0] + locd[0][1] + locd[0][2] ) / 24.0 - ( coord[0][1] - coord[0][0] ) * ( 2.0*locd[1][0] + locd[1][1] + locd[1][2] ) / 24.0;
  locmat[8][4] = -( coord[1][0] - coord[1][1] ) * ( locd[0][0] + 2.0*locd[0][1] + locd[0][2] ) / 24.0 - ( coord[0][1] - coord[0][0] ) * ( locd[1][0] + 2.0*locd[1][1] + locd[1][2] ) / 24.0;
  locmat[8][5] = -( coord[1][0] - coord[1][1] ) * ( locd[0][0] + locd[0][1] + 2.0*locd[0][2] ) / 24.0 - ( coord[0][1] - coord[0][0] ) * ( locd[1][0] + locd[1][1] + 2.0*locd[1][2] ) / 24.0;

  for (int j=0; j<3; j++)
    {
      locdat[j] = data->GetNodalEllipticCoeff(ind[j], 3);
    }

  locmat[6][6] += detJ * ( 3.0*locdat[0] + locdat[1] + locdat[2] ) / 60.0;
  locmat[6][7] += detJ * ( 2.0*locdat[0] + 2.0*locdat[1] + locdat[2] ) / 120.0;
  locmat[6][8] += detJ * ( 2.0*locdat[0] + locdat[1] + 2.0*locdat[2] ) / 120.0;
  locmat[7][6] = locmat[6][7];
  locmat[7][7] += detJ * ( locdat[0] + 3.0*locdat[1] + locdat[2] ) / 60.0;
  locmat[7][8] += detJ * ( locdat[0] + 2.0*locdat[1] + 2.0*locdat[2] ) / 120.0;
  locmat[8][6] = locmat[6][8];
  locmat[8][7] = locmat[7][8];
  locmat[8][8] += detJ * ( locdat[0] + locdat[1] + 3.0*locdat[2] ) / 60.0;
    
  /*
  cout << "locmat " << i << endl;
  for(int j=0; j<9; j++)
    {
      for(int k=0; k<9; k++)
	{
	  cout << locmat[j][k] << " ";
	}
      cout << endl;
    }
  cout<<endl<<endl;
*/

  
  for (int k=0; k<locd.Size(); k++) { delete []locd[k]; }
  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//=============================================================================

void LinearFEMGeoMech2D::ComputeReactionLocalSystem(int i, Array<int> &ind, 
						  Array<double *> &locmat)
{
  ;
}

//=============================================================================

void LinearFEMGeoMech2D::ComputeForceLocalSystem(int i, Array<int> &ind, 
					       Array<double> &locrhs)
{
  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];

  mesh->GetElementVerticesCoord(i, coord);
  double detJ = fabs((coord[0][1] - coord[0][0]) * (coord[1][2] 
		    - coord[1][0]) - (coord[1][1] - coord[1][0]) 
		    * (coord[0][2] - coord[0][0]));

  double f0,f1,f2,f3,f4,f5;
 
  f0 = data->GetNodalForceCoeff(ind[0], 0);
  f1 = data->GetNodalForceCoeff(ind[1], 0);
  f2 = data->GetNodalForceCoeff(ind[2], 0);
  f3 = data->GetNodalForceCoeff(ind[0], 1);
  f4 = data->GetNodalForceCoeff(ind[1], 1);
  f5 = data->GetNodalForceCoeff(ind[2], 1);

  locrhs[0] = (f2+f1+2.0*f0)/24.0;
  locrhs[1] = (f2+2.0*f1+f0)/24.0;
  locrhs[2] = (2.0*f2+f1+f0)/24.0;
  locrhs[3] = (f5+f4+2.0*f3)/24.0;
  locrhs[4] = (f5+2.0*f4+f3)/24.0;
  locrhs[5] = (2.0*f5+f4+f3)/24.0;

  for(int k=0; k<6; k++)
    locrhs[k] *= detJ;

  Array<double> locdat(3);

  for (int j=0; j<3; j++)
    locdat[j] = data->GetNodalForceCoeff(ind[j], 2);
  locrhs[6] = -locdat[0]/6.0 - locdat[1]/6.0 - locdat[2]/6.0;
  locrhs[7] = locdat[0]/6.0 + locdat[1]/6.0 + locdat[2]/6.0;

  for (int j=0; j<3; j++)
    locdat[j] = data->GetNodalForceCoeff(ind[j], 3);
  locrhs[6] += -locdat[0]/6.0 - locdat[1]/6.0 - locdat[2]/6.0;
  locrhs[8] = locdat[0]/6.0 + locdat[1]/6.0 + locdat[2]/6.0;

  //locrhs[6] = detJ * ( 2.0*locdat[0] + locdat[1] + locdat[2] ) / 24.0;
  //locrhs[7] = detJ * ( locdat[0] + 2.0*locdat[1] + locdat[2] ) / 24.0;
  //locrhs[8] = detJ * ( locdat[0] + locdat[1] + 2.0*locdat[2] ) / 24.0;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void LinearFEMGeoMech2D::ComputeNeumannLocalSystem(int i, Array<int> &ind, 
						 Array<double> &locrhs)
{
  ;
}

//==============================================================================

void LinearFEMGeoMech2D::ComputeNeumannLocalSystem(int b, int i, Array<int> &ind, 
					    Array<double> &locrhs)
{
  ;
}

//==============================================================================

void LinearFEMGeoMech2D::ComputeRobinLocalSystem(int i, int b, Array<int> &ind, 
					  Array<double *> &locmat, 
					  Array<double> &locrhs)
{
  ;
}

//==============================================================================

void LinearFEMGeoMech2D::ProjectAFunction(double (*func)(Array<double> &), 
					Vector &pval)
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

void LinearFEMGeoMech2D::SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
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
  //mat->Finalize();
  //mat->Print();
  //rhs.Print();
}

//==============================================================================

double LinearFEMGeoMech2D::ComputeL2Error(Function *exactsol, Vector &approxsol)
{ 
  double uh, tr = 0.0;
  double gp[8] = { 1.0/3.0, 1.0/3.0, 0.2, 0.6, 0.2, 0.2, 0.6, 0.2 };
  double gw[4] = { -27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0 };
  Array<double> x(3);
  Array<double> y(3);
   
  double err = 0.0;
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      mesh->GetElementVertices(i, ind);
       
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

      double gpcoord[2];
      for(int k=0; k<4; k++)
	{
	  double xx = gp[2*k];
	  double yy = gp[2*k+1];
	  gpcoord[0] = x[0] + t1*xx + t2*yy;
	  gpcoord[1] = y[0] + t3*xx + t4*yy;

	  uh = approxsol(ind[0]) * (1 - xx - yy) + approxsol(ind[1]) * xx + approxsol(ind[2]) * yy;
	  tr = exactsol->Eval(gpcoord);
	  err += gw[k] * (uh - tr)*(uh - tr) * detJ;	 
	}
    }  

  return sqrt(err);
  //cout << "The L2 norm error is: " << sqrt(err) << endl;
}
