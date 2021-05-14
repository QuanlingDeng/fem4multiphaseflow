#include "fem_header.h"

LinearFEMELAST2D::LinearFEMELAST2D(Mesh *_mesh, Data *_data) :
FEM(_mesh, _data)
{
    strcpy(Type, "LinearFEMELAST2D");
    NumGlobalDOF = 2*mesh->GetNV();
    NumLocalDOF = 6;
    rhs.SetSize(NumGlobalDOF);
    mat = new SparseMatrix(NumGlobalDOF, 25);
}

//============================================================================= 

void LinearFEMELAST2D::GetElementDOF(int i, Array<int> &ind)
{
  Array<int> tempind;
  GetElementVertices(i, tempind);
  ind.SetSize(6);  
  for(int i=0; i<3; i++)
    {
      ind[i] = tempind[i];
      ind[i+3] = tempind[i] + mesh->GetNV();
    }
}

//==============================================================================

void LinearFEMELAST2D::GetBdrElementDOF(int i, Array<int> &ind)
{
  Array<int> tempind;
  GetBdrElementVertices(i, tempind);
  ind.SetSize(4);
  for(int i=0; i<2; i++)
    {
      ind[i] = tempind[i];
      ind[i+2] = tempind[i] + mesh->GetNV();
    } 
}


//==============================================================================

void LinearFEMELAST2D::ComputeEllipticLocalSystem(int i, Array<int> &ind, 
						  Array<double *> &locmat)
{
 Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];

  mesh->GetElementVerticesCoord(i, coord);
  double det = (coord[0][1] - coord[0][0]) * (coord[1][2] 
		   - coord[1][0]) - (coord[1][1] - coord[1][0]) * 
		    (coord[0][2] - coord[0][0]);
  
  double intval = 1.0/(6.0*det);
 
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

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }

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

  /*
  cout << "locmat " << i << endl;
  for(int j=0; j<6; j++)
    {
      for(int k=0; k<6; k++)
	{
	  cout << locmat[j][k] << " ";
	}
      cout << endl;
    }
  */


}

//=============================================================================

void LinearFEMELAST2D::ComputeReactionLocalSystem(int i, Array<int> &ind, 
						  Array<double *> &locmat)
{
  ;
}

//=============================================================================

void LinearFEMELAST2D::ComputeForceLocalSystem(int i, Array<int> &ind, 
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
  
  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void LinearFEMELAST2D::ComputeNeumannLocalSystem(int i, Array<int> &ind, 
						 Array<double> &locrhs)
{

  double g00,g01,g10,g11;
  g00 = data->GetNodalNeumannCoeff(ind[0], 0);
  g01 = data->GetNodalNeumannCoeff(ind[1], 0);
  g10 = data->GetNodalNeumannCoeff(ind[0], 1);
  g11 = data->GetNodalNeumannCoeff(ind[1], 1);

  //cout << "g " << g00 << " " << g01 << " " << g10 << " " << g11 << endl;
  
  Array<double *> coord(2);
  coord[0] = new double[2];
  coord[1] = new double[2];
  
  mesh->GetBdrElementVerticesCoord(i, coord);

  double dx, dy, leng;
  dx = coord[0][0] - coord[0][1];
  dy = coord[1][0] - coord[1][1];

  leng = sqrt(pow(dx,2) + pow(dy,2));
  
  locrhs[0] = leng * ( 1.0/6.0 * g01 + 1.0/3.0 * g00 );
  locrhs[1] = leng * ( 1.0/3.0 * g01 + 1.0/6.0 * g00 );
  locrhs[2] = leng * ( 1.0/6.0 * g11 + 1.0/3.0 * g10 );
  locrhs[3] = leng * ( 1.0/3.0 * g11 + 1.0/6.0 * g10 );
  
  /*  cout << i << " ";
  for(int j=0; j<4; j++)
    cout << locrhs[j] << " ";
  cout << endl;
  */
  delete []coord[0];
  delete []coord[1];

}

//==============================================================================

void LinearFEMELAST2D::ComputeNeumannLocalSystem(int b, int i, Array<int> &ind, 
					    Array<double> &locrhs)
{
  int nb = mesh->GetNBdrs();

  double g00,g01,g10,g11;
  g00 = data->GetNeumannBdryVal(b, ind[0]);
  g01 = data->GetNeumannBdryVal(b, ind[1]);
  g10 = data->GetNeumannBdryVal(b+nb, ind[0]);
  g11 = data->GetNeumannBdryVal(b+nb, ind[1]);
  
  Array<double *> coord(2);
  coord[0] = new double[2];
  coord[1] = new double[2];
  
  mesh->GetBdrElementVerticesCoord(i, coord);

  double dx, dy, leng;
  dx = coord[0][0] - coord[0][1];
  dy = coord[1][0] - coord[1][1];

  leng = sqrt(pow(dx,2) + pow(dy,2));
  
  locrhs[0] = leng * ( 1.0/6.0 * g01 + 1.0/3.0 * g00 );
  locrhs[1] = leng * ( 1.0/3.0 * g01 + 1.0/6.0 * g00 );
  locrhs[2] = leng * ( 1.0/6.0 * g11 + 1.0/3.0 * g10 );
  locrhs[3] = leng * ( 1.0/3.0 * g11 + 1.0/6.0 * g10 );

 for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void LinearFEMELAST2D::ComputeRobinLocalSystem(int i, int b, Array<int> &ind, 
					  Array<double *> &locmat, 
					  Array<double> &locrhs)
{
  int nb = mesh->GetNBdrs();

  double g00,g01,g10,g11;
  g00 = -data->GetRobinCoeff(b, ind[0]);
  g01 = -data->GetRobinCoeff(b, ind[1]);

  g10 = -data->GetRobinCoeff(b+nb, ind[0]);
  g11 = -data->GetRobinCoeff(b+nb, ind[1]);

  double n00,n01,n10,n11;
  n00 = data->GetRobinBdryVal(b, ind[0]);
  n01 = data->GetRobinBdryVal(b, ind[1]);

  n10 = data->GetRobinBdryVal(b+nb, ind[0]);
  n11 = data->GetRobinBdryVal(b+nb, ind[1]);

  Array<double *> coord(2);
  coord[0] = new double[2];
  coord[1] = new double[2];
  
  mesh->GetBdrElementVerticesCoord(i, coord);

  double dx, dy, leng;
  dx = coord[0][0] - coord[0][1];
  dy = coord[1][0] - coord[1][1];
  leng = sqrt(pow(dx,2) + pow(dy,2));

  locmat[0][0] = leng * (g01 + 3.0*g00)/12.0;
  locmat[0][1] = leng * (g00 + g01)/12.0;

  locmat[0][2] = 0.0;
  locmat[0][3] = 0.0; 

  locmat[1][0] = leng * (g00 + g01)/12.0;
  locmat[1][1] = leng * (3.0*g01 + g00)/12.0;

  locmat[1][2] = 0.0;
  locmat[1][3] = 0.0;

  locmat[2][0] = 0.0;
  locmat[2][1] = 0.0;


  locmat[2][2] = leng * (g11 + 3.0*g10)/12.0;
  locmat[2][3] = leng * (g10 + g11)/12.0;

  locmat[3][0] = 0.0;
  locmat[3][1] = 0.0;

  locmat[3][2] = leng * (g10 + g11)/12.0;
  locmat[3][3] = leng * (3.0*g11 + g10)/12.0;

  
  locrhs[0] = leng * ( 1.0/6.0 * n01 + 1.0/3.0 * n00 );
  locrhs[1] = leng * ( 1.0/3.0 * n01 + 1.0/6.0 * n00 );

  locrhs[2] = leng * ( 1.0/6.0 * n11 + 1.0/3.0 * n10 );
  locrhs[3] = leng * ( 1.0/3.0 * n11 + 1.0/6.0 * n10 );

  delete []coord[0];
  delete []coord[1];
}

//==============================================================================

void LinearFEMELAST2D::ProjectAFunction(double (*func)(Array<double> &), 
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

void LinearFEMELAST2D::SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
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
  // mat->Print();
}

//==============================================================================

double LinearFEMELAST2D::ComputeL2Error(Function *exactsol, Vector &approxsol, double time)
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

      double gpcoord[3];
      gpcoord[2] = time;
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
