#include <iostream>
#include <fstream>
#include <cmath>
#include "../linalg/linalg_header.h"
#include "../mesh/mesh_header.h"
#include "../fem/fem_header.h"
#include "../timedependent/timedependent_header.h"

using namespace std;

//===========================================================

void ComputeConservativeFlux(DualMesh *dualmesh, Array<double *> &permdata, Array<double *> &divut, 
			     Array<double *> &forcedata, Array<double *> &QF, Vector &pressure, Array<double *> &flux)
{
  int n = dualmesh->GetNumLocalDOF();

  // Initialize flux to be zero
  for(int j=0; j<flux.Size(); j++)
    {
      for(int i=0; i<dualmesh->GetNE(); i++)
	{
	  flux[j][i] = 0.0;
	}
    }

  Vector cverr(dualmesh->GetNV()); // conservation error check
  cverr = 0.0;
  
  // Collect fluxes on each element
  for(int i=0; i<dualmesh->GetNE(); i++)
    {
      // Local flux and initialization
      Array<double> locflux(n);
      for(int k=0; k<locflux.Size(); k++) locflux[k] = 0.0;

      Array<int> ind;
      dualmesh->GetDOF(i, ind);
      
      // Get dualmesh info
      Array<double *> normals;
      Array<double> elengths;
      Array<double> areas;
      dualmesh->GetNormals(i, normals);
      dualmesh->GetEdgeLengths(i, elengths);
      dualmesh->GetAreas(i, areas);
      
      // Create Local Matrix
      SparseMatrix *locmat = new SparseMatrix(n,n);

      Array<double *> tem(n); // temperary array
      for(int k=0; k<n; k++) tem[k] = new double[n];
      for(int k=0; k<n; k++)
	for(int j=0; j<n; j++)
	  tem[k][j] = 0.0;
	
      Array<double> lock(n); // local elliptic coeffs
      Array<double> locmk(n); // local averaged elliptic coeff
      for(int j=0; j<n; j++) lock[j] = permdata[0][ind[j]];
      
      Array<double *> coord;
      dualmesh->GetElementVerticesCoord(i, coord);

      if(n==3) // local matrix for trianglular mesh
	{
	  locmk[0] = ( 5.0 * lock[0] + 5.0 * lock[1] + 2.0 * lock[2] ) / 12.0;
	  locmk[1] = ( 5.0 * lock[2] + 5.0 * lock[1] + 2.0 * lock[0] ) / 12.0;
	  locmk[2] = ( 5.0 * lock[0] + 5.0 * lock[2] + 2.0 * lock[1] ) / 12.0;
	  
	  double det = (coord[1][0] - coord[0][0]) * (coord[2][1] - coord[0][1])
	    - (coord[1][1] - coord[0][1]) * (coord[2][0] - coord[0][0]);
      
	  double dp0x, dp0y, dp1x, dp1y, dp2x, dp2y;
	  dp0x = (coord[2][1]-coord[1][1])/det;
	  dp0y = (coord[1][0]-coord[2][0])/det;
	  dp1x = (coord[2][1]-coord[0][1])/det;
	  dp1y = (coord[0][0]-coord[2][0])/det;
	  dp2x = (coord[0][1]-coord[1][1])/det;
	  dp2y = (coord[1][0]-coord[0][0])/det;

	  for(int k=0; k<3; k++)
	    {
	      tem[k][0] = locmk[0] * ( dp0x*normals[k][0] + dp0y*normals[k][1] ) * elengths[k];
	      tem[k][1] = -locmk[1] * ( dp1x*normals[k][0] + dp1y*normals[k][1] ) * elengths[k];
	      tem[k][2] = -locmk[2] * ( dp2x*normals[k][0] + dp2y*normals[k][1] ) * elengths[k];
	    }
	  
	  locmat->Elem(0,0) = tem[0][0] - tem[2][0];
	  locmat->Elem(0,1) = tem[0][1] - tem[2][1];
	  locmat->Elem(0,2) = tem[0][2] - tem[2][2];
	  
	  locmat->Elem(1,0) = tem[1][0] - tem[0][0];
	  locmat->Elem(1,1) = tem[1][1] - tem[0][1];
	  locmat->Elem(1,2) = tem[1][2] - tem[0][2];
	  
	  locmat->Elem(2,0) = tem[2][0] - tem[1][0];
	  locmat->Elem(2,1) = tem[2][1] - tem[1][1];
	  locmat->Elem(2,2) = tem[2][2] - tem[1][2];
	}
      else if(n==4) // local matrix for rectangular mesh
	{
	  locmk[0] = ( 3.0 * lock[0] + 3.0 * lock[1] + lock[2] + lock[3] ) / 8.0;
	  locmk[1] = ( 3.0 * lock[2] + 3.0 * lock[1] + lock[0] + lock[3] ) / 8.0;
	  locmk[2] = ( 3.0 * lock[2] + 3.0 * lock[3] + lock[0] + lock[1] ) / 8.0;
	  locmk[3] = ( 3.0 * lock[0] + 3.0 * lock[3] + lock[1] + lock[2] ) / 8.0;
      
	  double hx = coord[1][0] - coord[0][0]; 
	  double hy = coord[2][1] - coord[0][1];
	  double ryx = 0.5 * hy / hx;
	  double rxy = 0.5 * hx / hy;
	  tem[0][0] = locmk[0] * 0.75 * ryx;
	  tem[0][1] = -tem[0][0];
	  tem[0][2] = -locmk[0] * 0.25 * ryx;
	  tem[0][3] = -tem[0][2];
	  
	  tem[1][0] = locmk[1] * 0.25 * rxy;
	  tem[1][1] = locmk[1] * 0.75 * rxy;
	  tem[1][2] = -tem[1][1];
	  tem[1][3] = -tem[1][0];
	  
	  tem[2][0] = -locmk[2] * 0.25 * ryx;
	  tem[2][1] = -tem[2][0];
	  tem[2][2] = locmk[2] * 0.75 * ryx;
	  tem[2][3] = -tem[2][2];
	  
	  tem[3][0] = -locmk[3] * 0.75 * rxy;
	  tem[3][1] = -locmk[3] * 0.25 * rxy;
	  tem[3][2] = -tem[3][1];
	  tem[3][3] = -tem[3][0];
	  
	  locmat->Elem(0, 0) = tem[0][0] - tem[3][0];
	  locmat->Elem(0, 1) = tem[0][1] - tem[3][1];
	  locmat->Elem(0, 2) = tem[0][2] - tem[3][2];
	  locmat->Elem(0, 3) = tem[0][3] - tem[3][3];
	  
	  locmat->Elem(1, 0) = tem[1][0] - tem[0][0];
	  locmat->Elem(1, 1) = tem[1][1] - tem[0][1];
	  locmat->Elem(1, 2) = tem[1][2] - tem[0][2];
	  locmat->Elem(1, 3) = tem[1][3] - tem[0][3];
	  
	  locmat->Elem(2, 0) = tem[2][0] - tem[1][0];
	  locmat->Elem(2, 1) = tem[2][1] - tem[1][1];
	  locmat->Elem(2, 2) = tem[2][2] - tem[1][2];
	  locmat->Elem(2, 3) = tem[2][3] - tem[1][3];
	  
	  locmat->Elem(3, 0) = tem[3][0] - tem[2][0];
	  locmat->Elem(3, 1) = tem[3][1] - tem[2][1];
	  locmat->Elem(3, 2) = tem[3][2] - tem[2][2];
	  locmat->Elem(3, 3) = tem[3][3] - tem[2][3];
	}
      else { cout << "Mesh Error: Post-processing is valid for this mesh." << endl; }
      
      locmat->Finalize();
      
      // Create RHS
      Array<double> f(n);
      if(n==3)
	for(int k=0; k<n; k++) f[k] = ( forcedata[0][ind[k]] - divut[k][i] ) * areas[k];
      else if(n==4)
	for(int k=0; k<n; k++) f[k] = forcedata[0][ind[k]] * areas[k] - divut[k][i];
      else { cout << "Mesh Error: Post-processing is valid for this mesh." << endl; }
	
      Vector RHS(n);
      for(int k=0; k<n; k++) RHS(k) =  QF[k][i] + f[k];
	      
      // Solve the local system
      Vector s(n);
      s = 0.0;
      
      RectangularMatrix localconsmat(n);
      for(int j=0; j<n; j++)
	for(int k=0; k<n; k++)
	  localconsmat(j,k) = locmat->Elem(j,k);
      
      for(int k=1; k<n; k++) localconsmat(0, k) = 0.0;
      localconsmat(0,0) = 1.0;
      RHS(0) = pressure(ind[0]);
      
      DenseMatrixInverse invmat(localconsmat);
      invmat.Mult(RHS, s);

      for(int k=0; k<n; k++)
	{
	  //locflux[k] = s(k);
	  
	  for(int j=0; j<n; j++)
	    {
	      locflux[k] += tem[k][j] * s(j);
	    }
	}

      for(int k=0; k<n; k++) flux[k][i] = locflux[k];
      
      if(n==3)
	{
	  cverr(ind[0]) += -f[0] + locflux[0] - locflux[2];
	  cverr(ind[1]) += -f[1] + locflux[1] - locflux[0];
	  cverr(ind[2]) += -f[2] + locflux[2] - locflux[1];
	}
      else if(n==4)
	{
	  for(int k=0; k<n; k++) cverr(ind[k]) += -f[k] + locflux[k] - locflux[ (k+3)%4]; 
	}
      else { cout << "Mesh Error: Post-processing is valid for this mesh." << endl; }

      delete locmat;
      for(int k=0; k<n; k++) delete []tem[k];
    }

  cout << "Local Conservation Error for pressure" << endl;
  //for(int k=0; k<dualmesh->GetNV(); k++) if( fabs(cverr(k)) > 1.0e-7 ) cout << k << "\t" << cverr(k) << endl;
  //for(int k=0; k<dualmesh->GetNV(); k++) cout << k << "\t" << cverr(k) << endl;
}

//======================================================

void Computeuhdotn(DualMesh *dualmesh, Vector &ut, Array<double *> &uhdotn, Array<double> &buhdotn)
{
  int NE = dualmesh->GetNE();
  int NV = dualmesh->GetNV();
  int n = dualmesh->GetNumLocalDOF();
  Array<double *> normals(n);
  Array<double> lengths;
  Array<double *> ecoord(n+1);
  Array<double> areas;
  
  // Calculate uhdotn at internal boundaries
  for( int i=0; i<NE; i++)
    {
      dualmesh->GetCoordinates(i, ecoord);
      dualmesh->GetNormals(i, normals);
      dualmesh->GetEdgeLengths(i, lengths);
      dualmesh->GetAreas(i, areas);

      Array<int> ind;
      dualmesh->GetDOF(i, ind);
      
      if( n == 3)
	{
	  double utx = normals[0][0] * lengths[0] * ( 5.0 * ut(ind[0]) + 5.0 * ut(ind[1]) + 2.0 * ut(ind[2]) ) / 12.0;
	  double uty = normals[0][1] * lengths[0] * ( 5.0 * ut(ind[0]+NV) + 5.0 * ut(ind[1]+NV) + 2.0 * ut(ind[2]+NV) ) / 12.0; 
	  uhdotn[0][i] = utx + uty;

	  utx = normals[1][0] * lengths[1] * ( 5.0 * ut(ind[2]) + 5.0 * ut(ind[1]) + 2.0 * ut(ind[0]) ) / 12.0;
	  uty = normals[1][1] * lengths[1] * ( 5.0 * ut(ind[2]+NV) + 5.0 * ut(ind[1]+NV) + 2.0 * ut(ind[0]+NV) ) / 12.0; 
	  uhdotn[1][i] = utx + uty;

	  utx = normals[2][0] * lengths[2] * ( 5.0 * ut(ind[2]) + 5.0 * ut(ind[0]) + 2.0 * ut(ind[1]) ) / 12.0;
	  uty = normals[2][1] * lengths[2] * ( 5.0 * ut(ind[2]+NV) + 5.0 * ut(ind[0]+NV) + 2.0 * ut(ind[1]+NV) ) / 12.0; 
	  uhdotn[2][i] = utx + uty;
	}
      else if( n==4 )
	{
	  double utx = normals[0][0] * lengths[0] * ( 3.0 * ut(ind[0]) + 3.0 * ut(ind[1]) + ut(ind[2]) + ut(ind[3]) ) / 8.0;
	  double uty = normals[0][1] * lengths[0] * ( 3.0 * ut(ind[0]+NV) + 3.0 * ut(ind[1]+NV) + ut(ind[2]+NV) + ut(ind[3]+NV) ) / 8.0;
	  uhdotn[0][i] = utx + uty;

	  utx = normals[1][0] * lengths[1] * ( 3.0 * ut(ind[2]) + 3.0 * ut(ind[1]) + ut(ind[0]) + ut(ind[3]) ) / 8.0;
	  uty = normals[1][1] * lengths[1] * ( 3.0 * ut(ind[2]+NV) + 3.0 * ut(ind[1]+NV) + ut(ind[0]+NV) + ut(ind[3]+NV) ) / 8.0;
	  uhdotn[1][i] = utx + uty;

	  utx = normals[2][0] * lengths[2] * ( 3.0 * ut(ind[2]) + 3.0 * ut(ind[3]) + ut(ind[0]) + ut(ind[1]) ) / 8.0;
	  uty = normals[2][1] * lengths[2] * ( 3.0 * ut(ind[2]+NV) + 3.0 * ut(ind[3]+NV) + ut(ind[0]+NV) + ut(ind[1]+NV) ) / 8.0;
	  uhdotn[2][i] = utx + uty;

	  utx = normals[3][0] * lengths[3] * ( 3.0 * ut(ind[0]) + 3.0 * ut(ind[3]) + ut(ind[1]) + ut(ind[2]) ) / 8.0;
	  uty = normals[3][1] * lengths[3] * ( 3.0 * ut(ind[0]+NV) + 3.0 * ut(ind[3]+NV) + ut(ind[1]+NV) + ut(ind[2]+NV) ) / 8.0;
	  uhdotn[3][i] = utx + uty;
	}
      else { cout << "Mesh Error: Post-processing is valid for this mesh." << endl; }
    }

  // Calculate buhdotn at the global boundaries
  for(int i=0; i<NV; i++) buhdotn[i] = 0.0;

  Mesh *mesh = dualmesh->GetMesh();
  Array<double> tem(2*dualmesh->GetNBE());
  for( int i=0; i<mesh->GetNBE(); i++)
    {
      int nb = mesh->GetNBdrs();      

      Array<int> ind;
      dualmesh->GetBdrElementDOF(i, ind);

      Array<double *> coord(2);
      coord[0] = new double[2];
      coord[1] = new double[2];

      dualmesh->GetBdrElementVerticesCoord(i, coord);
      double len = sqrt( (coord[0][0] - coord[0][1])*(coord[0][0] - coord[0][1]) + (coord[1][0] - coord[1][1])*(coord[1][0] - coord[1][1]) ) / 2.0;
      
      if(i/nb==0)
	{
	  buhdotn[ind[0]] += -len * ( 3.0 * ut(ind[0]+NV) + ut(ind[1]+NV) ) / 4.0; 
	  buhdotn[ind[1]] += -len * ( 3.0 * ut(ind[1]+NV) + ut(ind[0]+NV) ) / 4.0; 
	}
      else if(i/nb==1)
	{
	  buhdotn[ind[0]] += len * ( 3.0 * ut(ind[0]) + ut(ind[1]) ) / 4.0; 
	  buhdotn[ind[1]] += len * ( 3.0 * ut(ind[1]) + ut(ind[0]) ) / 4.0; 
	}
      else if(i/nb==2)
	{
	  buhdotn[ind[0]] += len * ( 3.0 * ut(ind[0]+NV) + ut(ind[1]+NV) ) / 4.0; 
	  buhdotn[ind[1]] += len * ( 3.0 * ut(ind[1]+NV) + ut(ind[0]+NV) ) / 4.0; 
	}
      else
	{
	  buhdotn[ind[0]] += -len * ( 3.0 * ut(ind[0]) + ut(ind[1]) ) / 4.0; 
	  buhdotn[ind[1]] += -len * ( 3.0 * ut(ind[1]) + ut(ind[0]) ) / 4.0; 
	}
      
      delete []coord[0];
      delete []coord[1];
    }
}

//======================================================

void calculatedivu(Mesh *mesh, const Vector &u, Array<double *> &divu)
{
  int nv = mesh->GetNV();
  int ne = mesh->GetNE();

  // Calculate elemental divu
  for(int k=0; k<ne; k++)
    {
      Array<int> ind;
      mesh->GetElementVertices(k, ind);

      Array<double *> vcoords(mesh->GetDim());
      for (int i=0; i<vcoords.Size(); i++)
	{
	  vcoords[i] = new double[ind.Size()];
	}
      mesh->GetElementVerticesCoord(k, vcoords);

      Array<double> alpha(ind.Size());
      Array<double> beta(ind.Size());
      for(int i=0; i<ind.Size(); i++)
	{
	  alpha[i] = u(ind[i]);
	  beta[i] = u(ind[i] + nv);
	}

      if(mesh->GetMType()==Element::TRIANGLE)
	{
	  double x0 = vcoords[0][0];
	  double y0 = vcoords[1][0];
	  double x1 = vcoords[0][1];
	  double y1 = vcoords[1][1];
	  double x2 = vcoords[0][2];
	  double y2 = vcoords[1][2];
	  double det = -x1*y0 + x2*y0 + x0*y1 - x2*y1 - x0*y2 + x1*y2; 
	  
	  double du1x = (alpha[0]*(y1-y2) + alpha[1]*(y2-y0) + alpha[2]*(y0-y1))/det;
	  double du2y = (beta[0]*(x2-x1) + beta[1]*(x0-x2) + beta[2]*(x1-x0))/det;
	  double du = du1x + du2y;

	  for(int i=0; i<ind.Size(); i++) divu[i][k] = du;
	}

      else if(mesh->GetMType()==Element::QUADRILATERAL)
	{	
	  //Assuming rectangles	  
	  double dx = fabs(vcoords[0][1]-vcoords[0][0]);
	  double dy = fabs(vcoords[1][3]-vcoords[1][0]);

	  divu[0][k] = (alpha[1]-alpha[0])/dx;
	  divu[1][k] = (alpha[0]-alpha[1])/dx;
	  divu[2][k] = (alpha[3]-alpha[2])/dx;
	  divu[3][k] = (alpha[2]-alpha[3])/dx;

	  divu[0][k] = (beta[3]-beta[0])/dy;
	  divu[1][k] = (beta[2]-beta[1])/dy;
	  divu[2][k] = (beta[2]-beta[3])/dy;
	  divu[3][k] = (beta[0]-beta[0])/dy;
	}
	 
      for(int i=0; i<vcoords.Size(); i++) delete []vcoords[i];
    }
}

//======================================================

double ComputeDiscL2Error(Mesh *mesh, Function *exactsol, Array<double *> &sol)
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

      if(ind.Size()==3)
	{
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
	      
	      uh = sol[0][i] * (1 - xx - yy) + sol[1][i] * xx + sol[2][i] * yy;
	      tr = exactsol->Eval(gpcoord);
	      err += gw[k] * (uh - tr)*(uh - tr) * detJ;	 
	    }
	}
      else if(ind.Size()==4)
	{
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
	      
	      uh = sol[0][i] * (1 - xx - yy) + sol[1][i] * xx + sol[2][i] * yy;
	      tr = exactsol->Eval(gpcoord);
	      err += gw[k] * (uh - tr)*(uh - tr) * detJ;	 
	    }

	  for(int k=0; k<3; k++)
	    {
	      x[k] = mesh->GetVertex(ind[(k+2)%4], 0);
	      y[k] = mesh->GetVertex(ind[(k+2)%4], 1);	  
	    }
	  
	  t1 = x[1] - x[0];
	  t2 = x[2] - x[0];
	  t3 = y[1] - y[0];
	  t4 = y[2] - y[0];
	  detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<4; k++)
	    {
	      double xx = gp[2*k];
	      double yy = gp[2*k+1];
	      gpcoord[0] = x[0] + t1*xx + t2*yy;
	      gpcoord[1] = y[0] + t3*xx + t4*yy;
	      
	      uh = sol[2][i] * (1 - xx - yy) + sol[3][i] * xx + sol[0][i] * yy;
	      tr = exactsol->Eval(gpcoord);
	      err += gw[k] * (uh - tr)*(uh - tr) * detJ;	 
	    }
	}
      else { cout << "Mesh Error: Post-processing is invalid for this mesh." << endl; }
    }  

  return sqrt(err);
  //cout<<"The L2 norm error is: "<< sqrt(err)<<endl;
}

//======================================================

double ComputeDiscH1SemiNormError(Mesh *mesh, Function *dirx, Function *diry, Array<double *> &sol)
{ 
  double gp[8] = { 1.0/3.0, 1.0/3.0, 0.2, 0.6, 0.2, 0.2, 0.6, 0.2 };
  double gw[4] = { -27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0 };
  Array<double> x(3);
  Array<double> y(3);
  double coord[2];
   
  double err = 0.0;
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      mesh->GetElementVertices(i, ind);

      if(ind.Size()==3)
	{
	  for(int k=0; k<3; k++)
	    {
	      x[k] = mesh->GetVertex(ind[k], 0);
	      y[k] = mesh->GetVertex(ind[k], 1);	  
	    }
	  
	  double t1 = x[1] - x[0];
	  double t2 = x[2] - x[0];
	  double t3 = y[1] - y[0];
	  double t4 = y[2] - y[0];
	  
	  double s1 = y[1] - y[2];
	  double s2 = y[2] - y[0];
	  double s3 = y[0] - y[1];
	  double s4 = x[2] - x[1];
	  double s5 = x[0] - x[2];
	  double s6 = x[1] - x[0];
	  
	  double detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<4; k++)
	    {
	      double xx = gp[2*k];
	      double yy = gp[2*k+1];
	      coord[0] = x[0] + t1*xx + t2*yy;
	      coord[1] = y[0] + t3*xx + t4*yy;
	      
	      double uhx = ( sol[0][i]*s1 + sol[1][i]*s2 + sol[2][i]*s3 ) / detJ;
	      uhx -= dirx->Eval(coord);
	      uhx = uhx * uhx;
	      double uhy = ( sol[0][i]*s4 + sol[1][i]*s5 + sol[2][i]*s6 ) / detJ;
	      uhy -= diry->Eval(coord);
	      uhy = uhy * uhy;
	      err += ( uhx + uhy ) * gw[k] * detJ;
	    }
	}
      else if(ind.Size()==4)
	{
	  for(int k=0; k<3; k++)
	    {
	      x[k] = mesh->GetVertex(ind[k], 0);
	      y[k] = mesh->GetVertex(ind[k], 1);	  
	    }
	  
	  double t1 = x[1] - x[0];
	  double t2 = x[2] - x[0];
	  double t3 = y[1] - y[0];
	  double t4 = y[2] - y[0];
	  
	  double s1 = y[1] - y[2];
	  double s2 = y[2] - y[0];
	  double s3 = y[0] - y[1];
	  double s4 = x[2] - x[1];
	  double s5 = x[0] - x[2];
	  double s6 = x[1] - x[0];
	  
	  double detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<4; k++)
	    {
	      double xx = gp[2*k];
	      double yy = gp[2*k+1];
	      coord[0] = x[0] + t1*xx + t2*yy;
	      coord[1] = y[0] + t3*xx + t4*yy;
	      
	      double uhx = ( sol[0][i]*s1 + sol[1][i]*s2 + sol[2][i]*s3 ) / detJ;
	      uhx -= dirx->Eval(coord);
	      uhx = uhx * uhx;
	      double uhy = ( sol[0][i]*s4 + sol[1][i]*s5 + sol[2][i]*s6 ) / detJ;
	      uhy -= diry->Eval(coord);
	      uhy = uhy * uhy;
	      err += ( uhx + uhy ) * gw[k] * detJ;
	    }
	  
	  for(int k=0; k<3; k++)
	    {
	      x[k] = mesh->GetVertex(ind[(k+2)%4], 0);
	      y[k] = mesh->GetVertex(ind[(k+2)%4], 1);	  
	    }
	  
	  t1 = x[1] - x[0];
	  t2 = x[2] - x[0];
	  t3 = y[1] - y[0];
	  t4 = y[2] - y[0];
	  
	  s1 = y[1] - y[2];
	  s2 = y[2] - y[0];
	  s3 = y[0] - y[1];
	  s4 = x[2] - x[1];
	  s5 = x[0] - x[2];
	  s6 = x[1] - x[0];
	  
	  detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<4; k++)
	    {
	      double xx = gp[2*k];
	      double yy = gp[2*k+1];
	      coord[0] = x[0] + t1*xx + t2*yy;
	      coord[1] = y[0] + t3*xx + t4*yy;
	      
	      double uhx = ( sol[0][i]*s1 + sol[1][i]*s2 + sol[2][i]*s3 ) / detJ;
	      uhx -= dirx->Eval(coord);
	      uhx = uhx * uhx;
	      double uhy = ( sol[0][i]*s4 + sol[1][i]*s5 + sol[2][i]*s6 ) / detJ;
	      uhy -= diry->Eval(coord);
	      uhy = uhy * uhy;
	      err += ( uhx + uhy ) * gw[k] * detJ;
	    }
	}
      else { cout << "Mesh Error: Post-processing is invalid for this mesh." << endl; }
    }

  return sqrt(err);
}

//======================================================

double ComputeH1SemiNormError(Mesh *mesh, Function *dirx, Function *diry, Vector &sol)
{ 
  double gp[8] = { 1.0/3.0, 1.0/3.0, 0.2, 0.6, 0.2, 0.2, 0.6, 0.2 };
  double gw[4] = { -27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0 };
  Array<double> x(3);
  Array<double> y(3);
  double coord[2];
   
  double err = 0.0;
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      mesh->GetElementVertices(i, ind);

      if(ind.Size()==3)
	{
	  for(int k=0; k<3; k++)
	    {
	      x[k] = mesh->GetVertex(ind[k], 0);
	      y[k] = mesh->GetVertex(ind[k], 1);	  
	    }
	  
	  double t1 = x[1] - x[0];
	  double t2 = x[2] - x[0];
	  double t3 = y[1] - y[0];
	  double t4 = y[2] - y[0];
	  
	  double s1 = y[1] - y[2];
	  double s2 = y[2] - y[0];
	  double s3 = y[0] - y[1];
	  double s4 = x[2] - x[1];
	  double s5 = x[0] - x[2];
	  double s6 = x[1] - x[0];
	  
	  double detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<4; k++)
	    {
	      double xx = gp[2*k];
	      double yy = gp[2*k+1];
	      coord[0] = x[0] + t1*xx + t2*yy;
	      coord[1] = y[0] + t3*xx + t4*yy;
	      
	      double uhx = ( sol(ind[0])*s1 + sol(ind[1])*s2 + sol(ind[2])*s3 ) / detJ;
	      uhx -= dirx->Eval(coord);
	      uhx = uhx * uhx;
	      double uhy = ( sol(ind[0])*s4 + sol(ind[1])*s5 + sol(ind[2])*s6 ) / detJ;
	      uhy -= diry->Eval(coord);
	      uhy = uhy * uhy;
	      err += ( uhx + uhy ) * gw[k] * detJ;
	    }
	}
      else if(ind.Size()==4)
	{
	  for(int k=0; k<3; k++)
	    {
	      x[k] = mesh->GetVertex(ind[k], 0);
	      y[k] = mesh->GetVertex(ind[k], 1);	  
	    }
	  
	  double t1 = x[1] - x[0];
	  double t2 = x[2] - x[0];
	  double t3 = y[1] - y[0];
	  double t4 = y[2] - y[0];
	  
	  double s1 = y[1] - y[2];
	  double s2 = y[2] - y[0];
	  double s3 = y[0] - y[1];
	  double s4 = x[2] - x[1];
	  double s5 = x[0] - x[2];
	  double s6 = x[1] - x[0];
	  
	  double detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<4; k++)
	    {
	      double xx = gp[2*k];
	      double yy = gp[2*k+1];
	      coord[0] = x[0] + t1*xx + t2*yy;
	      coord[1] = y[0] + t3*xx + t4*yy;
	      
	      double uhx = ( sol(ind[0])*s1 + sol(ind[1])*s2 + sol(ind[2])*s3 ) / detJ;
	      uhx -= dirx->Eval(coord);
	      uhx = uhx * uhx;
	      double uhy = ( sol(ind[0])*s4 + sol(ind[1])*s5 + sol(ind[2])*s6 ) / detJ;
	      uhy -= diry->Eval(coord);
	      uhy = uhy * uhy;
	      err += ( uhx + uhy ) * gw[k] * detJ;
	    }
	  
	  for(int k=0; k<3; k++)
	    {
	      x[k] = mesh->GetVertex(ind[(k+2)%4], 0);
	      y[k] = mesh->GetVertex(ind[(k+2)%4], 1);	  
	    }
	  
	  t1 = x[1] - x[0];
	  t2 = x[2] - x[0];
	  t3 = y[1] - y[0];
	  t4 = y[2] - y[0];
	  
	  s1 = y[1] - y[2];
	  s2 = y[2] - y[0];
	  s3 = y[0] - y[1];
	  s4 = x[2] - x[1];
	  s5 = x[0] - x[2];
	  s6 = x[1] - x[0];
	  
	  detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<4; k++)
	    {
	      double xx = gp[2*k];
	      double yy = gp[2*k+1];
	      coord[0] = x[0] + t1*xx + t2*yy;
	      coord[1] = y[0] + t3*xx + t4*yy;
	      
	      double uhx = ( sol(ind[0])*s1 + sol(ind[1])*s2 + sol(ind[2])*s3 ) / detJ;
	      uhx -= dirx->Eval(coord);
	      uhx = uhx * uhx;
	      double uhy = ( sol(ind[0])*s4 + sol(ind[1])*s5 + sol(ind[2])*s6 ) / detJ;
	      uhy -= diry->Eval(coord);
	      uhy = uhy * uhy;
	      err += ( uhx + uhy ) * gw[k] * detJ;
	    }
	}
      else { cout << "Mesh Error: Post-processing is invalid for this mesh." << endl; }
    }

  return sqrt(err);
}

//======================================================

void calculategradp(Mesh *mesh, const Vector &p, Array<double *> &gradp)
{
  int nv = mesh->GetNV();
  int ne = mesh->GetNE();

  gradp.SetSize(mesh->GetDim());
  for(int i=0; i<mesh->GetDim(); i++)
    {
      gradp[i] = new double[nv];
      for(int k=0; k<nv; k++)
	{
	  gradp[i][k] = 0.0;
	}
    }

  Vector count(nv);
  count = 0.0;

  for(int k=0; k<ne; k++)
    {
      Array<int> ind;
      mesh->GetElementVertices(k, ind);
      Array<double *> vcoords(mesh->GetDim());
      for (int i=0; i<vcoords.Size(); i++)
	{
	  vcoords[i] = new double[ind.Size()];
	}
      mesh->GetElementVerticesCoord(k, vcoords);

      Array<double> alpha(ind.Size());
      for(int i=0; i<ind.Size(); i++)
	{
	  alpha[i] = p(ind[i]);
	}

      if(mesh->GetMType()==Element::TRIANGLE)
	{
	  double x0 = vcoords[0][0];
	  double y0 = vcoords[1][0];
	  double x1 = vcoords[0][1];
	  double y1 = vcoords[1][1];
	  double x2 = vcoords[0][2];
	  double y2 = vcoords[1][2];
	  double det = -x1*y0 + x2*y0 + x0*y1 - x2*y1 - x0*y2 + x1*y2; 
	  
	  double dpx = (alpha[0]*(y1-y2) + alpha[1]*(y2-y0) + alpha[2]*(y0-y1))/det;
	  double dpy = (alpha[0]*(x2-x1) + alpha[1]*(x0-x2) + alpha[2]*(x1-x0))/det;

	  for(int i=0; i<ind.Size(); i++)
	    {
	      gradp[0][ind[i]] += dpx;
	      gradp[1][ind[i]] += dpy;
	      count(ind[i]) += 1.0;
	    }
	}

      if(mesh->GetMType()==Element::QUADRILATERAL)
	{
	  double dx = fabs(vcoords[0][1]-vcoords[0][0]);
	  double dy = fabs(vcoords[1][3]-vcoords[1][0]);

	  gradp[0][ind[0]] += (alpha[1]-alpha[0])/dx;
	  gradp[0][ind[1]] += (alpha[0]-alpha[1])/dx;
	  gradp[0][ind[2]] += (alpha[3]-alpha[2])/dx;
	  gradp[0][ind[3]] += (alpha[2]-alpha[3])/dx;

	  gradp[1][ind[0]] += (alpha[3]-alpha[0])/dy;
	  gradp[1][ind[1]] += (alpha[2]-alpha[1])/dy;
	  gradp[1][ind[2]] += (alpha[2]-alpha[3])/dy;
	  gradp[1][ind[3]] += (alpha[0]-alpha[0])/dy;

	  for(int i=0; i<ind.Size(); i++)
	    {
	      count(ind[i]) += 1.0;
	    }
	}
  
      for(int i=0; i<vcoords.Size(); i++)
	delete []vcoords[i];
    }

  for(int k=0; k<nv; k++)
    {
      gradp[0][k] /= ((double) count(k)); 
      gradp[1][k] /= ((double) count(k)); 
    }

}

//======================================================

int main(int argc, char *argv[])
{

 //-----------------------------
  // Input File Information
  //-----------------------------
  char fname[100]; ifstream input_file;

  if( argc > 1 )
    input_file.open( argv[1] );
  else
    {
      cout << " Enter input file name:... ";
      cin >> fname;
      input_file.open( fname );
    }
  
  if( !input_file )
    {
      cerr << " Cannot open the specified file. \n \n" << endl;
      return 3;
    }

  const int bufflen = 256; char buffer[bufflen];

  //---------------------------------------------------------
  //  Get Time March Data 
  //---------------------------------------------------------

  double finaltime;
  int numfinesteps;
  
  input_file.getline( buffer, bufflen, ' ' );
  input_file >> finaltime >> numfinesteps;
  double dt = finaltime/(numfinesteps);

  cout << endl << "Time Discretiazations: " << endl;
  cout << "Final Time: " << finaltime << endl;
  cout << "Number of Fine Time Steps: " << numfinesteps << endl;
  cout << "Fine Time Step Length: " << dt << endl; 
  cout << endl;

  //-----------------------------
  // Build the Global Mesh
  //-----------------------------

  //Gather Mesh Data
  Array<double> L(4);
  input_file.getline( buffer, bufflen, ' ' );
  input_file >> L[0] >> L[1] >> L[2] >> L[3];

  Array<int> N(2);
  int meshtype;
  input_file.getline( buffer, bufflen, ' ' );
  input_file >> N[0] >> N[1] >> meshtype;

  cout << "Mesh Configuration: " << endl;
  cout << "Nx: " << N[0] << " " << "Ny: " << N[1] << endl;
  cout << "Lx: " << L[1]-L[0] << " " << "Ly: " << L[3]-L[2] << endl;
  if(meshtype==Element::TRIANGLE){cout << "meshtype: Triangle" << endl;}
  if(meshtype==Element::QUADRILATERAL){cout << "meshtype: Quadrilateral" << endl; }
  cout << endl;

  //build mesh
  Mesh *mesh;
  if(meshtype==Element::TRIANGLE){mesh = new TMesh<2>(N, L, Element::TRIANGLE);}
  if(meshtype==Element::QUADRILATERAL){ mesh = new TMesh<2>(N, L, Element::QUADRILATERAL); }
  cout << "Mesh Built" << endl;

  //print mesh
  mesh->ParaviewPrint(meshtype);

  int nv = mesh->GetNV(); //number of vertices
  double h = (L[1]-L[0])/N[0];
  double hx = (L[1]-L[0])/N[0];
  double hy = (L[3]-L[2])/N[1];
  //---------------------------------------------------------
  // Gather Data for Pressure Equation
  //--------------------------------------------------------

  Array<double *> permdata(1);
  permdata[0] = new double[nv];

  char permflag[bufflen]; 
  //input_file.getline( buffer, bufflen, ' ' );
  ReadIdentifier(input_file, buffer, bufflen);
  input_file >> permflag;

  if (!strcmp(permflag, "equation"))
    {
      Function *permfunc;
      ReadIdentifier(input_file, buffer, bufflen);
      permfunc = ReadFunction(input_file);

      for(int k=0; k<nv; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  permdata[0][k] = permfunc->Eval(coord);
	}      
      delete permfunc;
    }

  //Gather Force Data
  ReadIdentifier(input_file, buffer, bufflen);
  Function *force = ReadFunction(input_file);

  //Gather Boundary Data
  int nbdry = mesh->GetNBdrs();
  Array<double *> nbdryval(nbdry);
  Array<double *> rbdryval(nbdry);
  Array<double *> rcoeff(nbdry);
  Array<double *> dbdryval(nbdry);
  Array<bool> BdrNeumann(nbdry);
  Array<bool> BdrRobin(nbdry);
  Array<bool> BdrDirichlet(nbdry);

  char bdryflag[bufflen]; 
  for(int b=0; b<nbdry; b++)
    {
      int bnv = mesh->GetNBdryElem(b)+1;
      nbdryval[b] = new double[bnv];
      rbdryval[b] = new double[bnv];
      rcoeff[b] = new double[bnv];
      dbdryval[b] = new double[bnv];
      
      input_file.getline( buffer, bufflen, ' ' );
      input_file >> bdryflag;

      if (!strcmp(bdryflag, "neumann"))
	{
	  BdrNeumann[b] = true;
	  BdrRobin[b] = false;
	  BdrDirichlet[b] = false;

	  ReadIdentifier(input_file, buffer, bufflen);
	  Function *bdryval = ReadFunction(input_file);
	  for(int k=0; k<bnv; k++)
	    {
	      int indie = mesh->GetBdryGindex(b, k);
	      double* coord = mesh->GetVertex(indie);

	      nbdryval[b][k] = bdryval->Eval(coord);
	     
	      rbdryval[b][k] = 0.0;
	      rcoeff[b][k] = 0.0;
	      dbdryval[b][k] = 0.0;
	    }
	  delete bdryval;
	}

      if (!strcmp(bdryflag, "dirichlet"))
	{
	  BdrDirichlet[b] = true;
	  BdrNeumann[b] = false;
	  BdrRobin[b] = false;

	  ReadIdentifier(input_file, buffer, bufflen);
	  Function *bdryval = ReadFunction(input_file);
	  for(int k=0; k<bnv; k++)
	    {
	      int indie = mesh->GetBdryGindex(b, k);
	      double* coord = mesh->GetVertex(indie);

	      dbdryval[b][k] = bdryval->Eval(coord);
	     
	      rbdryval[b][k] = 0.0;
	      rcoeff[b][k] = 0.0;
	      nbdryval[b][k] = 0.0;
	    }
	  delete bdryval;
	}
      
      if (!strcmp(bdryflag, "robin"))
	{
	  BdrRobin[b] = true;
	  BdrNeumann[b] = false;
	  BdrDirichlet[b] = false;

	  ReadIdentifier(input_file, buffer, bufflen);
	  Function *bdryval = ReadFunction(input_file);
	  ReadIdentifier(input_file, buffer, bufflen);
	  Function *coeff = ReadFunction(input_file);
	  for(int k=0; k<bnv; k++)
	    {
	      int indie = mesh->GetBdryGindex(b, k);
	      double* coord = mesh->GetVertex(indie);

	      rbdryval[b][k] = bdryval->Eval(coord);
	      rcoeff[b][k] = coeff->Eval(coord);

	      dbdryval[b][k] = 0.0;
	      nbdryval[b][k] = 0.0;
	    }
	  delete bdryval;
	  delete coeff;
	}
      
    }

  //set data for pressure equation
  Data *pressuredata;
  pressuredata = new Data();
  pressuredata->SetEllipticData(permdata);
  pressuredata->SetNeumannData(nbdryval, BdrNeumann);
  pressuredata->SetRobinData(rcoeff, rbdryval, BdrRobin);
  pressuredata->SetDirichletData(dbdryval, BdrDirichlet);
  cout << "Pressure Data Set" << endl;


  //---------------------------------------------------------
  //  Gather Data for elasticity part
  //--------------------------------------------------------- 

  //Gather elasticitiy Data
  Array<double *> lamecoeff(2);
  lamecoeff[0] = new double[nv];
  lamecoeff[1] = new double[nv];
  Vector elmod(nv);


  char elmodflag[bufflen]; 
  ReadIdentifier(input_file, buffer, bufflen);
  input_file >> elmodflag;
  
  if (!strcmp(elmodflag, "equation"))
    {
      Function *modfunc;
      ReadIdentifier(input_file, buffer, bufflen);
      modfunc = ReadFunction(input_file);
      
      double prat = 0.2;
      for(int k=0; k<nv; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  double data = modfunc->Eval(coord);
	  elmod(k) = data;
	  lamecoeff[0][k] = data*prat/((1.0+prat)*(1.0-2.0*prat)); // -0.5;
	  lamecoeff[1][k] = data/(1+prat); // 1.0;
	}
      delete modfunc;
    }
  mesh->ParaviewPrint(elmod, meshtype);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *forcex = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *forcey = ReadFunction(input_file);

  //Set Boundary Data
  int nbredges = mesh->GetNBdrs();
  int nbdryel = 2*nbredges;
  Array<double *> elnbdryval(nbdryel);
  Array<double *> elrbdryval(nbdryel);
  Array<double *> elrcoeff(nbdryel);
  Array<double *> eldbdryval(nbdryel);
  Array<bool> elBdrNeumann(nbdryel);
  Array<bool> elBdrRobin(nbdryel);
  Array<bool> elBdrDirichlet(nbdryel);

  for(int b=0; b<nbredges; b++)
    {
      int bnv = mesh->GetNBdryElem(b)+1;
      for(int i=0; i<2; i++)
	{
	  int bb = b + i*nbredges;
	  elnbdryval[bb] = new double[bnv];
	  elrbdryval[bb] = new double[bnv];
	  elrcoeff[bb] = new double[bnv];
	  eldbdryval[bb] = new double[bnv];
	  
	  input_file.getline( buffer, bufflen, ' ' );
	  input_file >> bdryflag;
	  
	  if (!strcmp(bdryflag, "neumann"))
	    {
	      elBdrNeumann[bb] = true;
	      elBdrRobin[bb] = false;
	      elBdrDirichlet[bb] = false;
	      
	      ReadIdentifier(input_file, buffer, bufflen);
	      Function *bdryval = ReadFunction(input_file);
	      for(int k=0; k<bnv; k++)
		{
		  int indie = mesh->GetBdryGindex(b, k);
		  double* coord = mesh->GetVertex(indie);
	
		  elnbdryval[bb][k] = bdryval->Eval(coord);
		  
		  elrbdryval[bb][k] = 0.0;
		  elrcoeff[bb][k] = 0.0;
		  eldbdryval[bb][k] = 0.0;
		}

	      delete bdryval;
	    }
	  
	  if (!strcmp(bdryflag, "dirichlet"))
	    {
	      elBdrDirichlet[bb] = true;
	      elBdrNeumann[bb] = false;
	      elBdrRobin[bb] = false;

	      ReadIdentifier(input_file, buffer, bufflen);
	      Function *bdryval = ReadFunction(input_file);
	      for(int k=0; k<bnv; k++)
		{
		  int indie = mesh->GetBdryGindex(b, k);
		  double* coord = mesh->GetVertex(indie);

		  eldbdryval[bb][k] = bdryval->Eval(coord);
	     
		  elrbdryval[bb][k] = 0.0;
		  elrcoeff[bb][k] = 0.0;
		  elnbdryval[bb][k] = 0.0;
		}
	      delete bdryval;
	    }
      
	  if (!strcmp(bdryflag, "robin"))
	    {
	      elBdrRobin[bb] = true;
	      elBdrNeumann[bb] = false;
	      elBdrDirichlet[bb] = false;

	      ReadIdentifier(input_file, buffer, bufflen);
	      Function *bdryval = ReadFunction(input_file);
	      ReadIdentifier(input_file, buffer, bufflen);
	      Function *coeff = ReadFunction(input_file);
	      for(int k=0; k<bnv; k++)
		{
		  int indie = mesh->GetBdryGindex(b, k);
		  double* coord = mesh->GetVertex(indie);

		  elrbdryval[bb][k] = bdryval->Eval(coord);
		  elrcoeff[bb][k] = coeff->Eval(coord);

		  eldbdryval[bb][k] = 0.0;
		  elnbdryval[bb][k] = 0.0;
		}
	      delete bdryval;
	      delete coeff;
	    }
      
	}
    }

  //set elasticity data
  Data *elasticdata;
  elasticdata = new Data();
  elasticdata->SetNumberofBoundaryConditions(2);
  elasticdata->SetEllipticData(lamecoeff);
  elasticdata->SetNeumannData(elnbdryval, elBdrNeumann);
  elasticdata->SetRobinData(elrcoeff, elrbdryval, elBdrRobin);
  elasticdata->SetDirichletData(eldbdryval, elBdrDirichlet);

  cout << "Elastic Data Set" << endl;

  //---------------------------------------------------------
  //  Set initial profiles
  //---------------------------------------------------------

  //Set initial displacment and pressure
  ReadIdentifier(input_file, buffer, bufflen);
  Function *initu1 = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *initu2 = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *initp = ReadFunction(input_file);

  Vector displacementinit(2*nv);
  Vector pressureinit(nv);
  for(int k=0; k<nv; k++)
    {
      double* coord = mesh->GetVertex(k);
      displacementinit(k) = initu1->Eval(coord);
      displacementinit(k+nv) = initu2->Eval(coord);
      pressureinit(k) = initp->Eval(coord);
    }
  delete initu1;
  delete initu2;
  delete initp;

  
  //---------------------------------------------------------
  //  Set Exact Solutions
  //---------------------------------------------------------
  
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solp = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solu1 = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solu2 = ReadFunction(input_file);

  //---------------------------------------------------------
  //  Build Global Matrix
  //---------------------------------------------------------
  SparseMatrix *globalmat;
  globalmat = new SparseMatrix(3*nv, 50);
  Vector globalrhs(3*nv);
  globalrhs = 0.0;

  //Fill Ae to Global Matrix
  FEM *elasticfem;
  if(meshtype==Element::QUADRILATERAL){ elasticfem = new BilinearFEMELAST2D(mesh, elasticdata); }
  if(meshtype==Element::TRIANGLE){ elasticfem = new LinearFEMELAST2D(mesh, elasticdata); }
  elasticfem->AssembleMatrix();
  elasticfem->FinalizeMatrix();
  SparseMatrix elastmatrix = elasticfem->getMatrix();

  int *I = elastmatrix.GetI();
  int *J = elastmatrix.GetJ();
  for(int k=0; k<2*nv; k++)
    {
       for(int i=I[k]; i<I[k+1]; i++)
	{
	  globalmat->Elem(k, J[i]) = elastmatrix(k, J[i]);
	}
    }
  
  //Fill Ap to Global Matrix  
  FEM *pressurefem;
  if(meshtype==Element::QUADRILATERAL){ pressurefem = new BilinearFEM2D(mesh, pressuredata); }
  if(meshtype==Element::TRIANGLE){ pressurefem = new LinearFEM2D(mesh, pressuredata); }
  pressurefem->AssembleMatrix(); 
  pressurefem->FinalizeMatrix();
  SparseMatrix pressmatrix = pressurefem->getMatrix();

  I = pressmatrix.GetI();
  J = pressmatrix.GetJ();
  for(int k=0; k<nv; k++)
    {
      for(int i=I[k]; i<I[k+1]; i++)
	{
	  globalmat->Elem(k+2*nv, J[i]+2*nv) = pressmatrix(k, J[i]);
	}
    }

  
  Array<double *> corrector(1);
  corrector[0] = new double[nv];
  for(int i=0; i<nv; i++)
    {
      corrector[0][i] = 0.0; //h*h/(4*dt*(lamecoeff[0][i]+lamecoeff[1][i]));
    }
  Data *betadata;
  betadata = new Data();
  betadata->SetEllipticData(corrector);

  FEM *betafem;
  if(meshtype==Element::QUADRILATERAL){ betafem = new BilinearFEM2D(mesh, betadata); }
  if(meshtype==Element::TRIANGLE){ betafem = new LinearFEM2D(mesh, betadata); }
  betafem->AssembleMatrix();
  betafem->FinalizeMatrix(); 
  SparseMatrix betamatrix = betafem->getMatrix();
  
  I = pressmatrix.GetI();
  J = pressmatrix.GetJ();
  for(int k=0; k<nv; k++)
    {
      for(int i=I[k]; i<I[k+1]; i++)
	{
	  //globalmat->Elem(k+2*nv, J[i]+2*nv) += betamatrix(k, J[i]);
	}
    }
  
  // cout << endl;
  
  //Set off diagonal data
  Data *scalardata;
  scalardata = new Data();

  Array<double *> scalarcoeff(2);
  scalarcoeff[0] = new double[nv];
  scalarcoeff[1] = new double[nv];
  for(int i=0; i<nv; i++)
    {
      scalarcoeff[0][i] = 1.0; 
      scalarcoeff[1][i] = 0.0; 
    }
  scalardata->SetAdvectionData(scalarcoeff);

  FEM *scalarfemx;
  if(meshtype==Element::QUADRILATERAL){ scalarfemx = new BilinearFEM2D(mesh, scalardata); }
  if(meshtype==Element::TRIANGLE){ scalarfemx = new LinearFEM2D(mesh, scalardata); }
  scalarfemx->AssembleMatrix();
  scalarfemx->FinalizeMatrix(); 
  SparseMatrix scalarxmatrix = scalarfemx->getMatrix();

  I = scalarxmatrix.GetI();
  J = scalarxmatrix.GetJ();
  for(int k=0; k<nv; k++)
    {
       for(int i=I[k]; i<I[k+1]; i++)
	{
	  globalmat->Elem(k, J[i]+2*nv) = scalarxmatrix(k, J[i]);
	  globalmat->Elem(k+2*nv, J[i]) = scalarxmatrix(k, J[i])/dt;
	}
    }

  Data *scalardatay = new Data();
  Array<double *> scalarcoeffy(2);
  scalarcoeffy[0] = new double[nv];
  scalarcoeffy[1] = new double[nv];

  for(int i=0; i<nv; i++)
    {
      scalarcoeffy[0][i] = 0.0; 
      scalarcoeffy[1][i] = 1.0; 
    }
  scalardatay->SetAdvectionData(scalarcoeffy);

  FEM *scalarfemy;
  if(meshtype==Element::QUADRILATERAL){ scalarfemy = new BilinearFEM2D(mesh, scalardatay); }
  if(meshtype==Element::TRIANGLE){ scalarfemy = new LinearFEM2D(mesh, scalardatay); }
  scalarfemy->AssembleMatrix(); 
  scalarfemy->FinalizeMatrix(); 
  SparseMatrix scalarymatrix = scalarfemy->getMatrix();

  I = scalarymatrix.GetI();
  J = scalarymatrix.GetJ();
  for(int k=0; k<nv; k++)
    {
       for(int i=I[k]; i<I[k+1]; i++)
	{
	  globalmat->Elem(k+nv, J[i]+2*nv) = scalarymatrix(k, J[i]);
	  globalmat->Elem(k+2*nv, J[i]+nv) = scalarymatrix(k, J[i])/dt;
	}
    }
  
  
  //Dirichlet Boundaries Elastic
  if (elasticdata->BdrDirichletExists())
    {
      for(int j=0; j<elasticdata->GetNumberofBoundaryConditions(); j++)
	{
	  for(int b=0; b<nbdry; b++)
	    {
	      int bb = b + j*nbdry;
	      if(elasticdata->BdryDirichlet(bb))
		{
		  int numel = mesh->GetNBdryElem(b);
		  for(int i=0; i<=numel; i++)
		    {
		      int index = mesh->GetBdryGindex(b, i) + j*mesh->GetNV();
		      double dval = elasticdata->GetDirichletBdryVal(bb, i);
		      globalmat->ClearRow(index, dval, globalrhs); 
		    }	   
		}
	    }
	}
    }

  //Dirichlet Boundaries Pressure
  if (pressuredata->BdrDirichletExists())
    {
      for(int b=0; b<nbdry; b++)
	{
	  if(pressuredata->BdryDirichlet(b))
	    {
	      int numel = mesh->GetNBdryElem(b);
	      for(int i=0; i<=numel; i++)
		{
		  int index = mesh->GetBdryGindex(b, i)+2*nv;
		  double dval = pressuredata->GetDirichletBdryVal(b, i);
		  globalmat->ClearRow(index, dval, globalrhs); 
		}	   
	    }
	}
    }

  globalmat->Finalize();

  Vector u1(nv);
  Vector u2(nv);
  Vector pressure(nv);
  Vector displacement(2*nv);
  for(int i=0; i<nv; i++){u1(i) = displacementinit(i); u2(i) = displacementinit(i+nv); pressure(i) = pressureinit(i);}
  for(int i=0; i<2*nv; i++){displacement(i) = displacementinit(i);}

 //---------------------------------------------------------
  //  Build Global Matrix
  //---------------------------------------------------------

  for(int n=1; n<=numfinesteps; n++)
    {
      char filename[256];
      sprintf(filename, "pressure_%f.vtk", (n-1)*dt);
      ofstream output;
      output.open(filename);
      mesh->ParaviewPrint(output, pressure, meshtype); 
      output.close();  
  
      sprintf(filename, "u1_%f.vtk",(n-1)*dt);
      output.open(filename);
      mesh->ParaviewPrint(output, u1, meshtype); 
      output.close(); 

      sprintf(filename, "u2_%f.vtk", (n-1)*dt);
      output.open(filename);
      mesh->ParaviewPrint(output, u2, meshtype); 
      output.close();  
  
      sprintf(filename, "defmesh_%f.vtk", (n-1)*dt);
      output.open(filename);
      mesh->ParaviewPrintDeformedMesh(output, displacement, meshtype);
      output.close();     


      globalrhs = 0.0;
      //Add F1 and F2 to rhs
      Array<double *> elfdata(2);
      elfdata[0] = new double[nv];
      elfdata[1] = new double[nv];

      for(int k=0; k<nv; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  double coords[3];
	  coords[0] = coord[0];
	  coords[1] = coord[1];
	  coords[2] = n*dt;
	  elfdata[0][k] = forcex->Eval(coords);
	  elfdata[1][k] = forcey->Eval(coords);
	}

      Data *elastdata;
      elastdata = new Data();
      elasticdata->SetNumberofBoundaryConditions(2);
      elastdata->SetForceData(elfdata);
      
      FEM *elastfem;
      if(meshtype==Element::QUADRILATERAL){ elastfem = new BilinearFEMELAST2D(mesh, elastdata); }
      if(meshtype==Element::TRIANGLE){ elastfem = new LinearFEMELAST2D(mesh, elastdata); }
      elastfem->AssembleMatrix();
      Vector urhs(2*nv);
      elastfem->getRHS(urhs);

      for(int j=0; j<2*nv; j++)
	{
	  globalrhs(j) += urhs(j);
	}

      //Add divu and F3 to rhs
      Array<double *> forcedata(1);
      forcedata[0] = new double[nv];

      for(int k=0; k<nv; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  double coords[3];
	  coords[0] = coord[0];
	  coords[1] = coord[1];
	  coords[2] = n*dt;
	  forcedata[0][k] = force->Eval(coords);
	}

      Data *pressdata;
      pressdata = new Data();
      pressdata->SetForceData(forcedata);
  
      FEM *pressfem;
      if(meshtype==Element::QUADRILATERAL){ pressfem = new BilinearFEM2D(mesh, pressdata); }
      if(meshtype==Element::TRIANGLE){ pressfem = new LinearFEM2D(mesh, pressdata); }
      pressfem->AssembleMatrix(); 
      Vector prhs(nv);
      pressfem->getRHS(prhs);

      for(int j=0; j<nv; j++)
	{
	  globalrhs(j+2*nv) += prhs(j);
	}

     Vector Qu1(nv);
     Vector Ru2(nv);
     scalarxmatrix.Mult(u1, Qu1);
     scalarymatrix.Mult(u2, Ru2);

     for(int j=0; j<nv; j++)
       {
	 globalrhs(j+2*nv) += (Qu1(j)+Ru2(j))/dt;
       }

      
      //Add Corrector to RHS 
      Vector betap0(nv);
      betamatrix.Mult(pressure, betap0);
      for(int k=0; k<nv; k++)
	{
	  //globalrhs(k+2*nv) += betap0(k);
	}
      

      //Dirichlet Boundaries Elastic
      if (elasticdata->BdrDirichletExists())
	{
	  for(int j=0; j<elasticdata->GetNumberofBoundaryConditions(); j++)
	    {
	      for(int b=0; b<nbdry; b++)
		{
		  int bb = b + j*nbdry;
		  if(elasticdata->BdryDirichlet(bb))
		    {
		      int numel = mesh->GetNBdryElem(b);
		      for(int i=0; i<=numel; i++)
			{
			  int index = mesh->GetBdryGindex(b, i) + j*mesh->GetNV();
			  double dval = elasticdata->GetDirichletBdryVal(bb, i);
			  globalrhs(index) = dval; 
			}	   
		    }
		}
	    }
	}

      //Dirichlet Boundaries Pressure
      if (pressuredata->BdrDirichletExists())
	{
	  for(int b=0; b<nbdry; b++)
	    {
	      if(pressuredata->BdryDirichlet(b))
		{
		  int numel = mesh->GetNBdryElem(b);
		  for(int i=0; i<=numel; i++)
		    {
		      int index = mesh->GetBdryGindex(b, i)+2*nv;
		      double dval = pressuredata->GetDirichletBdryVal(b, i);
		      globalrhs(index) = dval; 
		    }	   
		}
	    }
	}

      Vector sol(3*nv);
      sol = 0.0;

      int mem_fact = 100;
      MultiFrontalMatrixInverse *InvMat = new MultiFrontalMatrixInverse(*globalmat, mem_fact);
      InvMat->Mult(globalrhs, sol);
      delete InvMat;

      Vector displacementold(2*nv);
      displacementold = displacement;
      Vector u1old(nv);
      Vector u2old(nv);
      u1old = u1;
      u2old = u2;

      for(int i=0; i<nv; i++){u1(i) = sol(i); u2(i) = sol(i+nv); pressure(i) = sol(i+2*nv);}
      for(int i=0; i<2*nv; i++){displacement(i) = sol(i);}

      cout << "solved" << endl;
      
      //sol.Print();

      // Post-processing for fluxes
      cout << endl << "Post-processing..." << endl;
      DualMesh *dualmesh = new DualMesh(mesh);

      Array<double *> flux(dualmesh->GetNumLocalDOF());
      for(int i=0; i<flux.Size(); i++ ) flux[i] = new double[dualmesh->GetNE()];

      // Collect Q, F
      Array<double *> QF(dualmesh->GetNumLocalDOF());
      for(int i=0; i<QF.Size(); i++ ) QF[i] = new double[dualmesh->GetNE()];

      for(int i=0; i<dualmesh->GetNE(); i++)
	{
	  Array<int> ind;
	  dualmesh->GetDOF(i, ind);

	  int locdof = dualmesh->GetNumLocalDOF();

	  Vector Q(locdof);
	  Q = 0.0;
	  Array<double *> tempmat(locdof);
	  for(int k=0; k<locdof; k++) tempmat[k] = new double[locdof];
	  
	  pressurefem->GetEllipticLocalSystem(i, ind, tempmat);
	  
	  for(int k=0; k<locdof; k++)
	    for(int j=0; j<locdof; j++)
	      Q(k) += tempmat[k][j] * pressure(ind[j]);

	  Array<double> F(locdof);
	  pressfem->GetForceLocalSystem(i, ind, F);

	  // Local scalarfemx mult u1
	  scalarfemx->GetAdvectionLocalSystem(i, ind, tempmat);
	  
	  for(int k=0; k<locdof; k++)
	    for(int j=0; j<locdof; j++)
	      Q(k) += tempmat[k][j] * ( u1(ind[j]) - u1old(ind[j]) ) / dt;
	  
	  // Local scalarfemy mult u2
	  scalarfemy->GetAdvectionLocalSystem(i, ind, tempmat);
	  
	  for(int k=0; k<locdof; k++)
	    for(int j=0; j<locdof; j++)
	      Q(k) += tempmat[k][j] * ( u2(ind[j]) - u2old(ind[j]) ) / dt;

	  for(int k=0; k<locdof; k++) QF[k][i] = Q(k) - F[k];

	  for(int k=0; k<tempmat.Size(); k++) delete []tempmat[k];
	}

      // Gathering global neumann data locally
      Array<double *> ndata(dualmesh->GetNumLocalDOF());
      for(int i=0; i<ndata.Size(); i++ ) ndata[i] = new double[dualmesh->GetNE()];
      for(int i=0; i<ndata.Size(); i++ ) 
	for(int k=0; k<dualmesh->GetNE(); k++)
	  ndata[i][k] = 0.0;

      int bdryindex = 0;
      if(pressuredata->BdrNeumannExists())
	for(int b=0; b<mesh->GetNBdrs(); b++)
	  {
	    // Go to each boundary of the domain
	    if(pressuredata->BdryNeumann(b))
	      {
		// Go to Bdryelement calculation if its a Neumann boundary 
		for(int i=0; i<mesh->GetNBdryElem(b); i++)
		  {
		    Array<int> ind;
		    mesh->GetBdrElementVertices(bdryindex+i, ind);
		    
		    // Figure out which element the bdryelement belongs to
		    int cvne = dualmesh->Getcvnumel(ind[1]);
		    Array<int> cv2el(cvne);
		    for(int k=0; k<cvne; k++)
		      cv2el[k] = dualmesh->Getcvelindex(ind[1], k);

		    int elnum = -1;
		    for(int k=0; k<cvne; k++)
		      {
			Array<int> tempind;
			mesh->GetElementVertices(cv2el[k], tempind);
			for(int i=0; i<tempind.Size(); i++)
			  {
			    if(ind[0]==tempind[i])
			      elnum = cv2el[k];
			  }
		      }
		    if(elnum == -1) cout << "Element associated the boundary is not found" << endl;
		    
		    // Get the local Neumann data
		    Array<double> locndata(2);
		    pressfem->GetNeumannLocalSystem(0, bdryindex+i, ind, locndata);
		    
		    if(b<2)
		      {
			ndata[0][elnum] = locndata[0];
			ndata[1][elnum] = locndata[1];
		      }
		    else if(b==2)
		      {
			ndata[1][elnum] = locndata[0];
			ndata[2][elnum] = locndata[1];
		      }
		    else
		      {
			ndata[2][elnum] = locndata[0];
			ndata[0][elnum] = locndata[1];
		      }
		  }
	      }
	    bdryindex += mesh->GetNBdryElem(b);
	  }

      // Accumulating the global neumann data to F
      for(int i=0; i<QF.Size(); i++ ) 
	for(int k=0; k<dualmesh->GetNE(); k++)
	  QF[i][k] += ndata[i][k];
      

      // Calculate divut
      Array<double *> divut(dualmesh->GetNumLocalDOF());
      for(int i=0; i<divut.Size(); i++ ) divut[i] = new double[dualmesh->GetNE()];

      Array<double *> divu(dualmesh->GetNumLocalDOF());
      for(int i=0; i<divu.Size(); i++ ) divu[i] = new double[dualmesh->GetNE()];
      Array<double *> divuold(dualmesh->GetNumLocalDOF());
      for(int i=0; i<divuold.Size(); i++ ) divuold[i] = new double[dualmesh->GetNE()];

      calculatedivu(mesh, displacementold, divuold);
      calculatedivu(mesh, displacement, divu);

      for(int i=0; i<dualmesh->GetNE(); i++)
	{
	  if(dualmesh->GetNumLocalDOF() == 3)
	    {
	      for(int k=0; k<dualmesh->GetNumLocalDOF(); k++)
		divut[k][i] = ( divu[k][i] - divuold[k][i] ) / dt;
	    }
	  else if(dualmesh->GetNumLocalDOF() == 4)
	    {
	      Array<int> ind;
	      dualmesh->GetDOF(i, ind);

	      double tx[4];
	      double ty[4];
	      for(int k=0; k<dualmesh->GetNumLocalDOF(); k++)
		{
		  tx[k] = ( u1(ind[k]) - u1old(ind[k]) ) / ( 16.0*dt );
		  ty[k] = ( u2(ind[k]) - u2old(ind[k]) ) / ( 16.0*dt );
		}

	      divut[0][i] = ( 3.0*tx[1] - 3.0*tx[0] + tx[2] - tx[3] ) * hy + ( 3.0*ty[3] - 3.0*ty[0] + ty[2] - ty[1] ) * hx;
	      divut[1][i] = ( 3.0*tx[1] - 3.0*tx[0] + tx[2] - tx[3] ) * hy + ( 3.0*ty[2] - 3.0*ty[1] + ty[3] - ty[0] ) * hx;
	      divut[2][i] = ( 3.0*tx[2] - 3.0*tx[3] + tx[1] - tx[0] ) * hy + ( 3.0*ty[2] - 3.0*ty[1] + ty[3] - ty[0] ) * hx;
	      divut[3][i] = ( 3.0*tx[2] - 3.0*tx[3] + tx[1] - tx[0] ) * hy + ( 3.0*ty[3] - 3.0*ty[0] + ty[2] - ty[1] ) * hx;
	    }
	  else { cout << "Mesh Error: Post-processing is valid for this mesh." << endl; }
	}

      ComputeConservativeFlux(dualmesh, permdata, divut, forcedata, QF, pressure, flux);

      // Calculate the post-processed L2 norm error
      double perr = ComputeDiscL2Error(mesh, solp, flux);
      cout << "The L2 error for the post-processed pressure is: " << perr << endl;

      ReadIdentifier(input_file, buffer, bufflen);
      Function *solpx = ReadFunction(input_file);
      ReadIdentifier(input_file, buffer, bufflen);
      Function *solpy = ReadFunction(input_file);
      double pherr = ComputeDiscH1SemiNormError(mesh, solpx, solpy, flux);
      cout << "The H1 SemiNorm error for the post-processed pressure is: " << pherr << endl;

      double preherr = ComputeH1SemiNormError(mesh, solpx, solpy, pressure);
      cout << "The H1 SemiNorm error for the pressure is: " << preherr << endl;

      cout << "Post-processing for flux is done." << endl << endl;

      // Post-processing for displacements
      Vector ut(2*nv);
      for(int i=0; i<2*nv; i++) ut(i) = ( displacement(i) - displacementold(i) ) / dt;

      Array<double *> uhdotn(dualmesh->GetNumLocalDOF());
      for(int i=0; i<uhdotn.Size(); i++ ) uhdotn[i] = new double[dualmesh->GetNE()];

      Array<double> buhdotn(nv);

      Computeuhdotn(dualmesh, ut, uhdotn, buhdotn);

      cout << "Post-processing for uhdotn is done." << endl << endl;
      
      cout << "Local Conservation Checking..." << endl;

      Vector cverr(nv);
      cverr = 0.0;
      for(int k=0; k<mesh->GetNE(); k++)
	{
	  Array<int> ind;
	  dualmesh->GetDOF(k, ind);

	  if(dualmesh->GetNumLocalDOF()==3)
	    {
	      cverr(ind[0]) += flux[0][k] - flux[2][k] + uhdotn[0][k] - uhdotn[2][k];
	      cverr(ind[1]) += flux[1][k] - flux[0][k] + uhdotn[1][k] - uhdotn[0][k];
	      cverr(ind[2]) += flux[2][k] - flux[1][k] + uhdotn[2][k] - uhdotn[1][k];
	    }
	  else if(dualmesh->GetNumLocalDOF()==4)
	    {
	      cverr(ind[0]) += flux[0][k] - flux[3][k] + uhdotn[0][k] - uhdotn[3][k];
	      cverr(ind[1]) += flux[1][k] - flux[0][k] + uhdotn[1][k] - uhdotn[0][k];
	      cverr(ind[2]) += flux[2][k] - flux[1][k] + uhdotn[2][k] - uhdotn[1][k];
	      cverr(ind[3]) += flux[3][k] - flux[2][k] + uhdotn[3][k] - uhdotn[2][k];
	    }
	  else { cout << "Mesh Error: Post-processing is valid for this mesh." << endl; }
	}

      cout << "Local conservation error for displacement" << endl;
      //for(int k=0; k<nv; k++) cout << k << "\t" << cverr(k) << endl;

      for(int k=0; k<ndata.Size(); k++) delete []ndata[k];
      for(int k=0; k<uhdotn.Size(); k++) delete []uhdotn[k];
      for(int k=0; k<flux.Size(); k++) delete []flux[k];
      for(int k=0; k<QF.Size(); k++) delete []QF[k];
      for(int k=0; k<divut.Size(); k++) delete []divut[k];
      for(int k=0; k<divu.Size(); k++) delete []divu[k];
      for(int k=0; k<divuold.Size(); k++) delete []divuold[k];
      
      delete dualmesh;
      delete elastfem;
      delete pressfem;
      delete elastdata;
      delete pressdata;
      delete []forcedata[0];
      for(int k=0; k<elfdata.Size(); k++)
	delete []elfdata[k]; 	  
    }
  
  //Estimate errors
  cout << "L2error pressure: " << pressurefem->ComputeL2Error(solp, pressure) << endl;	 
  cout << "L2error u1: " << elasticfem->ComputeL2Error(solu1, u1) << endl;
  cout << "L2error u2: " << elasticfem->ComputeL2Error(solu2, u2) << endl;

  DualMesh *dualmesh = new DualMesh(mesh);
  delete dualmesh;
  //---------------------------------------------------------
  //  Freedom!!!!
  //---------------------------------------------------------
  delete betafem;
  delete betadata;
  delete force;
  delete forcex;
  delete forcey;
  delete elasticfem;
  delete pressurefem;
  delete solp;
  delete solu1;
  delete solu2;
  delete globalmat;

  delete mesh;
  delete pressuredata;
  delete elasticdata;
  delete scalardata;
  delete scalardatay;
  delete scalarfemx;
  delete scalarfemy;
  delete []permdata[0];
  delete []corrector[0];

  for(int k=0; k<scalarcoeffy.Size(); k++)
    delete []scalarcoeffy[k];
  
  for(int k=0; k<scalarcoeff.Size(); k++)
    delete []scalarcoeff[k];

  for(int k=0; k<lamecoeff.Size(); k++)
    delete []lamecoeff[k];

  for(int k=0; k<dbdryval.Size(); k++)
    delete []dbdryval[k];

  for(int k=0; k<nbdryval.Size(); k++)
    delete []nbdryval[k];

  for(int k=0; k<rbdryval.Size(); k++)
    delete []rbdryval[k];

  for(int k=0; k<rcoeff.Size(); k++)
    delete []rcoeff[k];

  for(int k=0; k<eldbdryval.Size(); k++)
    delete []eldbdryval[k];

  for(int k=0; k<elnbdryval.Size(); k++)
    delete []elnbdryval[k];

  for(int k=0; k<elrbdryval.Size(); k++)
    delete []elrbdryval[k];

  for(int k=0; k<elrcoeff.Size(); k++)
    delete []elrcoeff[k];


  cout << endl << "fin" << endl;
}
