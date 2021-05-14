#include <iostream>
#include <fstream>
#include <cmath>
#include "../linalg/linalg_header.h"
#include "../mesh/mesh_header.h"
#include "../fem/fem_header.h"
#include "../timedependent/timedependent_header.h"

using namespace std;

//=====================================================================================================
void GetPPDarcy(DualMesh *dualmesh, Data *pressuredata, Vector &sol, Vector &solold, FEM *pressurefem, 
		FEM *scalarfemx, FEM *scalarfemy, FEM *betafem, double dt, Array<double *> &flux, Array<double *> &ppsol)
{
  cout << endl;
  int nlocdof = dualmesh->GetNumLocalDOF();
  int ne = dualmesh->GetNE();
  int nv = dualmesh->GetNV();

  flux.SetSize(nlocdof);
  ppsol.SetSize(nlocdof);
  Array<double*> ndata(nlocdof);
  for(int i=0; i<nlocdof; i++)
    {
      flux[i] = new double[ne];
      ppsol[i] = new double[ne]; 
      ndata[i] = new double[ne];
      for(int k=0; k<ne; k++)
	{
	  flux[i][k] = 0.0;
	  ppsol[i][k] = 0.0;
	  ndata[i][k] = 0.0;
	}
    }

  // Gathering global neumann data locally
  int bdryelindex = 0;
  if(pressuredata->BdrNeumannExists())
    {
      for(int b=0; b<dualmesh->GetNBdrs(); b++)
	{
	  // Go to each boundary of the domain
	  if(pressuredata->BdryNeumann(b))
	    {
	      // Go to Bdryelement calculation if its a Neumann boundary 
	      for(int i=0; i<dualmesh->GetNBdryElem(b); i++)
		{
		  Array<int> ind;
		  dualmesh->GetBdrElementVertices(bdryelindex+i, ind);
		  
		  // Figure out which element the bdryelement belongs to
		  int cvne0 = dualmesh->Getcvnumel(ind[0]);
		  int cvne1 = dualmesh->Getcvnumel(ind[1]);
		  Array<int> cvelind0(cvne0);
		  Array<int> cvelind1(cvne1);
		  
		  for(int k=0; k<cvne0; k++)
		    {
		      cvelind0[k] = dualmesh->Getcvelindex(ind[0], k);
		    }
		  for(int k=0; k<dualmesh->Getcvnumel(ind[1]); k++)
		    {
		      cvelind1[k] = dualmesh->Getcvelindex(ind[1], k);
		    }
		  
		  int el;
		  int mincvne = (cvne0 < cvne1) ? cvne0 : cvne1;
		  for(int k=0; k<mincvne; k++)
		    {
		      if(cvelind0[k] == cvelind1[k])
			{
			  el = cvelind0[k];
			  break;
			}
		    }
		  
		  // Get the local Neumann data
		  Array<double> locndata(2);
		  Array<int> localindex(2);
		  localindex[0] = i;
		  localindex[1] = i+1;
		  pressurefem->GetNeumannLocalSystem(b, i, localindex, locndata);	
		  
		  if(b==0 && nlocdof==4)
		    {
		      ndata[0][el] += locndata[0];
		      ndata[1][el] += locndata[1];
		    }
		  else if(b==1 && nlocdof==4)
		    {
		      ndata[1][el] += locndata[0];
		      ndata[2][el] += locndata[1];
		    }
		  else if(b==2 && nlocdof==4)
		    {
		      ndata[3][el] += locndata[0];
		      ndata[2][el] += locndata[1];
		    }
		  else if(b==3 && nlocdof==4)
		    {
		      ndata[0][el] += locndata[0];
		      ndata[3][el] += locndata[1];
		    }		
		  else if(b<2 && nlocdof==3)
		    {
		      ndata[0][el] += locndata[0];
		      ndata[1][el] += locndata[1];
		    }
		  else if(b==2  && nlocdof==3)
		    {
		      ndata[1][el] += locndata[1];
		      ndata[2][el] += locndata[0];
		    }
		  else if(b==3 && nlocdof==3)
		    {
		      ndata[2][el] += locndata[1];
		      ndata[0][el] += locndata[0];
		    }	
		}
	      bdryelindex += dualmesh->GetNBdryElem(b);
	    }
	}
    }
      

  for(int i=0; i<ne; i++)
    {
      Array<double *> normals;
      Array<double> elengths;
      Array<double> areas;     
      dualmesh->GetNormals(i, normals);
      dualmesh->GetEdgeLengths(i, elengths);
      dualmesh->GetAreas(i, areas);

      //=============== Build Local RHS ==========
      Array<int> ind;
      dualmesh->GetDOF(i, ind);

      Vector test(nlocdof);
      test = 0.0;

      Vector locrhs(nlocdof);
    
      locrhs = 0.0;
      Array<double *> tempmat(nlocdof);
      for(int k=0; k<nlocdof; k++)
	{
	  tempmat[k] = new double[nlocdof];
	}

      Vector locp(nlocdof);
      Vector locpt(nlocdof);
      Vector locu1t(nlocdof);
      Vector locu2t(nlocdof);
      double tau = 1.0/((double) dt);
      for(int k=0; k<nlocdof; k++)
	{
	  locu1t(k) = tau*(sol(ind[k])-solold(ind[k]));
	  locu2t(k) = tau*(sol(ind[k]+nv)-solold(ind[k]+nv));
	  locpt(k) = (sol(ind[k]+2*nv)-solold(ind[k]+2*nv));
	  locp(k) = sol(ind[k]+2*nv);
	}

      pressurefem->GetEllipticLocalSystem(i, ind, tempmat);
      for(int k=0; k<nlocdof; k++)
	{
	  for(int j=0; j<nlocdof; j++)	      
	    {
	      locrhs(k) += tempmat[k][j] * locp(j);
	    }
	}

      Array<double> F(nlocdof);
      pressurefem->GetForceLocalSystem(i, ind, F);
      for(int k=0; k<nlocdof; k++)	      
	{
	  locrhs(k) -= F[k];
	  locrhs(k) += pressuredata->GetNodalForceCoeff(ind[k])*areas[k];
	}
    

      betafem->GetEllipticLocalSystem(i, ind, tempmat);
      for(int k=0; k<nlocdof; k++)
	{
	  for(int j=0; j<nlocdof; j++)	      
	    {
	      locrhs(k) += tempmat[k][j] * locpt(j);
	    }
	}
	 
      scalarfemx->GetAdvectionLocalSystem(i, ind, tempmat);	  
      for(int k=0; k<nlocdof; k++)
	{
	  for(int j=0; j<nlocdof; j++)
	    {
	      locrhs(k) += tempmat[k][j] * locu1t(j);
	      test(k) += tempmat[k][j] * locu1t(j);
	    }
	}

      scalarfemy->GetAdvectionLocalSystem(i, ind, tempmat);	  
      for(int k=0; k<nlocdof; k++)
	{
	  for(int j=0; j<nlocdof; j++)
	    {
	      locrhs(k) += tempmat[k][j] * locu2t(j);
	      test(k) += tempmat[k][j] * locu2t(j);
	    }
	}

      for(int k=0; k<nlocdof; k++)
	{
	  locrhs(k) += ndata[k][i];
	}

 


      Array<double *> coord;
      dualmesh->GetElementVerticesCoord(i, coord);
      if(nlocdof==3)
	{
	  double x0 = coord[0][0];
	  double y0 = coord[0][1];
	  double x1 = coord[1][0];
	  double y1 = coord[1][1];
	  double x2 = coord[2][0];
	  double y2 = coord[2][1];
	  double det = -x1*y0 + x2*y0 + x0*y1 - x2*y1 - x0*y2 + x1*y2; 
	  
	  double du1tx = (locu1t(0)*(y1-y2) + locu1t(1)*(y2-y0) + locu1t(2)*(y0-y1))/det;
	  double du2ty = (locu2t(0)*(x2-x1) + locu2t(1)*(x0-x2) + locu2t(2)*(x1-x0))/det;
	  double dut = du1tx + du2ty;

	  for(int k=0; k<nlocdof; k++)
	    {
	      locrhs(k) -= dut * areas[k];
	      test(k) -= dut * areas[k];
	    }
	}
      else if(nlocdof == 4)
	{

	  double dx = coord[1][0] - coord[0][0];
	  double dy = coord[2][1] - coord[0][1];
	  locrhs(0) -= ((-3*locu2t(0) - locu2t(1) + locu2t(2) + 3*locu2t(3))*dx + (-3*locu1t(0) + 3*locu1t(1) + locu1t(2) - locu1t(3))*dy)/16.;
	  locrhs(1) -= ((-locu2t(0) - 3*locu2t(1) + 3*locu2t(2) + locu2t(3))*dx + (-3*locu1t(0) + 3*locu1t(1) + locu1t(2) - locu1t(3))*dy)/16.;
	  locrhs(2) -= ((-locu2t(0) - 3*locu2t(1) + 3*locu2t(2) + locu2t(3))*dx + (-locu1t(0) + locu1t(1) + 3*locu1t(2) - 3*locu1t(3))*dy)/16.;
	  locrhs(3) -= ((-3*locu2t(0) - locu2t(1) + locu2t(2) + 3*locu2t(3))*dx + (-locu1t(0) + locu1t(1) + 3*locu1t(2) - 3*locu1t(3))*dy)/16.;

	  test(0) -= ((-3*locu2t(0) - locu2t(1) + locu2t(2) + 3*locu2t(3))*dx + (-3*locu1t(0) + 3*locu1t(1) + locu1t(2) - locu1t(3))*dy)/16.;
	  test(1) -= ((-locu2t(0) - 3*locu2t(1) + 3*locu2t(2) + locu2t(3))*dx + (-3*locu1t(0) + 3*locu1t(1) + locu1t(2) - locu1t(3))*dy)/16.;
	  test(2) -= ((-locu2t(0) - 3*locu2t(1) + 3*locu2t(2) + locu2t(3))*dx + (-locu1t(0) + locu1t(1) + 3*locu1t(2) - 3*locu1t(3))*dy)/16.;
	  test(3) -= ((-3*locu2t(0) - locu2t(1) + locu2t(2) + 3*locu2t(3))*dx + (-locu1t(0) + locu1t(1) + 3*locu1t(2) - 3*locu1t(3))*dy)/16.;
	}

      /*    
      double utestval=0.0;
      double rhstestval = 0;
      for(int j=0; j<nlocdof; j++)
	{
	  utestval += test(j);
	  rhstestval += locrhs(j);
	}
      cout << endl << i << " " << utestval << " " << rhstestval << endl;
      */
      
      //=============== Create Local Matrix ==========
      RectangularMatrix locmat(nlocdof);
        
      Array<double *> tem(nlocdof);
      for(int k=0; k<nlocdof; k++)
	{
	  tem[k] = new double[nlocdof];
	  for(int j=0; j<nlocdof; j++)
	    {
	      tem[k][j] = 0.0;
	    }
	}

      Array<double> lock(nlocdof);
      Array<double> locmk(nlocdof);
      for(int j=0; j<nlocdof; j++){ lock[j] = pressuredata->GetNodalEllipticCoeff(ind[j], 0); }
  
      if(nlocdof==3)
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

	  locmat(0,0) = tem[0][0] - tem[2][0];
	  locmat(0,1) = tem[0][1] - tem[2][1];
	  locmat(0,2) = tem[0][2] - tem[2][2];
      
	  locmat(1,0) = tem[1][0] - tem[0][0];
	  locmat(1,1) = tem[1][1] - tem[0][1];
	  locmat(1,2) = tem[1][2] - tem[0][2];
      
	  locmat(2,0) = tem[2][0] - tem[1][0];
	  locmat(2,1) = tem[2][1] - tem[1][1];
	  locmat(2,2) = tem[2][2] - tem[1][2];
	}
      else
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
      
	  locmat(0, 0) = tem[0][0] - tem[3][0];
	  locmat(0, 1) = tem[0][1] - tem[3][1];
	  locmat(0, 2) = tem[0][2] - tem[3][2];
	  locmat(0, 3) = tem[0][3] - tem[3][3];

	  locmat(1, 0) = tem[1][0] - tem[0][0];
	  locmat(1, 1) = tem[1][1] - tem[0][1];
	  locmat(1, 2) = tem[1][2] - tem[0][2];
	  locmat(1, 3) = tem[1][3] - tem[0][3];

	  locmat(2, 0) = tem[2][0] - tem[1][0];
	  locmat(2, 1) = tem[2][1] - tem[1][1];
	  locmat(2, 2) = tem[2][2] - tem[1][2];
	  locmat(2, 3) = tem[2][3] - tem[1][3];

	  locmat(3, 0) = tem[3][0] - tem[2][0];
	  locmat(3, 1) = tem[3][1] - tem[2][1];
	  locmat(3, 2) = tem[3][2] - tem[2][2];
	  locmat(3, 3) = tem[3][3] - tem[2][3];
	}

       locmat(0,0) += 1.0;  
      
      /*
      for(int k=1; k<nlocdof; k++) locmat(0, k) = 0.0;
      locmat(0,0) = 1.0;
      locrhs(0) = locp(0);
      */
      //=============== Solve Local System ===================
      Vector s(nlocdof);
      s = 0.0;
      DenseMatrixInverse invmat(locmat);
      invmat.Mult(locrhs, s);
      
      for(int k=0; k<nlocdof; k++)
	{      
	  ppsol[k][i] = s(k);
	  for(int j=0; j<nlocdof; j++)
	    {
	      //flux[k][i] += tem[k][j] * s(j);
	      flux[k][i] += tem[k][j] * locp(j);
	    }
	}
      
      for(int k=0; k<nlocdof; k++){ delete []tem[k]; }
      for(int k=0; k<nlocdof; k++){ delete []tempmat[k]; }  
    }
  
  /*
  cout << endl;
  for(int i=0; i<ne; i++)
    {
      cout << i << " ";
      for(int k=0; k<nlocdof; k++){ cout << flux[k][i] << "\t";  }
      cout << endl;
    }
  cout << endl;
  */

 for(int i=0; i<nlocdof; i++)
    {
      delete []ndata[i];
    }
 /*
 Vector cverr(nv);
 cverr = 0.0;
 for(int k=0; k<ne; k++)
   {
     Array<int> ind;
     dualmesh->GetDOF(k, ind);
    Array<double> areas;     
      dualmesh->GetAreas(k, areas);

     if(dualmesh->GetNumLocalDOF()==3)
       {
	 cverr(ind[0]) += flux[0][k] - flux[2][k];
	 cverr(ind[1]) += flux[1][k] - flux[0][k];
	 cverr(ind[2]) += flux[2][k] - flux[1][k];
       }
     else if(dualmesh->GetNumLocalDOF()==4)
       {
	 cverr(ind[0]) += flux[0][k] - flux[3][k];
	 cverr(ind[1]) += flux[1][k] - flux[0][k];
	 cverr(ind[2]) += flux[2][k] - flux[1][k];
	 cverr(ind[3]) += flux[3][k] - flux[2][k];

	 for(int j=0; j<nlocdof; j++){ cverr(ind[j]) -= pressuredata->GetNodalForceCoeff(ind[j])*areas[j]; }
       }
     else { cout << "Mesh Error: Post-processing is not valid for this mesh." << endl; }
   }

 for(int k=0; k<dualmesh->GetNV(); k++)
   {
     cout << k << " " << cverr(k) << endl;	
   } 
 */  
}


//==================================================================================================================

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
  int bdryelindex = 0;
  for(int b=0; b<dualmesh->GetNBdrs(); b++)
    {
      for(int i=0; i<dualmesh->GetNBdryElem(b); i++)
	{
	  Array<int> ind;
	  dualmesh->GetBdrElementVertices(bdryelindex+i, ind);

	  Array<double *> coord(2);
	  coord[0] = new double[2];
	  coord[1] = new double[2];
	  
	  dualmesh->GetBdrElementVerticesCoord(i, coord);
	  double lenx = (coord[0][0] - coord[0][1]);
	  double leny = (coord[1][0] - coord[1][1]);
	  double len = sqrt(lenx*lenx + leny*leny );
	  
	  if(b==0)
	    {
	      buhdotn[ind[0]] += -len * ( 3.0 * ut(ind[0]+NV) + ut(ind[1]+NV) ) / 8.0; 
	      buhdotn[ind[1]] += -len * ( 3.0 * ut(ind[1]+NV) + ut(ind[0]+NV) ) / 8.0; 
	    }
	  else if(b==1)
	    {
	      buhdotn[ind[0]] += len * ( 3.0 * ut(ind[0]) + ut(ind[1]) ) / 8.0; 
	      buhdotn[ind[1]] += len * ( 3.0 * ut(ind[1]) + ut(ind[0]) ) / 8.0; 
	    }
	  else if(b==2)
	    {
	      buhdotn[ind[0]] += len * ( 3.0 * ut(ind[0]+NV) + ut(ind[1]+NV) ) / 8.0; 
	      buhdotn[ind[1]] += len * ( 3.0 * ut(ind[1]+NV) + ut(ind[0]+NV) ) / 8.0; 
	    }
	  else if(b==3)
	    {
	      buhdotn[ind[0]] += -len * ( 3.0 * ut(ind[0]) + ut(ind[1]) ) / 8.0; 
	      buhdotn[ind[1]] += -len * ( 3.0 * ut(ind[1]) + ut(ind[0]) ) / 8.0; 
	    }
      
	  delete []coord[0];
	  delete []coord[1];
	}
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

double ComputeDiscH1SemiNormError(Mesh *mesh, Function *dirx, Function *diry, Array<double *> &sol, double currtime)
{ 
  double gp[8] = { 1.0/3.0, 1.0/3.0, 0.2, 0.6, 0.2, 0.2, 0.6, 0.2 };
  double gw[4] = { -27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0 };
  Array<double> x(3);
  Array<double> y(3);
  double coord[3];
  coord[2] = currtime;
   
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

double ComputeH1SemiNormError(Mesh *mesh, Function *dirx, Function *diry, Vector &sol, double currtime)
{ 
  double gp[8] = { 1.0/3.0, 1.0/3.0, 0.2, 0.6, 0.2, 0.2, 0.6, 0.2 };
  double gw[4] = { -27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0 };
  Array<double> x(3);
  Array<double> y(3);
  double coord[3];
  coord[2] = currtime;
   
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

void calculatedivu(Mesh *mesh, const Vector &u, Vector &divu)
{
  int nv = mesh->GetNV();
  int ne = mesh->GetNE();

  divu.SetSize(nv);
  divu = 0.0;

  Vector count(nv);
  count = 0.0;

  //Accumulate elemental divu value to each node
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

	  for(int i=0; i<ind.Size(); i++)
	    {
	      divu(ind[i]) += du;
	      count(ind[i]) += 1.0;
	    }
	}

      if(mesh->GetMType()==Element::QUADRILATERAL)
	{	
	  //Assuming rectangles	  
	  double dx = fabs(vcoords[0][1]-vcoords[0][0]);
	  double dy = fabs(vcoords[1][3]-vcoords[1][0]);

	  divu(ind[0]) += (alpha[1]-alpha[0])/dx;
	  divu(ind[1]) += (alpha[0]-alpha[1])/dx;
	  divu(ind[2]) += (alpha[3]-alpha[2])/dx;
	  divu(ind[3]) += (alpha[2]-alpha[3])/dx;

	  divu(ind[0]) += (beta[3]-beta[0])/dy;
	  divu(ind[1]) += (beta[2]-beta[1])/dy;
	  divu(ind[2]) += (beta[2]-beta[3])/dy;
	  divu(ind[3]) += (beta[0]-beta[0])/dy;

	  for(int i=0; i<ind.Size(); i++)
	    {
	      count(ind[i]) += 1.0;
	    }
	}
	 

      for(int i=0; i<vcoords.Size(); i++)
	delete []vcoords[i];
    }

  //Averageing
  for(int k=0; k<nv; k++){ divu(k) /= ((double) count(k)); }

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

  //--------------------------------------------------------------------------------
  // Set time data.
  //--------------------------------------------------------------------------------
  double finaltime;
  int numtimesteps;
  
  input_file.getline( buffer, bufflen, ' ' );
  input_file >> finaltime >> numtimesteps;
  double dt = finaltime/(numtimesteps);

  cout << endl << "Time Discretiazations: " << endl;
  cout << "Final Time: " << finaltime << endl;
  cout << "Number of Time Steps: " << numtimesteps << endl;
  cout << "Time Step Length: " << dt << endl; 
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

  double hx = (L[1]-L[0])/N[0];
  double hy = (L[3]-L[2])/N[1];
  double h = sqrt(hx*hx + hy*hy); 

  //print mesh
  mesh->ParaviewPrint(meshtype);
  int nv = mesh->GetNV(); //number of vertices

  //---------------------------------------------------------
  // Gather Data for Hydrodynamic Pressure Equation
  //---------------------------------------------------------
  //Set Permeability Data
  Array<double *> permdata(1);
  permdata[0] = new double[nv];

  char permflag[bufflen]; 
  ReadIdentifier(input_file, buffer, bufflen);
  input_file >> permflag;
  cout << permflag << endl;
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
  if(!strcmp(permflag, "data"))
    {
      char permfile[bufflen];
      double sigma;
      double a;
      ifstream permin;
      input_file >> permfile >> sigma >> a;
      permin.open(permfile);
      if(!permin) { cout << permfile << " does not exist.\n";  exit(1); }

      double val;
      double datamin = 10000.0;
      double datamax = 0.0;
      for (int k=0; k<nv; k++)
	{
	  permin >> val;
	  permdata[0][k] = a*exp(sigma * val);
	  if(permdata[0][k] > datamax){datamax = permdata[0][k]; }
	  if(permdata[0][k] < datamin){datamin = permdata[0][k]; }
	}
      cout << "Permeability Range: [" << datamin << ", " << datamax << "]" << endl;
      permin.close();
    }



  //Gather Force Data
  ReadIdentifier(input_file, buffer, bufflen);
  Function *force = ReadFunction(input_file);

  //Gather Boundary Data
  int nbdry = mesh->GetNBdrs();
  Array<double *> nbdryval(nbdry);
  Array<double *> dbdryval(nbdry);
  Array<bool> BdrNeumann(nbdry);
  Array<bool> BdrDirichlet(nbdry);
  Array<Function *> pressureboundaryfuncs(nbdry);

  char bdryflag[bufflen]; 
  for(int b=0; b<nbdry; b++)
    {
      int nvb = mesh->GetNBdryElem(b)+1;
      nbdryval[b] = new double[nvb];
      dbdryval[b] = new double[nvb];
      
      input_file.getline( buffer, bufflen, ' ' );
      input_file >> bdryflag;
      ReadIdentifier(input_file, buffer, bufflen);
      pressureboundaryfuncs[b] = ReadFunction(input_file);

      if (!strcmp(bdryflag, "neumann"))
	{
	  BdrNeumann[b] = true;
	  BdrDirichlet[b] = false;
	}
      if (!strcmp(bdryflag, "dirichlet"))
	{
	  BdrDirichlet[b] = true;
	  BdrNeumann[b] = false;
	}
      for(int k=0; k<nvb; k++)
	{
	  nbdryval[b][k] = 0.0;
	  dbdryval[b][k] = 0.0;
	}  
    }
  
  //set data for pressure equation
  Data *pressuredata;
  pressuredata = new Data();
  pressuredata->SetEllipticData(permdata);
  cout << "Pressure Data Set" << endl;


  //---------------------------------------------------------
  //  Gather Data for Poroelasticity Problem
  //--------------------------------------------------------- 

  int nnn = 2*nv; //number of dof
  //Gather elasticitiy Data
  Array<double *> lamecoeff(2);
  lamecoeff[0] = new double[nv];
  lamecoeff[1] = new double[nv];
  Vector elmod(nv);

  char elmodflag[bufflen]; 
  ReadIdentifier(input_file, buffer, bufflen);
  input_file >> elmodflag;
  double prat = 0.25;

  if (!strcmp(elmodflag, "equation"))
    {
      Function *modfunc;
      ReadIdentifier(input_file, buffer, bufflen);
      modfunc = ReadFunction(input_file);
      
      for(int k=0; k<nv; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  double data = modfunc->Eval(coord);
	  elmod(k) = log(data);
	  lamecoeff[0][k] = data*prat/((1.0+prat)*(1.0-2.0*prat));
	  lamecoeff[1][k] = data/(1+prat);
	}
      delete modfunc;
    }
  if(!strcmp(elmodflag, "data"))
    {
      char modfile[bufflen];
      double sigma;
      double a;
      ifstream modin;
      input_file >> modfile >> sigma >> a;
      modin.open(modfile);
      if(!modin) { cout << modfile << " does not exist.\n";  exit(1); }

      double val;
      double data;
      double datamin = a*1000.0;
      double datamax = 0.0;
      for (int k=0; k<nv; k++)
	{
	  modin >> val;
	  data = a*exp(sigma * val);
	  elmod(k) = log(data);
	  lamecoeff[0][k] = data*prat/((1.0+prat)*(1.0-2.0*prat));
	  lamecoeff[1][k] = data/(1+prat);
	  if(data > datamax){datamax = data; }
	  if(data < datamin){datamin = data; }
	}
      cout << "E Range: [" << datamin << ", " << datamax << "]" << endl;
      modin.close();
    }


  char filename[256];
  ofstream output;
  sprintf(filename, "elmod.vtk");
  output.open(filename);
  mesh->ParaviewPrint(output, elmod);
  output.close();

  cout << elmodflag << endl;
  ReadIdentifier(input_file, buffer, bufflen);
  Function *forcex = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *forcey = ReadFunction(input_file);

  //Set Boundary Data
  int nbredges = mesh->GetNBdrs();
  int nbdryel = 2*nbredges;
  Array<double *> elnbdryval(nbdryel);
  Array<double *> eldbdryval(nbdryel);
  Array<bool> elBdrNeumann(nbdryel);
  Array<bool> elBdrDirichlet(nbdryel);
  Array<Function *> elasticboundaryfuncs(nbdryel);

  for(int b=0; b<nbredges; b++)
    {
      int nvb = mesh->GetNBdryElem(b)+1;
      for(int i=0; i<2; i++)
	{
	  int bb = b + i*nbredges;
	  elnbdryval[bb] = new double[nvb];
	  eldbdryval[bb] = new double[nvb];
	  
	  input_file.getline( buffer, bufflen, ' ' );
	  input_file >> bdryflag;
	  ReadIdentifier(input_file, buffer, bufflen);
	  elasticboundaryfuncs[bb] = ReadFunction(input_file);
	  cout << bdryflag << endl;
	  if (!strcmp(bdryflag, "neumann"))
	    {
	      elBdrNeumann[bb] = true;
	      elBdrDirichlet[bb] = false; 
	    }
	  if (!strcmp(bdryflag, "dirichlet"))
	    {
	      elBdrDirichlet[bb] = true;
	      elBdrNeumann[bb] = false;
	    }
	  
	  for(int k=0; k<nvb; k++)
	    {
	      elnbdryval[bb][k] = 0.0;
	      eldbdryval[bb][k] = 0.0;
	    }  
	}
    }

  //set elasticity data
  Data *elasticdata;
  elasticdata = new Data();
  elasticdata->SetNumberofBoundaryConditions(2);
  elasticdata->SetEllipticData(lamecoeff);

  cout << "Elastic Data Set" << endl;   
  //---------------------------------------------------------
  //  Set Exact Solutions
  //---------------------------------------------------------
  
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solp = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solu1 = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solu2 = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solpx = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solpy = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solu1t = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solu2t = ReadFunction(input_file);

  cout <<buffer << endl;
  //---------------------------------------------------------
  //  Time March 
  //---------------------------------------------------------
  Vector u1(nv);
  Vector u1old(nv);
  Vector u2(nv);
  Vector u2old(nv);
  Vector pressure(nv);
  Vector pressureold(nv);
  Vector displacement(2*nv);
  Vector displacementold(2*nv);

  Vector sol(3*nv);
  Vector solold(3*nv);
  u1old = 0.0; 
  u2old = 0.0;
  pressureold = 0.0;
  displacementold = 0.0;
  sol = 0.0;
  solold = 0.0;

  cout << endl;

  ofstream fileout("timesumcverrns.out");
  
  double currtime = 0.0;
  for (int n=1; n<=numtimesteps; n++) 
    { 
      DualMesh *dualmesh = new DualMesh(mesh);
      cout << "//********************************************" << endl;
      cout << "Current Time: " << endl << currtime/(24.0*60.0*60.0) << " Days, " << 100*currtime/finaltime << "%" << endl; 

      currtime += dt;
      //---------------------------------------------------------
      //  Solve Hydrodynamics and Poroelastic System
      //---------------------------------------------------------
      SparseMatrix *globalmat;
      globalmat = new SparseMatrix(3*nv, 50);
      Vector globalrhs(3*nv);
      globalrhs = 0.0;

      cout << endl << "Building Global System:" << endl;
      cout << "Elasticity Block..." << flush;
      //Fill Ae to Global Matrix
      Array<double *> elfdata(2);
      elfdata[0] = new double[nv];
      elfdata[1] = new double[nv];

      for(int k=0; k<nv; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  double coords[3];
	  coords[0] = coord[0];
	  coords[1] = coord[1];
	  coords[2] = currtime;
	  elfdata[0][k] = forcex->Eval(coords);
	  elfdata[1][k] = forcey->Eval(coords);
	}
      elasticdata->SetForceData(elfdata);

      for(int b=0; b<mesh->GetNBdrs(); b++)
	{
	  int nvb = mesh->GetNBdryElem(b)+1;
	  for(int i=0; i<2; i++)
	    {
	      int bb = b + i*mesh->GetNBdrs();
	      for(int k=0; k<nvb; k++)
		{
		  int indie = mesh->GetBdryGindex(b, k);
		  double* coord = mesh->GetVertex(indie);
		  double coords[3];
		  coords[0] = coord[0];
		  coords[1] = coord[1];
		  coords[2] = currtime;
		  if(elBdrDirichlet[bb]==true){ eldbdryval[bb][k] = elasticboundaryfuncs[bb]->Eval(coords); }
		  else if(elBdrNeumann[bb]==true){ elnbdryval[bb][k] = elasticboundaryfuncs[bb]->Eval(coords); }
		}
	    }
	}
      elasticdata->SetDirichletData(eldbdryval, elBdrDirichlet);
      elasticdata->SetNeumannData(elnbdryval, elBdrNeumann);

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

      Vector urhs(2*nv);
      elasticfem->getRHS(urhs);      
      for(int j=0; j<2*nv; j++)
	{
	  globalrhs(j) += urhs(j);
	}
      cout << "done" << endl;

      //Fill Ap to Global Matrix  
      cout << "Pressure Block..." << flush;
      Array<double *> forcedata(1);
      forcedata[0] = new double[nv];

      for(int k=0; k<nv; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  double coords[3];
	  coords[0] = coord[0];
	  coords[1] = coord[1];
	  coords[2] = currtime;
	  forcedata[0][k] = force->Eval(coords);
	}
      pressuredata->SetForceData(forcedata);

      for(int b=0; b<nbdry; b++)
	{
	  int nvb = mesh->GetNBdryElem(b)+1;
	  for(int k=0; k<nvb; k++)
	    {
	      int indie = mesh->GetBdryGindex(b, k);
	      double* coord = mesh->GetVertex(indie);
	      double coords[3];
	      coords[0] = coord[0];
	      coords[1] = coord[1];
	      coords[2] = currtime;
	      if(BdrDirichlet[b]==true){ dbdryval[b][k] = pressureboundaryfuncs[b]->Eval(coords); }
	      else if(BdrNeumann[b]==true){ nbdryval[b][k] = pressureboundaryfuncs[b]->Eval(coords); }
	    }
	}     
      pressuredata->SetDirichletData(dbdryval, BdrDirichlet);
      pressuredata->SetNeumannData(nbdryval, BdrNeumann);
      
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

      Vector prhs(nv);
      pressurefem->getRHS(prhs);      
      for(int j=0; j<nv; j++)
	{
	  globalrhs(j+2*nv) += prhs(j);
	}
      
      cout << "done" << endl;
      cout << "Stabilization Block..." << flush;
      Array<double *> corrector(1);
      corrector[0] = new double[nv];
      for(int i=0; i<nv; i++)
	{
	  corrector[0][i] = 0.0; //0.25*h*h/(double(dt));
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
	      globalmat->Elem(k+2*nv, J[i]+2*nv) += betamatrix(k, J[i]);
	    }
	}

      //Add Corrector to RHS 
      Vector betap0(nv);
      betamatrix.Mult(pressureold, betap0);
      for(int k=0; k<nv; k++)
	{
	  globalrhs(k+2*nv) += betap0(k);
	}

      cout << "done" << endl;
      cout << "Off Diagonal Blocks..." << flush;      
      //Set off diagonal data
      Data *scalardatax;
      scalardatax = new Data();
      
      Array<double *> scalarcoeffx(2);
      scalarcoeffx[0] = new double[nv];
      scalarcoeffx[1] = new double[nv];
      for(int i=0; i<nv; i++)
	{
	  scalarcoeffx[0][i] = 1.0; 
	  scalarcoeffx[1][i] = 0.0; 
	}
      scalardatax->SetAdvectionData(scalarcoeffx);
      
      FEM *scalarfemx;
      if(meshtype==Element::QUADRILATERAL){ scalarfemx = new BilinearFEM2D(mesh, scalardatax); }
      if(meshtype==Element::TRIANGLE){ scalarfemx = new LinearFEM2D(mesh, scalardatax); }
      scalarfemx->AssembleMatrix();
      scalarfemx->FinalizeMatrix(); 
      SparseMatrix Q = scalarfemx->getMatrix();

      I = Q.GetI();
      J = Q.GetJ();
      for(int k=0; k<nv; k++)
	{
	  for(int i=I[k]; i<I[k+1]; i++)
	    {
	      globalmat->Elem(k, J[i]+2*nv) = Q(k, J[i]);
	      globalmat->Elem(k+2*nv, J[i]) = Q(k, J[i])/dt;
	    }
	}
      Data *scalardatay;
      scalardatay = new Data();
      
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
      SparseMatrix R = scalarfemy->getMatrix();

      I = R.GetI();
      J = R.GetJ();
      for(int k=0; k<nv; k++)
	{
	  for(int i=I[k]; i<I[k+1]; i++)
	    {
	      globalmat->Elem(k+nv, J[i]+2*nv) = R(k, J[i]);
	      globalmat->Elem(k+2*nv, J[i]+nv) = R(k, J[i])/dt;
	    }
	}


      Vector Qu1(nv);
      Vector Ru2(nv);
      Q.Mult(u1old, Qu1);
      R.Mult(u2old, Ru2);
      
      for(int j=0; j<nv; j++)
	{
	  globalrhs(j+2*nv) += (Qu1(j)+Ru2(j))/dt;
	}

      cout << "done" << endl;
      cout << "Boundary Conditions..." << flush;  
  
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
      cout << "done" << endl << endl;
      cout << "Solving..." << flush;  

      globalmat->Finalize();
      int mem_fact = 100;
      MultiFrontalMatrixInverse *InvMat = new MultiFrontalMatrixInverse(*globalmat, mem_fact);
      InvMat->Mult(globalrhs, sol);
      delete InvMat;
      cout << "done" << endl;

      for(int i=0; i<nv; i++){u1(i) = sol(i); u2(i) = sol(i+nv); pressure(i) = sol(i+2*nv);}
      for(int i=0; i<2*nv; i++){displacement(i) = sol(i); }
      
      //---------------------------------------------------------
      //  Post-Proccessing
      //---------------------------------------------------------
      Array<double *> flux;
      Array<double *> ppsol;
      GetPPDarcy(dualmesh, pressuredata, sol, solold, pressurefem, scalarfemx, scalarfemy, betafem, dt, flux, ppsol);

      double pherr = ComputeDiscH1SemiNormError(mesh, solpx, solpy, ppsol, currtime);
      cout << "The H1 SemiNorm error for the post-processed pressure is: " << pherr << endl;

      double preherr = ComputeH1SemiNormError(mesh, solpx, solpy, pressure, currtime);
      cout << "The H1 SemiNorm error for the pressure is: " << preherr << endl;

      // Post-processing for displacements
      Vector ut(2*nv);
      for(int i=0; i<2*nv; i++) ut(i) = ( displacement(i) - displacementold(i) ) / dt;

      Array<double *> uhdotn(dualmesh->GetNumLocalDOF());
      for(int i=0; i<uhdotn.Size(); i++ ) uhdotn[i] = new double[dualmesh->GetNE()];

      Array<double> buhdotn(nv);
      Computeuhdotn(dualmesh, ut, uhdotn, buhdotn);

      Array<double> bdrflux(nv);
      for(int i=0; i<bdrflux.Size(); i++){ bdrflux[i] = 0.0; }
      for(int b=0; b<mesh->GetNBdrs(); b++)
	{
	  for(int i=0; i<=mesh->GetNBdryElem(b); i++)
	    {
	      int numlocaldof = dualmesh->GetNumLocalDOF();
	      int vindex = mesh->GetBdryGindex(b, i);
	      double internalut = 0.0;
	      for(int k=0; k<dualmesh->Getcvnumel(vindex); k++)
		{
		  int el = dualmesh->Getcvelindex(vindex, k);
		  int locvert = dualmesh->Getcvlocalindex(vindex,k);
		  internalut += uhdotn[locvert][el];

		  int locbdryind = (locvert+numlocaldof-1)%numlocaldof; 
		  internalut -= uhdotn[locbdryind][el];
		}
	       bdrflux[vindex] = -internalut; 
	    }
	}    
  
      Vector force(nv);
      for(int k=0; k<nv; k++){ force(k) = pressuredata->GetNodalForceCoeff(k); }
      Array<double> tempbdrflux;
      dualmesh->GetFluxonBdrCVEdges(flux, tempbdrflux, force);
      for(int i=0; i<nv; i++){ bdrflux[i] += tempbdrflux[i]-buhdotn[i]; }
      
      //Data for upwinding
      Array<FData*> fluxdata(2);
      fluxdata[0] = new FData(flux, bdrflux);
      fluxdata[1] = new FData(uhdotn, buhdotn);
      
      //---------------------------------------------------------
      //  Check Local Conservation
      //---------------------------------------------------------
      // cout << "Local Conservation Checking..." << flush;
      Vector cverr(nv);
      cverr = 0.0;
      for(int k=0; k<mesh->GetNE(); k++)
	{
	  Array<int> ind;
	  dualmesh->GetDOF(k, ind);

	  Array<double> areas;     
	  dualmesh->GetAreas(k, areas);
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

	  for(int j=0; j<ind.Size(); j++){ cverr(ind[j]) -= pressuredata->GetNodalForceCoeff(ind[j])*areas[j]; }
	}
    
      double maxcverror= 0.0;
      int maxindex;
      for(int k=0; k<nv; k++)
	{
	  cverr(k) +=  bdrflux[k] + buhdotn[k]; 
	  // cout << k << " " << cverr(k) << endl;
	  if(maxcverror<fabs(cverr(k))){ maxcverror = fabs(cverr(k)); maxindex = k;}
	} 
      cout << "Local conservation error: " << maxcverror << " " << "in cv " << maxindex << endl;

      double sumcverr = 0.0;
      for(int k=0; k<nv; k++)
	{
	  sumcverr += fabs(cverr(k));
	}
      cout << "Sum of Local Conservation Error: " << sumcverr << endl;

      fileout << dt*n << "\t" << sumcverr << endl;

      //ofstream fileout("lce.out");
      //for(int k=0; k<nv; k++)
      //fileout << k << "\t" << cverr(k) << endl;
      //fileout.close(); 

      //sprintf(filename, "localconservationerror_%f.vtk", currtime);
      //output.open(filename);
      //mesh->ParaviewPrint(output, cverr); 
      //output.close();    

      for(int i=0; i<nv; i++){ ; }
  
      Vector u1t(nv);
      Vector u2t(nv);
      for(int k=0; k<nv; k++)
	{
	  u1t(k) = (u1(k)-u1old(k))/dt;
	  u2t(k) = (u2(k)-u2old(k))/dt;
	} 
      //Estimate errors
      cout << "L2error pressure: " << pressurefem->ComputeL2Error(solp, pressure, currtime) << endl;	  
      cout << "L2error u1: " << elasticfem->ComputeL2Error(solu1, u1, currtime) << endl;
      cout << "L2error u2: " << elasticfem->ComputeL2Error(solu2, u2, currtime) << endl;
      cout << "L2error u1t: " << pressurefem->ComputeL2Error(solu1t, u1t, currtime) << endl;
      cout << "L2error u2t: " << elasticfem->ComputeL2Error(solu2t, u2t, currtime) << endl;

   
      pressureold = pressure;
      displacementold = displacement;
      u1old = u1;
      u2old = u2;
      solold = sol;

      //---------------------------------------------------------
      //  Print
      //---------------------------------------------------------
      cout << "Printing..." << flush;  
      sprintf(filename, "pressure_%f.vtk", currtime);
      output.open(filename);
      mesh->ParaviewPrint(output, pressureold); 
      output.close();  

      sprintf(filename, "uone_%f.vtk", currtime);
      output.open(filename);
      mesh->ParaviewPrint(output, u1old); 
      output.close();  

      sprintf(filename, "utwo_%f.vtk", currtime);
      output.open(filename);
      mesh->ParaviewPrint(output, u2old); 
      output.close();  
      
      sprintf(filename, "defmesh_%f.vtk", currtime);
      output.open(filename);
      mesh->ParaviewPrintDeformedMesh(output, displacementold);
      output.close();
      
      for(int i=0; i<uhdotn.Size(); i++){ delete []uhdotn[i]; } 
      for(int i=0; i<fluxdata.Size(); i++){ delete fluxdata[i]; }
      for(int i=0; i<ppsol.Size(); i++){ delete []ppsol[i]; }
      for(int i=0; i<flux.Size(); i++){ delete []flux[i]; }
      for(int k=0; k<elfdata.Size(); k++){ delete []elfdata[k]; } 
      for(int k=0; k<scalarcoeffx.Size(); k++){ delete []scalarcoeffx[k]; }
      for(int k=0; k<scalarcoeffy.Size(); k++){ delete []scalarcoeffy[k]; };
      delete []forcedata[0];	  
      delete []corrector[0];
      delete globalmat;
      delete betafem;
      delete betadata;
      delete scalardatax;
      delete scalardatay;
      delete pressurefem;
      delete elasticfem;
      delete scalarfemx;
      delete scalarfemy;
      delete dualmesh;
      cout << endl;
    }

  fileout.close();

  //---------------------------------------------------------
  //  Freedom!!!!
  //---------------------------------------------------------
  delete mesh;
  delete pressuredata;
  delete elasticdata;
  delete force;
  delete forcex;
  delete forcey;
  delete []permdata[0];

  delete solp;
  delete solu1;
  delete solu2;
  delete solpx;
  delete solpy;
  delete solu1t;
  delete solu2t;

  for(int k=0; k<lamecoeff.Size(); k++)
    delete []lamecoeff[k];

  for(int k=0; k<dbdryval.Size(); k++)
    delete []dbdryval[k];

  for(int k=0; k<nbdryval.Size(); k++)
    delete []nbdryval[k];

  for(int k=0; k<eldbdryval.Size(); k++)
    delete []eldbdryval[k];

  for(int k=0; k<elnbdryval.Size(); k++)
    delete []elnbdryval[k];
 
  for (int k=0; k<elasticboundaryfuncs.Size(); k++)
    delete elasticboundaryfuncs[k]; 

  for (int k=0; k<pressureboundaryfuncs.Size(); k++)
    delete pressureboundaryfuncs[k]; 


  cout << endl << "fin" << endl;
}
