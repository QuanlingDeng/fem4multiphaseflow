#include <iostream>
#include <fstream>
#include <cmath>
#include "../linalg/linalg_header.h"
#include "../mesh/mesh_header.h"
#include "../fem/fem_header.h"
#include "../timedependent/timedependent_header.h"

using namespace std;

//======================================================
void ComputePPuhdotn(DualMesh *dual, double dt, Vector &displacementold, Vector &displacement,
		     Array<double *> &flux, Array<double *> &ppuhdotn)
{
  int NE = dual->GetNE();
  int NV = dual->GetNV();
  int n = dual->GetNumLocalDOF();
  Array<double *> normals(n);
  Array<double> lengths;
  Array<double *> ecoord(n+1);  
  Vector dpo(2*n);
  Vector dp(2*n);
  
  for( int i=0; i<NE; i++)
    {
      dual->GetCoordinates(i, ecoord);
      dual->GetNormals(i, normals);
      dual->GetEdgeLengths(i, lengths);

      Array<int> ind;
      dual->GetDOF(i, ind);
      for(int k=0; k<n; k++)
	{
	  dpo(k) = displacementold( ind[k] );
	  dpo(k+n) = displacementold( ind[k] + NV );
	  dp(k) = displacement( ind[k] );
	  dp(k+n) = displacement( ind[k] + NV );
	}
	
      if( n == 3)
	{
	  double utx = normals[0][0] * lengths[0] * ( 5.0 * (dpo(0) - dp(0)) + 5.0 * (dpo(1) - dp(1)) + 2.0 * (dpo(2) - dp(2)) ) / (12.0 * dt);
	  double uty = normals[0][1] * lengths[0] * ( 5.0 * (dpo(3) - dp(3)) + 5.0 * (dpo(4) - dp(4)) + 2.0 * (dpo(5) - dp(5)) ) / (12.0 * dt);      
	  double tem = utx + uty;

	  ppuhdotn[0][i] = tem;
	  ppuhdotn[1][i] = tem + flux[1][i] - flux[0][i];
	  ppuhdotn[2][i] = tem + flux[2][i] - flux[0][i];	  
	}
      else
	{
	  double utx = normals[0][0] * lengths[0] * ( 3.0 * (dpo(0) - dp(0)) + 3.0 * (dpo(1) - dp(1)) + (dpo(2) - dp(2)) + (dpo(3) - dp(3)) ) / (8.0 * dt);
	  double uty = normals[0][1] * lengths[0] * ( 3.0 * (dpo(4) - dp(4)) + 3.0 * (dpo(5) - dp(5)) + (dpo(6) - dp(6)) + (dpo(7) - dp(7)) ) / (8.0 * dt);      
	  double tem = utx + uty;

	  ppuhdotn[0][i] = tem;
	  ppuhdotn[1][i] = tem + flux[1][i] - flux[0][i];
	  ppuhdotn[2][i] = tem + flux[2][i] - flux[0][i];	  
	  ppuhdotn[3][i] = tem + flux[3][i] - flux[0][i];	  
	}
    }

  for( int i=0; i<NE; i++ )
    for( int j=0; j<n; j++ ) ppuhdotn[j][i] = -ppuhdotn[j][i];
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
  Array<double *> forcedata(1);
  forcedata[0] = new double[nv];

  ReadIdentifier(input_file, buffer, bufflen);
  Function *force = ReadFunction(input_file);
  for(int k=0; k<nv; k++)
    {
      double* coord = mesh->GetVertex(k);
      forcedata[0][k] = force->Eval(coord);
    }
  delete force;

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
      int n = mesh->GetNBdryElem(b)+1;
      nbdryval[b] = new double[n];
      rbdryval[b] = new double[n];
      rcoeff[b] = new double[n];
      dbdryval[b] = new double[n];
      
      input_file.getline( buffer, bufflen, ' ' );
      input_file >> bdryflag;

      if (!strcmp(bdryflag, "neumann"))
	{
	  BdrNeumann[b] = true;
	  BdrRobin[b] = false;
	  BdrDirichlet[b] = false;

	  ReadIdentifier(input_file, buffer, bufflen);
	  Function *bdryval = ReadFunction(input_file);
	  for(int k=0; k<n; k++)
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
	  for(int k=0; k<n; k++)
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
	  for(int k=0; k<n; k++)
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
  pressuredata->SetForceData(forcedata);
  pressuredata->SetNeumannData(nbdryval, BdrNeumann);
  pressuredata->SetRobinData(rcoeff, rbdryval, BdrRobin);
  pressuredata->SetDirichletData(dbdryval, BdrDirichlet);
  cout << "Pressure Data Set" << endl;


  //---------------------------------------------------------
  //  Gather Data for elasticity part
  //--------------------------------------------------------- 

  int nnn = 2*nv; //number of dof

  //Gather elasticitiy Data
  Array<double *> lamecoeff(2);
  lamecoeff[0] = new double[nv];
  lamecoeff[1] = new double[nv];
  Vector elmod(nv);
  Array<double *> elforcedata(2);
  elforcedata[0] = new double[nv];
  elforcedata[1] = new double[nv];

  char elmodflag[bufflen]; 
  ReadIdentifier(input_file, buffer, bufflen);
  input_file >> elmodflag;
  
  if (!strcmp(elmodflag, "equation"))
    {
      Function *modfunc;
      ReadIdentifier(input_file, buffer, bufflen);
      modfunc = ReadFunction(input_file);
      
      double prat = 0.25;
      for(int k=0; k<nv; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  double data = modfunc->Eval(coord);
	  elmod(k) = data;
	  lamecoeff[0][k] = data*prat/((1.0+prat)*(1.0-2.0*prat));
	  lamecoeff[1][k] = data/(1+prat);

	}
      delete modfunc;
    }
  mesh->ParaviewPrint(elmod, meshtype);

  for(int i=0; i<elforcedata.Size(); i++)
    {
      ReadIdentifier(input_file, buffer, bufflen);
      Function *force = ReadFunction(input_file);
      for(int k=0; k<nv; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  elforcedata[i][k] = force->Eval(coord);
	}
      delete force;
    }

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
      int nv = mesh->GetNBdryElem(b)+1;
      for(int i=0; i<2; i++)
	{
	  int bb = b + i*nbredges;
	  elnbdryval[bb] = new double[nv];
	  elrbdryval[bb] = new double[nv];
	  elrcoeff[bb] = new double[nv];
	  eldbdryval[bb] = new double[nv];
	  
	  input_file.getline( buffer, bufflen, ' ' );
	  input_file >> bdryflag;
	  
	  if (!strcmp(bdryflag, "neumann"))
	    {
	      elBdrNeumann[bb] = true;
	      elBdrRobin[bb] = false;
	      elBdrDirichlet[bb] = false;
	      
	      ReadIdentifier(input_file, buffer, bufflen);
	      Function *bdryval = ReadFunction(input_file);
	      for(int k=0; k<nv; k++)
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
	      for(int k=0; k<nv; k++)
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
	      for(int k=0; k<nv; k++)
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
  elasticdata->SetForceData(elforcedata);
  elasticdata->SetNeumannData(elnbdryval, elBdrNeumann);
  elasticdata->SetRobinData(elrcoeff, elrbdryval, elBdrRobin);
  elasticdata->SetDirichletData(eldbdryval, elBdrDirichlet);

  cout << "Elastic Data Set" << endl;

  //---------------------------------------------------------
  //  Set initial profiles
  //---------------------------------------------------------

  //Set initial displacment
  ReadIdentifier(input_file, buffer, bufflen);
  cout << buffer << endl;
  Function *initu1 = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *initu2 = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *initp = ReadFunction(input_file);

  Vector displacementold(2*nv);
  Vector pold(nv);
  displacementold = 0.0;
  for(int k=0; k<nv; k++)
    {
      double* coord = mesh->GetVertex(k);
      displacementold(k) = initu1->Eval(coord);
      displacementold(k+nv) = initu2->Eval(coord);
      pold(k) = initp->Eval(coord);
    }
  delete initu1;
  delete initu2;

  
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
  //  Time March 
  //---------------------------------------------------------
  SparseMatrix *globalmat;
  globalmat = new SparseMatrix(3*nv, 3*nv);
  Vector globalrhs(3*nv);
  globalrhs = 0.0;

  //Fill Ae to Global Matrix
  FEM *elasticfem;
  if(meshtype==Element::QUADRILATERAL){ elasticfem = new BilinearFEMELAST2D(mesh, elasticdata); }
  if(meshtype==Element::TRIANGLE){ elasticfem = new LinearFEMELAST2D(mesh, elasticdata); }
  elasticfem->AssembleMatrix();
  elasticfem->FinalizeMatrix();
  SparseMatrix elastmatrix = elasticfem->getMatrix();
  Vector urhs;
  elasticfem->getRHS(urhs);

  int *I = elastmatrix.GetI();
  int *J = elastmatrix.GetJ();
  for(int k=0; k<2*nv; k++)
    {
       for(int i=I[k]; i<I[k+1]; i++)
	{
	  globalmat->Elem(k, J[i]) = elastmatrix(k, J[i]);
	}
    }
 
  for(int j=0; j<2*nv; j++)
    {
      globalrhs(j) += urhs(j);
    }
  
  //Fill Ap to Global Matrix
  Vector divuold;
  calculatedivu(mesh, displacementold, divuold);
  for(int i=0; i<nv; i++)
    {
      forcedata[0][i] += divuold(i)/dt;
      cout << forcedata[0][i] << endl;
    }
  pressuredata->SetForceData(forcedata);
  cout << endl;
  FEM *pressurefem;
  if(meshtype==Element::QUADRILATERAL){ pressurefem = new BilinearFEM2D(mesh, pressuredata); }
  if(meshtype==Element::TRIANGLE){ pressurefem = new LinearFEM2D(mesh, pressuredata); }
  pressurefem->AssembleMatrix(); 
  pressurefem->FinalizeMatrix();
  SparseMatrix pressmatrix = pressurefem->getMatrix();
  Vector prhs;
  pressurefem->getRHS(prhs);


  I = pressmatrix.GetI();
  J = pressmatrix.GetJ();
  for(int k=0; k<nv; k++)
    {
      for(int i=I[k]; i<I[k+1]; i++)
	{
	  globalmat->Elem(k+2*nv, J[i]+2*nv) = pressmatrix(k, J[i]);
	}
    }

  for(int j=0; j<nv; j++)
    {
      globalrhs(j+2*nv) += prhs(j);
      cout << prhs(j) << endl;
    }
  
  cout << endl;
  Array<double *> corrector(1);
  corrector[0] = new double[nv];
  for(int i=0; i<nv; i++)
    {
      corrector[0][i] = h*h/(4*dt*(lamecoeff[0][i]+lamecoeff[1][i]));
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

  Vector betap0(nv);
  betamatrix.Mult(pold, betap0);
  for(int k=0; k<nv; k++)
    {
      globalrhs(k+2*nv) += betap0(k);
      cout << betap0(k) << endl;
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
	  globalmat->Elem(J[i]+2*nv, k) = -scalarxmatrix(k, J[i])/dt;
	}
    }
  
  /*
  Vector u1old(nv);
  Vector Atu1(nv);
  for(int i=0; i<nv; i++){u1old(i) = displacementold(i);}
  scalarxmatrix.MultTranspose(u1old, Atu1);
  for(int k=0; k<nv; k++)
    {
      globalrhs(k+2*nv) += -Atu1(k)/dt;
    }
  */
  delete scalarfemx;


  for(int i=0; i<nv; i++)
    {
      scalarcoeff[0][i] = 0.0; 
      scalarcoeff[1][i] = 1.0; 
    }
  scalardata->SetAdvectionData(scalarcoeff);

  FEM *scalarfemy;
  if(meshtype==Element::QUADRILATERAL){ scalarfemy = new BilinearFEM2D(mesh, scalardata); }
  if(meshtype==Element::TRIANGLE){ scalarfemy = new LinearFEM2D(mesh, scalardata); }
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
	  globalmat->Elem(J[i]+2*nv, k+nv) = -scalarymatrix(k, J[i])/dt;
	}
    }
  
  /*
  Vector u2old(nv);
  Vector Atu2(nv);
  for(int i=0; i<nv; i++){u2old(i) = displacementold(i+nv);}
  scalarxmatrix.MultTranspose(u2old, Atu2);
  for(int k=0; k<nv; k++)
    {
      globalrhs(k+2*nv) += -Atu2(k)/dt;
    }
  */
  delete scalarfemy;


  /*
  for(int i=0; i<nv; i++)
    {
      scalarcoeff[0][i] = 0.0; 
      scalarcoeff[1][i] = 0.0; 
    }
  scalardata->SetAdvectionData(scalarcoeff);

  Array<double *> beta(1);
  beta[0] = new double[nv];

  for(int k=0; k<nv; k++)
    {
      beta[0][k] = 0.0/dt; //2.0/(dt*(2.0*lamecoeff[0][k]+lamecoeff[1][k]));
    }   
  scalardata->SetReactionData(beta);

  FEM *scalarfemreaction;
  if(meshtype==Element::QUADRILATERAL){ scalarfemreaction = new BilinearFEM2D(mesh, scalardata); }
  if(meshtype==Element::TRIANGLE){ scalarfemreaction = new LinearFEM2D(mesh, scalardata); }
  scalarfemreaction->AssembleMatrix(); 
  scalarfemreaction->FinalizeMatrix(); 
  SparseMatrix reactionmatrix = scalarfemreaction->getMatrix();

  I = scalarymatrix.GetI();
  J = scalarymatrix.GetJ();
  for(int k=0; k<nv; k++)
    {
       for(int i=I[k]; i<I[k+1]; i++)
	{
	  globalmat->Elem(k+2*nv, J[i]+2*nv) += reactionmatrix(k, J[i]);
	}
    }
*/
  /*
  double error = 0.0;
  for(int j=0; j<nv; j++)
    {
      error += fabs( -(Atu1(j)+Atu2(j))/dt - ukmterms(j));
      cout << j << "\t" << -(Atu1(j)+Atu2(j))/dt << "\t \t" << ukmterms(j) << endl;
    }
  cout << error << endl;
  cout << endl;
  */
  
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
  globalmat->Print();
  cout << endl;

  globalrhs.Print();
  Vector pressure(nv);
  Vector sol(3*nv);
  pressure = 1.0;
  sol = 0.0;

  
  int mem_fact = 100;
  MultiFrontalMatrixInverse *InvMat = new MultiFrontalMatrixInverse(*globalmat, mem_fact);
  InvMat->Mult(globalrhs, sol);
  delete InvMat;
  
  //update pressure
  for(int i=0; i<nv; i++){ pressure(i) = sol(i+2*nv); }



  /*
  //Iteration
  int MaxNumIter = 1000; double tol = 1.0e-6; double maxerr = 0.0; double maxerrprev=10.0;
  for (int k=0; k<MaxNumIter; k++)
    {
      maxerr = 0.0;
      Vector Ap(nv);
      Ap = 0.0;
      reactionmatrix.Mult(pressure, Ap);
      
      for(int j=0; j<nv; j++)
	{
	  globalrhs(j+2*nv) = ukmterms(j) + Ap(j);
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
		      int index = mesh->GetBdryGindex(b, i);
		      double dval = pressuredata->GetDirichletBdryVal(b, i);
		      globalrhs(index+2*nv) = dval;
		    }	   
		}
	    }
	}

      //Solve
      sol = 0.0;
      int mem_fact = 100;
      MultiFrontalMatrixInverse *InvMat = new MultiFrontalMatrixInverse(*globalmat, mem_fact);
      InvMat->Mult(globalrhs, sol);
      delete InvMat;
      
      
      //compute max-norm error the difference between iterations
      for(int i=0; i<nv; i++)
	{
	  double temperror = fabs(sol(i+2*nv) - pressure(i));
	  if( maxerr < temperror ){  maxerr = temperror; }
	}
      cout << maxerr << endl;

      //update pressure
      for(int i=0; i<nv; i++){ pressure(i) = sol(i+2*nv); }
      
      if( maxerr > maxerrprev )
	{
	  cout << "Divergent!!!" << endl;
	  break;
	}
      maxerrprev = maxerr;
      
      //check if convergence is reached
      if( maxerr < tol )
	{
	  cout << "done." << endl << "Convergence reached at iteration: " << k << " with maxerr: " << maxerr << endl;
	  break;
	}
      if(k==MaxNumIter-1)
	{
	  cout << "Convergence not reached. " << endl;
	  break;
	} 
      
    }
  */
  //Estimate errors
  Vector u1(nv);
  Vector u2(nv);
  Vector displacement(2*nv);

  u1 = 0.0;
  u2 = 0.0;
  for(int i=0; i<nv; i++){u1(i) = sol(i); u2(i) = sol(i+nv); pressure(i) = sol(i+2*nv);}
  for(int i=0; i<2*nv; i++){displacement(i) = sol(i);}

  //for(int i=0; i<3*nv; i++){cout << i%nv << "\t" << sol(i) << endl; }
  char filename[256];
  sprintf(filename, "pressure.vtk");
  ofstream output;
  output.open(filename);
  mesh->ParaviewPrint(output, pressure, meshtype); 
  output.close();  
  
  sprintf(filename, "u1.vtk");
  output.open(filename);
  mesh->ParaviewPrint(output, u1, meshtype); 
  output.close(); 

  sprintf(filename, "u2.vtk");
  output.open(filename);
  mesh->ParaviewPrint(output, u2, meshtype); 
  output.close();  
  
  sprintf(filename, "defmesh.vtk");
  output.open(filename);
  mesh->ParaviewPrintDeformedMesh(output, displacement, meshtype);
  output.close();  
     
  cout << "L2error pressure: " << pressurefem->ComputeL2Error(solp, pressure) << endl;	 
  cout << "L2error u1: " << elasticfem->ComputeL2Error(solu1, u1) << endl;
  cout << "L2error u2: " << elasticfem->ComputeL2Error(solu2, u2) << endl;
	  
  //---------------------------------------------------------
  //  Freedom!!!!
  //---------------------------------------------------------

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

  delete []permdata[0];
  delete []forcedata[0];

 
  for(int k=0; k<scalarcoeff.Size(); k++)
    delete []scalarcoeff[k];

  for(int k=0; k<lamecoeff.Size(); k++)
    delete []lamecoeff[k];

  for(int k=0; k<elforcedata.Size(); k++)
    delete []elforcedata[k];

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

