
#include <cmath>
#include "../fem/fem_header.h"
#include "../timedependent/timedependent_header.h"

using namespace std;

double L2Error(int nx, int ny, Vector &sol, Vector &rsol)
{
  double err = 0.0;
  double gp[8] = { 1.0/3.0, 1.0/3.0, 0.2, 0.6, 0.2, 0.2, 0.6, 0.2 };
  double gw[4] = { -27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0 };
  Array<double> x(3);
  Array<double> y(3);
  double dx = 1.0/nx;
  double dy = 1.0/ny;
   
  for(int i=0; i<nx; i++)
    {
      for(int j=0; j<ny; j++)
	{
	  //---low triangle----------
	  Array<int> ind(3);
	  ind[0] = j*(nx+1) + i;
	  ind[1] = ind[0] + 1;
	  ind[2] = ind[0] + nx + 1;

	  x[0] = i*dx; 
	  y[0] = j*dy;
	  x[1] = (i+1)*dx; 
	  y[1] = j*dy;
	  x[2] = i*dx; 
	  y[2] = (j+1)*dy;
      
	  double t1 = x[1] - x[0];
	  double t2 = x[2] - x[0];
	  double t3 = y[1] - y[0];
	  double t4 = y[2] - y[0];
	  double detJ = t1 * t4 - t2 * t3; 

	  for(int k=0; k<4; k++)
	    {
	      double xx = gp[2*k];
	      double yy = gp[2*k+1];

	      double uh = sol(ind[0]) * (1 - xx - yy) + sol(ind[1]) * xx + sol(ind[2]) * yy;
	      double ruh = rsol(ind[0]) * (1 - xx - yy) + rsol(ind[1]) * xx + rsol(ind[2]) * yy;

	      double tem = uh - ruh;
	      err += gw[k] * tem * tem * detJ;	 
	    }

	  //---upper triangle-----------
	  ind[0] = ind[1];
	  ind[1] = ind[0] + nx + 1;
	  ind[2] = ind[1] - 1;

	  x[0] = (i+1)*dx; 
	  y[0] = j*dy;
	  x[1] = (i+1)*dx; 
	  y[1] = (j+1)*dy;
	  x[2] = i*dx; 
	  y[2] = (j+1)*dy;
      
	  t1 = x[1] - x[0];
	  t2 = x[2] - x[0];
	  t3 = y[1] - y[0];
	  t4 = y[2] - y[0];
	  detJ = t1 * t4 - t2 * t3; 

	  for(int k=0; k<4; k++)
	    {
	      double xx = gp[2*k];
	      double yy = gp[2*k+1];

	      double uh = sol(ind[0]) * (1 - xx - yy) + sol(ind[1]) * xx + sol(ind[2]) * yy;
	      double ruh = rsol(ind[0]) * (1 - xx - yy) + rsol(ind[1]) * xx + rsol(ind[2]) * yy;

	      double tem = uh - ruh;
	      err += gw[k] * tem * tem * detJ;	 
	    }
	}
    }

  return sqrt(err);
}

//======================================================

int main(int argc, char *argv[])
{
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

  Array<int> aN(2);
  aN[0] = 2*N[0];
  aN[1] = 2*N[1];
  
  Mesh *mesh = new TMesh<2>(N, L, Element::TRIANGLE);
  Mesh *amesh = new TMesh<2>(aN, L, Element::TRIANGLE);

  cout << endl << "Mesh Size: " << N[0] << "\t" << N[1] << endl;

  int nv = mesh->GetNV();
  //int nbe = mesh->GetNBE();
  int dof = nv + mesh->GetNEdges();
   
  //---------------------------------------------------------
  // Gather Data for Diffusion Equation
  //---------------------------------------------------------
  //Set Permeability Data
  Array<double *> permdata(1);
  permdata[0] = new double[dof];

  char permflag[bufflen]; 
  ReadIdentifier(input_file, buffer, bufflen);
  input_file >> permflag;
  cout << permflag << endl;
  Function *permfunc;
  ReadIdentifier(input_file, buffer, bufflen);
  permfunc = ReadFunction(input_file);

  //Gather Force Data
  ReadIdentifier(input_file, buffer, bufflen);
  Function *force = ReadFunction(input_file);
  Array<double *> forcedata(1);
  forcedata[0] = new double[dof];

  Vector truepressure(dof);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *tpressure = ReadFunction(input_file);
  
  Array<int> ind;
  for(int i=0; i<mesh->GetNE(); i++)
    {
      mesh->GetElementVertices(i, ind);
      
      Array<int> edges;
      mesh->GetElementEdges(i, edges);
      for (int k = 0; k < edges.Size(); k++)
	ind.Append (mesh->GetNV() + edges[k]);

      Array<double> x(6);
      Array<double> y(6);
      for(int k=0; k<3; k++)
	{
	  x[k] = mesh->GetVertex(ind[k], 0);
	  y[k] = mesh->GetVertex(ind[k], 1);	  
	}

      x[3] = 0.5 * ( x[0] + x[1] );
      y[3] = 0.5 * ( y[0] + y[1] );
      x[4] = 0.5 * ( x[1] + x[2] );
      y[4] = 0.5 * ( y[1] + y[2] );
      x[5] = 0.5 * ( x[2] + x[0] );
      y[5] = 0.5 * ( y[2] + y[0] );
      
      for(int k=0; k<6; k++)
	{
	  double coord[2];
	  coord[0] = x[k];
	  coord[1] = y[k];
	  permdata[0][ind[k]] = permfunc->Eval(coord);
	  forcedata[0][ind[k]] = force->Eval(coord);	
	  truepressure(ind[k]) = tpressure->Eval(coord);  
	}
    }

  //Gather Boundary Data
  int nbdry = mesh->GetNBdrs();
  Array<double *> nbdryval(nbdry);
  Array<double *> dbdryval(nbdry);
  Array<bool> BdrNeumann(nbdry);
  Array<bool> BdrDirichlet(nbdry);
  Array<Function *> boundaryfuncs(nbdry);

  char bdryflag[bufflen]; 
  for(int b=0; b<nbdry; b++)
    {
      int nvb = mesh->GetNBdryElem(b)+1;
      nbdryval[b] = new double[2*nvb-1];
      dbdryval[b] = new double[2*nvb-1];
      
      input_file.getline( buffer, bufflen, ' ' );
      input_file >> bdryflag;
      ReadIdentifier(input_file, buffer, bufflen);
      boundaryfuncs[b] = ReadFunction(input_file);

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
	  int indie = mesh->GetBdryGindex(b, k);
	  double* coord = mesh->GetVertex(indie);
	  if(BdrDirichlet[b]==true){ dbdryval[b][k] = boundaryfuncs[b]->Eval(coord); }
	  else if(BdrNeumann[b]==true){ nbdryval[b][k] = boundaryfuncs[b]->Eval(coord); }
	}
      
      for(int k=0; k<nvb-1; k++)
	{
	  int indl = mesh->GetBdryGindex(b, k);
	  int indr = mesh->GetBdryGindex(b, k+1);
	  double* coordl = mesh->GetVertex(indl);
	  double* coordr = mesh->GetVertex(indr);

	  double coord[2];
	  coord[0] = 0.5 * ( coordl[0] + coordr[0] );
	  coord[1] = 0.5 * ( coordl[1] + coordr[1] );
	  	  
	  if(BdrDirichlet[b]==true){ dbdryval[b][k+nvb] = boundaryfuncs[b]->Eval(coord); }
	  else if(BdrNeumann[b]==true){ nbdryval[b][k+nvb] = boundaryfuncs[b]->Eval(coord); }
	} 
    }

  //set data
  Data *data;
  data = new Data();
  data->SetEllipticData(permdata);
  data->SetEllipticFunction(permfunc);
  data->SetForceData(forcedata);
  data->SetForceFunction(force);
  data->SetDirichletData(dbdryval, BdrDirichlet);
  data->SetNeumannData(nbdryval, BdrNeumann);
  cout << "Data Set" << endl;

  double final_time;
  int num_timesteps;
  int num_finetimesteps;
  input_file.getline( buffer, bufflen, ' ' );
  input_file >> final_time >> num_timesteps >> num_finetimesteps;

  ReadIdentifier(input_file, buffer, bufflen);
  Function *totalmobility = ReadFunction(input_file);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *fluxfuncs = ReadFunction(input_file);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *truepressurex = ReadFunction(input_file);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *truepressurey = ReadFunction(input_file);

  Array<Function *> dersol(2);
  dersol[0] = truepressurex;
  dersol[1] = truepressurey;
  
  ReadIdentifier(input_file, buffer, bufflen);
  Function *permy = ReadFunction(input_file);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *timesat = ReadFunction(input_file);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *finalsat = ReadFunction(input_file);

  // Close the main input file.
  input_file.close();

  double dt = final_time/(double)num_timesteps;
  double fdt = dt/(double)num_finetimesteps;
  Vector satold(dof);
  Vector sat(dof);
  Vector sol(dof);
  Vector tsat(dof);

  Array<double *> ttsat(3);
  for(int i=0; i<3; i++) ttsat[i] = new double[amesh->GetNE()];
  
  for(int k=0; k<dof; k++)
    {
      double* coord = amesh->GetVertex(k);
      
      double tcoord[3];
      tcoord[0] = coord[0];
      tcoord[1] = coord[1];
      tcoord[2] = 0.05;
      tsat(k) = timesat->Eval(tcoord);

      tcoord[2] = 0.0;
      satold(k) = timesat->Eval(tcoord);

      if( (coord[0] - permy->Eval(coord)*0.05) < 0.0 )
	tsat(k) = 1.0;
    }
  
  for(int i=0; i<amesh->GetNE(); i++)
    {
      amesh->GetElementVertices(i, ind);
      for(int j=0; j<3; j++)
	ttsat[j][i] = tsat(ind[j]);
    }
  

  //satold = 0.0;
  //for(int j=0; j<=aN[1]; j++){ satold(j*(aN[0]+1)) = 1.0; }
  sat = satold;
  sol = sat;

  ofstream out; ofstream output;
  char filename[256];
  sprintf(filename, "truesatat0.05.vtk");
  out.open(filename);
  amesh->ParaviewPrint(out, tsat); 
  out.close();

  //--- quadraticfem -------------
  FEM *fem = new QuadraticFEM2D(mesh, data);
  fem->Assemble();
  fem->FinalizeMatrix();
  fem->Solve(sol, "multifrontal");
  cout << "FEM Problem Solved" << endl;  

  //---------------------------------------------------------
  //  Post-Process
  //---------------------------------------------------------
  DualMesh *adual = new DualMesh(amesh); 
  DualMesh *dual = new DualMesh(mesh, 2); 

  Array<double *> flux(9);
  Array<double *> ppsol(6);
  Array<double *> ppflux(9);
  for(int i=0; i<flux.Size(); i++) flux[i] = new double[mesh->GetNE()];
  for(int i=0; i<ppsol.Size(); i++) ppsol[i] = new double[mesh->GetNE()];
  for(int i=0; i<ppflux.Size(); i++) ppflux[i] = new double[mesh->GetNE()];

  Array<double> ppbdrflux(dof);
  
  cout << "Computing Conservative Flux..." << endl;
  Postprocessing *fpp = new Postprocessing(dual, fem, data);

  fpp->ComputeConservativeFlux(truepressure, sol, flux, ppsol, ppflux, ppbdrflux);

  Array<double *> appflux(3);
  for(int i=0; i<appflux.Size(); i++) appflux[i] = new double[amesh->GetNE()];
  Array<double> appbdrflux(dof);
  for(int j=0; j<N[1]; j++)
    {
      for(int i=0; i<N[0]; i++)
	{
	  int rind = j*2*N[0] + 2*i;
	  int aind = j*8*N[0] + 4*i;

	  // lower triangle
	  appflux[0][aind] = ppflux[0][rind];
	  appflux[1][aind] = 0.5*ppflux[1][rind];
	  appflux[2][aind] = ppflux[2][rind];

	  appflux[0][aind+1] = -0.5*ppflux[4][rind];
	  appflux[1][aind+1] = -0.5*ppflux[7][rind];
	  appflux[2][aind+1] = -0.5*ppflux[1][rind];

	  appflux[0][aind+2] = ppflux[5][rind];
	  appflux[1][aind+2] = ppflux[3][rind];
	  appflux[2][aind+2] = 0.5*ppflux[4][rind];

	  appflux[0][aind+4*N[0]] = 0.5*ppflux[7][rind];
	  appflux[1][aind+4*N[0]] = ppflux[8][rind];
	  appflux[2][aind+4*N[0]] = ppflux[6][rind];

	  //upper triangle
	  appflux[0][aind+3] = ppflux[0][rind+1];
	  appflux[1][aind+3] = 0.5*ppflux[1][rind+1];
	  appflux[2][aind+3] = ppflux[2][rind+1];

	  appflux[0][aind+4*N[0]+2] = -0.5*ppflux[1][rind+1];
	  appflux[1][aind+4*N[0]+2] = -0.5*ppflux[4][rind+1];
	  appflux[2][aind+4*N[0]+2] = -0.5*ppflux[7][rind+1];

	  appflux[0][aind+4*N[0]+3] = ppflux[5][rind+1];
	  appflux[1][aind+4*N[0]+3] = ppflux[3][rind+1];
	  appflux[2][aind+4*N[0]+3] = 0.5*ppflux[4][rind+1];

	  appflux[0][aind+4*N[0]+1] = 0.5*ppflux[7][rind+1];
	  appflux[1][aind+4*N[0]+1] = ppflux[8][rind+1];
	  appflux[2][aind+4*N[0]+1] = ppflux[6][rind+1];

	  Array<int> lind;
	  Array<int> uind;
	  fem->getElementDOF(rind, lind);
	  fem->getElementDOF(rind+1, uind);
	  aind = j*2*(2*N[0]+1) + 2*i;
	  appbdrflux[aind] = ppbdrflux[lind[0]];
	  appbdrflux[aind+1] = ppbdrflux[lind[3]];
	  appbdrflux[aind+2] = ppbdrflux[lind[1]];
	  appbdrflux[aind+aN[0]+1] = ppbdrflux[lind[5]];
	  appbdrflux[aind+2*aN[0]+2] = ppbdrflux[lind[2]];

	  appbdrflux[aind+aN[0]+2] = ppbdrflux[lind[4]];
	  appbdrflux[aind+2*aN[0]+3] = ppbdrflux[uind[4]];
	  appbdrflux[aind+aN[0]+3] = ppbdrflux[uind[3]];
	  appbdrflux[aind+2*aN[0]+4] = ppbdrflux[uind[1]];
	}

    }

  Array<double> areas(dof);
  for(int i=0; i<dof; i++) areas[i] = 0.0;
  for(int i=0; i<amesh->GetNE(); i++)
    {
      Array<double> locareas(3);  
      Array<int> ind;
      amesh->GetElementVertices(i, ind);
	  
      adual->GetElemDualAreas(i, locareas);
	  
      for(int j=0; j<3; j++) areas[ind[j]] += locareas[j];
    }
      
  sprintf(filename, "saturation_0.vtk");
  out.open(filename);
  amesh->ParaviewPrint(out, satold); 
  out.close();

  int count = 0;
  while(count < num_timesteps)
    {
      //cout << "Solving Hyperbolic Problem: Percentage: " << 100*(count+1)*dt/final_time << "%" << endl;
    
      //-------time marching---------     
      TPFTimeMarchSolve(fdt, num_finetimesteps, adual, areas, appflux, appbdrflux, fluxfuncs, satold, sat);

      count ++;
      satold = sat;
      
      if(count%100==0)
	{
	  cout << "Solving Hyperbolic Problem: Percentage: " << 100*(count)*dt/final_time << "%" << endl;
	  sprintf(filename, "qsat_%d.vtk", count/100);
	  out.open(filename);
	  amesh->ParaviewPrint(out, satold); 
	  out.close();
	}
      
    }

  /*
  sprintf(filename, "ex13refq1sat_%d", aN[0]);
  out.open(filename);
  for(int i=0; i<dof; i++)
    {
      out << sat(i) << endl; 
    }
  out.close(); 
  */
  
  int nx = 256;
  int ny = 256;
  int rx = nx/aN[0];
  int ry = ny/aN[1];

  int rdof = (nx+1)*(ny+1);
  Vector rsat(rdof);
  Vector usat(rdof);
  rsat = 0.0;
  usat = 0.0;

  char refsatfile[256];
  sprintf(refsatfile, "ex13refq1sat_256");
  input_file.open( refsatfile );
  for(int i=0; i<rdof; i++)
    {
      input_file >> rsat(i);
    }
  input_file.close();

  int rind, rind0;
  for(int j=0; j<aN[1]; j++)
    {
      for(int i=0; i<aN[0]; i++)
	{
	  Array<int> ind(4);
	  ind[0] = j*(aN[0]+1) + i;
	  ind[1] = ind[0] + 1;
	  ind[2] = ind[1] + aN[0] + 1;
	  ind[3] = ind[0] + aN[0] + 1;

	  rind0 = j*ry*(nx+1) + i*rx;
	  for(int jj=0; jj<ry+1; jj++)
	    {
	      for(int ii=0; ii<rx+1; ii++)
		{
		  rind = rind0 + jj*(nx+1) + ii;
		  usat(rind) = sat(ind[0]) * ( 1.0 - ii/double(rx) ) * ( 1.0 - jj/double(ry) ) 
		    + sat(ind[1]) * ( ii/double(rx) ) * ( 1.0 - jj/double(ry) ) 
		    + sat(ind[2]) * ( ii/double(rx) ) * ( jj/double(ry) ) 
		    + sat(ind[3]) * ( 1.0 - ii/double(rx) ) * ( jj/double(ry) );
		}
	    }
	}
    }

  double maxerr = 0.0; int nn = 0;
  for(int k=0; k<rdof; k++)
    {
      //cout << k << "\t" << usat(k) << "\t" << rsat(k) << endl;
      if(fabs(rsat(k)-usat(k))>maxerr)
	{ maxerr = fabs(rsat(k)-usat(k)); nn = k; }
    }
  cout << nn << ": max saturation err: " << maxerr << endl;
  cout << "L2 saturation err: " << L2Error(nx, ny, usat, rsat) << endl;

  N[0] = 256;
  N[1] = 256;
  Mesh *mmesh = new TMesh<2>(N, L, Element::TRIANGLE);
  sprintf(filename, "refsat.vtk");
  out.open(filename);
  mmesh->ParaviewPrint(out, rsat); 
  out.close();
  
  
  Data *adata = new Data();
  FEM *afem = new LinearFEM2D(amesh, adata);
  Array<double> l2errs(3);
  Array<double> h1errs(3);
  //afem->ComputePPL2Error(l2errs, finalsat, sat, ttsat);
  //afem->ComputePPH1Error(h1errs, dersol, sat, ttsat);

  /*
  cout << "Sat L2err: " << l2errs[2] << endl;

  double maxerr = 0.0; int nn = 0;
  for(int k=0; k<dof; k++)
    if(fabs(tsat(k)-sat(k))>maxerr)
      { maxerr = fabs(tsat(k)-sat(k)); nn = k; }
  cout << nn << ": maxerr: " << maxerr << endl;
  */

  //-------------free storage--------------------------
  for(int i=0; i<ttsat.Size(); i++) delete []ttsat[i];
  for(int i=0; i<flux.Size(); i++) delete []flux[i];
  for(int i=0; i<appflux.Size(); i++) delete []appflux[i];
  for(int i=0; i<ppflux.Size(); i++) delete []ppflux[i];
  for(int i=0; i<ppsol.Size(); i++) delete []ppsol[i];

  delete afem;
  delete fem;
  delete fpp;
  delete mesh;
  delete []permdata[0];
  delete []forcedata[0];
  delete totalmobility;
  delete fluxfuncs;

  delete adata;
  delete adual;
  delete timesat;
  delete finalsat;
  delete permy;
  delete truepressurex;
  delete truepressurey;
  delete data; 
  delete dual;
  delete tpressure;
  delete force;
  delete permfunc;
  for(int k=0; k<nbdryval.Size(); k++) delete []nbdryval[k];
  for(int k=0; k<dbdryval.Size(); k++) delete []dbdryval[k];
  cout << endl << "fin" << endl;
}



/*
#include <cmath>
#include "../fem/fem_header.h"
#include "../timedependent/timedependent_header.h"

using namespace std;

//======================================================

int main(int argc, char *argv[])
{
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

  //build
  Mesh *mesh = new TMesh<2>(N, L, Element::TRIANGLE);

  cout << endl << "Mesh Size: " << N[0] << "\t" << N[1] << endl;

  int nv = mesh->GetNV();
  //int nbe = mesh->GetNBE();
  int dof = nv + mesh->GetNEdges();
   
  //---------------------------------------------------------
  // Gather Data for Diffusion Equation
  //---------------------------------------------------------
  //Set Permeability Data
  Array<double *> permdata(1);
  permdata[0] = new double[dof];

  char permflag[bufflen]; 
  ReadIdentifier(input_file, buffer, bufflen);
  input_file >> permflag;
  cout << permflag << endl;
  Function *permfunc;
  ReadIdentifier(input_file, buffer, bufflen);
  permfunc = ReadFunction(input_file);

  //Gather Force Data
  ReadIdentifier(input_file, buffer, bufflen);
  Function *force = ReadFunction(input_file);
  Array<double *> forcedata(1);
  forcedata[0] = new double[dof];

  Vector truepressure(dof);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *tpressure = ReadFunction(input_file);
  
  Array<int> ind;
  for(int i=0; i<mesh->GetNE(); i++)
    {
      mesh->GetElementVertices(i, ind);
      
      Array<int> edges;
      mesh->GetElementEdges(i, edges);
      for (int k = 0; k < edges.Size(); k++)
	ind.Append (mesh->GetNV() + edges[k]);

      Array<double> x(6);
      Array<double> y(6);
      for(int k=0; k<3; k++)
	{
	  x[k] = mesh->GetVertex(ind[k], 0);
	  y[k] = mesh->GetVertex(ind[k], 1);	  
	}

      x[3] = 0.5 * ( x[0] + x[1] );
      y[3] = 0.5 * ( y[0] + y[1] );
      x[4] = 0.5 * ( x[1] + x[2] );
      y[4] = 0.5 * ( y[1] + y[2] );
      x[5] = 0.5 * ( x[2] + x[0] );
      y[5] = 0.5 * ( y[2] + y[0] );
      
      for(int k=0; k<6; k++)
	{
	  double coord[2];
	  coord[0] = x[k];
	  coord[1] = y[k];
	  permdata[0][ind[k]] = permfunc->Eval(coord);
	  forcedata[0][ind[k]] = force->Eval(coord);	
	  truepressure(ind[k]) = tpressure->Eval(coord);  
	}
    }

  //Gather Boundary Data
  int nbdry = mesh->GetNBdrs();
  Array<double *> nbdryval(nbdry);
  Array<double *> dbdryval(nbdry);
  Array<bool> BdrNeumann(nbdry);
  Array<bool> BdrDirichlet(nbdry);
  Array<Function *> boundaryfuncs(nbdry);

  char bdryflag[bufflen]; 
  for(int b=0; b<nbdry; b++)
    {
      int nvb = mesh->GetNBdryElem(b)+1;
      nbdryval[b] = new double[2*nvb-1];
      dbdryval[b] = new double[2*nvb-1];
      
      input_file.getline( buffer, bufflen, ' ' );
      input_file >> bdryflag;
      ReadIdentifier(input_file, buffer, bufflen);
      boundaryfuncs[b] = ReadFunction(input_file);

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
	  int indie = mesh->GetBdryGindex(b, k);
	  double* coord = mesh->GetVertex(indie);
	  if(BdrDirichlet[b]==true){ dbdryval[b][k] = boundaryfuncs[b]->Eval(coord); }
	  else if(BdrNeumann[b]==true){ nbdryval[b][k] = boundaryfuncs[b]->Eval(coord); }
	}
      
      for(int k=0; k<nvb-1; k++)
	{
	  int indl = mesh->GetBdryGindex(b, k);
	  int indr = mesh->GetBdryGindex(b, k+1);
	  double* coordl = mesh->GetVertex(indl);
	  double* coordr = mesh->GetVertex(indr);

	  double coord[2];
	  coord[0] = 0.5 * ( coordl[0] + coordr[0] );
	  coord[1] = 0.5 * ( coordl[1] + coordr[1] );
	  	  
	  if(BdrDirichlet[b]==true){ dbdryval[b][k+nvb] = boundaryfuncs[b]->Eval(coord); }
	  else if(BdrNeumann[b]==true){ nbdryval[b][k+nvb] = boundaryfuncs[b]->Eval(coord); }
	} 
    }

  //set data
  Data *data;
  data = new Data();
  data->SetEllipticData(permdata);
  data->SetEllipticFunction(permfunc);
  data->SetForceData(forcedata);
  data->SetForceFunction(force);
  data->SetDirichletData(dbdryval, BdrDirichlet);
  data->SetNeumannData(nbdryval, BdrNeumann);
  cout << "Data Set" << endl;

  double final_time;
  int num_timesteps;
  int num_finetimesteps;
  input_file.getline( buffer, bufflen, ' ' );
  input_file >> final_time >> num_timesteps >> num_finetimesteps;

  ReadIdentifier(input_file, buffer, bufflen);
  Function *totalmobility = ReadFunction(input_file);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *fluxfuncs = ReadFunction(input_file);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *truepressurex = ReadFunction(input_file);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *truepressurey = ReadFunction(input_file);

  Array<Function *> dersol(2);
  dersol[0] = truepressurex;
  dersol[1] = truepressurey;
  
  ReadIdentifier(input_file, buffer, bufflen);
  Function *permy = ReadFunction(input_file);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *timesat = ReadFunction(input_file);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *finalsat = ReadFunction(input_file);

  // Close the main input file.
  input_file.close();

  double dt = final_time/(double)num_timesteps;
  double fdt = dt/(double)num_finetimesteps;
  Vector satold(dof);
  Vector sat(dof);
  Vector sol(dof);
  Vector tsat(dof);

  Array<double *> ttsat(6);
  for(int i=0; i<6; i++) ttsat[i] = new double[mesh->GetNE()];
  for(int i=0; i<mesh->GetNE(); i++)
    {
      mesh->GetElementVertices(i, ind);
      
      Array<int> edges;
      mesh->GetElementEdges(i, edges);
      for (int k = 0; k < edges.Size(); k++)
	ind.Append (mesh->GetNV() + edges[k]);

      Array<double> x(6);
      Array<double> y(6);
      for(int k=0; k<3; k++)
	{
	  x[k] = mesh->GetVertex(ind[k], 0);
	  y[k] = mesh->GetVertex(ind[k], 1);	  
	}

      x[3] = 0.5 * ( x[0] + x[1] );
      y[3] = 0.5 * ( y[0] + y[1] );
      x[4] = 0.5 * ( x[1] + x[2] );
      y[4] = 0.5 * ( y[1] + y[2] );
      x[5] = 0.5 * ( x[2] + x[0] );
      y[5] = 0.5 * ( y[2] + y[0] );
      
      for(int k=0; k<6; k++)
	{
	  double coord[3];
	  coord[0] = x[k];
	  coord[1] = y[k];
	  coord[2] = 0.05;
	  tsat(ind[k]) = timesat->Eval(coord);

	  if( (coord[0] - permy->Eval(coord)*0.05) < 0.0 )
	    tsat(ind[k]) = 1.0;

	  ttsat[k][i] = tsat(ind[k]);

	  coord[2] = 0.0;
	  satold(ind[k]) = timesat->Eval(coord);
	}
    }

  /*
  for(int j=0; j<=N[1]; j++){ satold(j*(N[0]+1)) = 1.0; }
  for(int j=0; j<nbe/4; j++)
    {
      int temp = mesh->GetBdrElementEdgeIndex(j+3*nbe/4);
      satold(temp+nv) = 1.0;
    }
  

  sat = satold;
  sol = sat;

  ofstream out; ofstream output;
  char filename[256];
  sprintf(filename, "truesatat0.05.vtk");
  out.open(filename);
  mesh->ParaviewPrint(out, tsat); 
  out.close();

  sprintf(filename, "initialsat.vtk");
  out.open(filename);
  mesh->ParaviewPrint(out, sat); 
  out.close();

  //--- quadraticfem -------------
  FEM *fem = new QuadraticFEM2D(mesh, data);
  fem->Assemble();
  fem->FinalizeMatrix();
  cout << "FEM Built" << endl;

  //---------------------------------------------------------
  //  Solve
  //---------------------------------------------------------

  fem->Solve(sol, "multifrontal");
  cout << "FEM Problem Solved" << endl;  

  //---------------------------------------------------------
  //  Post-Process
  //---------------------------------------------------------
  DualMesh *dual = new DualMesh(mesh, 2); 

  Array<double *> flux(9);
  Array<double *> ppsol(6);
  Array<double *> ppflux(9);
  for(int i=0; i<flux.Size(); i++) flux[i] = new double[mesh->GetNE()];
  for(int i=0; i<ppsol.Size(); i++) ppsol[i] = new double[mesh->GetNE()];
  for(int i=0; i<ppflux.Size(); i++) ppflux[i] = new double[mesh->GetNE()];

  Array<double> ppbdrflux(dof);
  
  cout << "Computing Conservative Flux..." << endl;
  Postprocessing *fpp = new Postprocessing(dual, fem, data);

  fpp->ComputeConservativeFlux(truepressure, sol, flux, ppsol, ppflux, ppbdrflux);

  Array<double> areas(dof);
  for(int i=0; i<dof; i++) areas[i] = 0.0;
      
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<double> locareas(1);  
      Array<int> ind;
      fem->getElementDOF(i, ind);
	  
      dual->GetElemDualAreas(i, locareas);
	  
      for(int j=0; j<3; j++) areas[ind[j]] += locareas[0];
      for(int j=3; j<6; j++) areas[ind[j]] += 3.0*locareas[0];
    }

  sprintf(filename, "saturation_0.vtk");
  out.open(filename);
  mesh->ParaviewPrint(out, satold); 
  out.close();

  int count = 0;
  while(count < num_timesteps)
    {
      cout << "Solving Hyperbolic Problem: Percentage: " << 100*(count+1)*dt/final_time << "%" << endl;
    
      //-------time marching---------     
      TPFTimeMarchSolve(fdt, num_finetimesteps, dual, areas, ppflux, ppbdrflux, fluxfuncs, satold, sat);

      count ++;
      satold = sat;

      if(count%10==0)
	{
	  sprintf(filename, "saturation_%d.vtk", count/10);
	  out.open(filename);
	  mesh->ParaviewPrint(out, satold); 
	  out.close();
	}
    }

  Array<double> l2errs(3);
  Array<double> h1errs(3);
  fem->ComputePPL2Error(l2errs, finalsat, sat, ttsat);
  fem->ComputePPH1Error(h1errs, dersol, sat, ttsat);

  cout << "Sat L2err: " << l2errs[2] << "\t" << h1errs[2] << endl;

  double maxerr = 0.0; int nn = 0;
  for(int k=0; k<dof; k++)
    if(fabs(tsat(k)-sat(k))>maxerr)
      { maxerr = fabs(tsat(k)-sat(k)); nn = k; }
  cout << nn << ": maxerr: " << maxerr << endl;
      

  //-------------free storage--------------------------
  for(int i=0; i<ttsat.Size(); i++) delete []ttsat[i];
  for(int i=0; i<flux.Size(); i++) delete []flux[i];
  for(int i=0; i<ppflux.Size(); i++) delete []ppflux[i];
  for(int i=0; i<ppsol.Size(); i++) delete []ppsol[i];
  
  delete fem;
  delete fpp;
  delete mesh;
  delete []permdata[0];
  delete []forcedata[0];
  delete totalmobility;
  delete fluxfuncs;

  delete timesat;
  delete finalsat;
  delete permy;
  delete truepressurex;
  delete truepressurey;
  delete data; 
  delete dual;
  delete tpressure;
  delete force;
  delete permfunc;
  for(int k=0; k<nbdryval.Size(); k++) delete []nbdryval[k];
  for(int k=0; k<dbdryval.Size(); k++) delete []dbdryval[k];
  cout << endl << "fin" << endl;
}
*/
