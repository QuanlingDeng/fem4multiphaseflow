#include <cmath>
#include <fstream>
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

//==============================================================

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
  Mesh *mesh = new TMesh<2>(N, L, Element::TRIANGLE, true);

  cout << endl << "Mesh: " << "\t" << N[0] << "\t" << N[1] << endl;

  //print
  ofstream out;
  char filename[256];
  //out.open("mesh.out"); mesh->Print(out, meshtype);  out.close();  mesh->ParaviewPrint(meshtype);
   
  //---------------------------------------------------------
  // Build Data
  //--------------------------------------------------------- 

  int nv = mesh->GetNV(); //number of vertices

  //Gather Permeability Data
  Array<double *> permdata(1);
  permdata[0] = new double[nv];

  char permflag[bufflen]; 
  ReadIdentifier(input_file, buffer, bufflen);
  input_file >> permflag;
  Function *permfunc;
  
  if (!strcmp(permflag, "equation"))
    {
      ReadIdentifier(input_file, buffer, bufflen);
      permfunc = ReadFunction(input_file);

      for(int k=0; k<nv; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  permdata[0][k] = permfunc->Eval(coord);
	}      
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
  Array<double *> forcedata(1);
  forcedata[0] = new double[nv];

  ReadIdentifier(input_file, buffer, bufflen);
  Function *force = ReadFunction(input_file);

  Vector truepressure(nv);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *tpressure = ReadFunction(input_file);

  for(int k=0; k<nv; k++)
    {
      double* coord = mesh->GetVertex(k);
      forcedata[0][k] = force->Eval(coord);
      truepressure(k) = tpressure->Eval(coord);
    }

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

  double final_time;
  int num_timesteps;
  int num_finetimesteps;
  input_file.getline( buffer, bufflen, ' ' );
  input_file >> final_time >> num_timesteps >> num_finetimesteps;

  cout << endl << "Final Time: " << final_time << endl;
  cout << "Course and Fine Time Steps: " << num_timesteps << "\t" << num_finetimesteps << endl << endl;

  ReadIdentifier(input_file, buffer, bufflen);
  Function *totalmobility = ReadFunction(input_file);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *fluxfuncs = ReadFunction(input_file);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *timesat = ReadFunction(input_file);

  Data *data;
  data = new Data();
  data->SetEllipticData(permdata);
  data->SetEllipticFunction(permfunc);
  data->SetForceData(forcedata);
  data->SetForceFunction(force);
  data->SetNeumannData(nbdryval, BdrNeumann);
  data->SetRobinData(rcoeff, rbdryval, BdrRobin);
  data->SetDirichletData(dbdryval, BdrDirichlet);
  cout << "Data Set" << endl;
  
  ReadIdentifier(input_file, buffer, bufflen);
  Function *truepressurex = ReadFunction(input_file);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *truepressurey = ReadFunction(input_file);

  Array<Function *> dersol(2);
  dersol[0] = truepressurex;
  dersol[1] = truepressurey;
  
  // Close the main input file.
  input_file.close();

  Vector initperm(nv);
  for(int i=0; i<nv; i++){ initperm(i) =  permdata[0][i]; }

  Vector perm(nv);   
  for(int i=0; i<nv; i++){ perm(i) = log(initperm(i)); }   

  sprintf(filename, "perm.vtk");
  out.open(filename);
  mesh->ParaviewPrint(out, perm); 
  out.close();  

  /*
  sprintf(filename, "perm.z");
  out.open(filename);
  out << "! nx " << N[0]+1 << " ny " << N[1]+1 << " xmin " << 0 << " xmax " << L[1] << " ymin " << 0 << " ymax " << L[3] << endl;
  out << setprecision(6) << setiosflags(ios::scientific | ios::showpos);
  for(int l=0; l<(N[1]+1); l++)
    {
      for(int k=0; k<(N[0]+1); k++) 
	{
	  int index = k+l*(N[0]+1);
	  out << log(permdata[0][index]) << " "; 
	}
      out << endl;
    }
  out.close(); 
  */

  int nl = 3;
  Array<double *> flux(nl);
  Array<double *> ppsol(nl);
  Array<double *> ppflux(nl);
  for(int i=0; i<nl; i++) flux[i] = new double[mesh->GetNE()];
  for(int i=0; i<nl; i++) ppsol[i] = new double[mesh->GetNE()];
  for(int i=0; i<nl; i++) ppflux[i] = new double[mesh->GetNE()];
  Array<double> ppbdrflux(nv);

  double dt = final_time/(double)num_timesteps;
  double fdt = dt/(double)num_finetimesteps;
  Vector satold(nv);
  Vector sat(nv);
  Vector sol(nv);
  sol = 0.0;
  satold = 0.0;
  for(int k=0; k<nv; k++)
    {
      double* coord = mesh->GetVertex(k);      
      double tcoord[2];
      tcoord[0] = coord[0];
      tcoord[1] = coord[1];
      satold(k) = timesat->Eval(tcoord);
    }
  //for(int j=0; j<=N[1]; j++){ satold(j*(N[0]+1)) = 1.0; }
  sat = satold;

  DualMesh *dual = new DualMesh(mesh); 
  Array<double> areas(nv);
  for(int i=0; i<nv; i++) areas[i] = 0.0;
      
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<double> locareas(3);  
      Array<int> ind;
      mesh->GetElementVertices(i, ind);
	  
      dual->GetElemDualAreas(i, locareas);
	  
      for(int j=0; j<3; j++) areas[ind[j]] += locareas[j];
    }
      
  int count = 1;
  while(count <= num_timesteps)
    {
      //cout << "Two Phase Flow Simulation Time Marching: " << 100*count*dt/final_time << "%" << endl;
      
      /*/---------------------------------------------------------
      //  Print
      //---------------------------------------------------------
      if(count%1==0)
	{
	  sprintf(filename, "saturation_%d.vtk", count/1);
	  out.open(filename);
	  mesh->ParaviewPrint(out, satold); 
	  out.close();

	  sprintf(filename, "pressure_%d.vtk", count/1);
	  out.open(filename);
	  mesh->ParaviewPrint(out, sol); 
	  out.close();  
	}
      */
      
      //---------------------------------------------------------
      //  LinearFEM for Pressure
      //--------------------------------------------------------- 
      for(int i=0; i<nv; i++)
	{
	  double sval[1];
	  sval[0] = sat(i);
	  double lambda = totalmobility->Eval(sval);
	  permdata[0][i] = lambda*initperm(i);
	}   
      data->SetEllipticData(permdata);
    
      FEM *fem = new LinearFEM2D(mesh, data);
      fem->Assemble();
      fem->FinalizeMatrix();
      fem->Solve(sol, "multifrontal");

      //---------------------------------------------------------
      //  Post-Process
      //---------------------------------------------------------
      Postprocessing *fpp = new Postprocessing(dual, fem, data);

      fpp->ComputeConservativeFlux(truepressure, sol, flux, ppsol, ppflux, ppbdrflux);

      //---------------------------------------------------------
      //  Upwind
      //---------------------------------------------------------

      TPFTimeMarchSolve(fdt, num_finetimesteps, dual, areas, ppflux, ppbdrflux, fluxfuncs, satold, sat);
            
      satold = sat;

      /*
      if(count%2000==0)
	{
	  sprintf(filename, "sat_%f.z", currtime);
	  out.open(filename);
	  out << "! nx " << N[0]+1 << " ny " << N[1]+1 << " xmin " << 0 << " xmax " << L[1]/660 << " ymin " << 0 << " ymax " << L[3]/180 << endl;
	  out << setprecision(6) << setiosflags(ios::scientific | ios::showpos);
	  for(int l=0; l<(N[1]+1); l++)
	    {
	      for(int k=0; k<(N[0]+1); k++) 
		{
		  int index = k+l*(N[0]+1);
		  out << sat(index) << " "; 
		}
	      out << endl;
	    }
	  out.close(); 
	}
      */

      count ++;
      delete fpp;
      delete fem;
    }
  cout << "--------Time Marching Done------------" << endl << endl;

  /*
  sprintf(filename, "ex21refsat_%d", N[0]);
  out.open(filename);
  for(int i=0; i<nv; i++)
    {
      out << sat(i) << endl; 
    }
  out.close(); 
  */

  
  int nx = 256;
  int ny = 256;
  int rx = nx/N[0];
  int ry = ny/N[1];

  int rdof = (nx+1)*(ny+1);
  Vector rsat(rdof);
  Vector usat(rdof);
  rsat = 0.0;
  usat = 0.0;

  char refsatfile[256];
  sprintf(refsatfile, "ex22srefqsat_256");
  input_file.open( refsatfile );
  for(int i=0; i<rdof; i++)
    {
      input_file >> rsat(i);
    }
  input_file.close();

  int rind, rind0;
  for(int j=0; j<N[1]; j++)
    {
      for(int i=0; i<N[0]; i++)
	{
	  Array<int> ind(4);
	  ind[0] = j*(N[0]+1) + i;
	  ind[1] = ind[0] + 1;
	  ind[2] = ind[1] + N[0] + 1;
	  ind[3] = ind[0] + N[0] + 1;

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
  

  //---------------------------------------------------------
  //  Freedom!!!!
  //---------------------------------------------------------
  for(int i=0; i<flux.Size(); i++) delete []flux[i];
  for(int i=0; i<ppflux.Size(); i++) delete []ppflux[i];
  for(int i=0; i<ppsol.Size(); i++) delete []ppsol[i];
  
  delete mesh;
  delete []permdata[0];
  delete []forcedata[0];
  delete totalmobility;
  delete fluxfuncs;
  delete timesat;
  delete truepressurex;
  delete truepressurey;
  delete data; 
  delete dual;
  delete tpressure;
  delete force;
  delete permfunc;
  for(int k=0; k<nbdryval.Size(); k++) delete []nbdryval[k];
  for(int k=0; k<dbdryval.Size(); k++) delete []dbdryval[k];
  for(int k=0; k<rbdryval.Size(); k++) delete []rbdryval[k];
  for(int k=0; k<rcoeff.Size(); k++) delete []rcoeff[k];
  cout << endl << "fin" << endl;
}
