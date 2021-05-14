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
  ofstream out; char filename[256];

  //-----------------------------
  // Build the Global Mesh
  //-----------------------------
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
  int anx = aN[0] + 1;
  
  Mesh *mesh = new TMesh<2>(N, L, Element::TRIANGLE);
  Mesh *amesh = new TMesh<2>(aN, L, Element::TRIANGLE);

  cout << endl << "Mesh Size: " << N[0] << "\t" << N[1] << endl;

  int nv = mesh->GetNV(); //number of vertices
  int nbe = mesh->GetNBE();
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

  cout << endl << "Final Time: " << final_time << endl;
  cout << "Course and Fine Time Steps: " << num_timesteps << "\t" << num_finetimesteps << endl << endl;

  ReadIdentifier(input_file, buffer, bufflen);
  Function *totalmobility = ReadFunction(input_file);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *fluxfuncs = ReadFunction(input_file);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *timesat = ReadFunction(input_file);

  // Close the main input file.
  input_file.close();

  double dt = final_time/(double)num_timesteps;
  double fdt = dt/(double)num_finetimesteps;
  Vector satold(dof);
  Vector sat(dof);
  Vector temsat(dof);
  Vector sol(dof);
  sol = 0.0;
  satold = 0.0;
  for(int k=0; k<dof; k++)
    {
      double* coord = amesh->GetVertex(k);
      double tcoord[2];
      tcoord[0] = coord[0];
      tcoord[1] = coord[1];
      satold(k) = timesat->Eval(tcoord);
    }
  //for(int j=0; j<=aN[1]; j++){ satold(j*(aN[0]+1)) = 1.0; }
  sat = satold;
  temsat = satold;

  Vector initperm(dof);
  for(int i=0; i<dof; i++){ initperm(i) =  permdata[0][i]; }

  Vector perm(nv);   
  for(int i=0; i<nv; i++){perm(i) = log(initperm(i)); }   
  sprintf(filename, "perm.vtk");
  out.open(filename);
  mesh->ParaviewPrint(out, perm); 
  out.close();  

  DualMesh *adual = new DualMesh(amesh); 
  DualMesh *dual = new DualMesh(mesh, 2); 
  Array<double> areas(dof);
  for(int i=0; i<dof; i++) areas[i] = 0.0;
      
  FEM *temfem = new QuadraticFEM2D(mesh, data);
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<double> locareas(1);  
      Array<int> ind;
      temfem->getElementDOF(i, ind);
	  
      dual->GetElemDualAreas(i, locareas);
	  
      for(int j=0; j<3; j++) areas[ind[j]] += locareas[0];
      for(int j=3; j<6; j++) areas[ind[j]] += 3.0*locareas[0];
    } 

  Array<double *> flux(9);
  Array<double *> ppsol(6);
  Array<double *> ppflux(9);
  for(int i=0; i<flux.Size(); i++) flux[i] = new double[mesh->GetNE()];
  for(int i=0; i<ppsol.Size(); i++) ppsol[i] = new double[mesh->GetNE()];
  for(int i=0; i<ppflux.Size(); i++) ppflux[i] = new double[mesh->GetNE()];

  Array<double> ppbdrflux(dof);

  Array<double *> appflux(3);
  for(int i=0; i<appflux.Size(); i++) appflux[i] = new double[amesh->GetNE()];
  Array<double> appbdrflux(dof);

  int count = 0;
  while(count < num_timesteps)
    {
      //cout << "Two Phase Flow Simulation Time Marching: " << 100*count*dt/final_time << "%" << endl;

      //---QuadraticFEM for Pressure-----------------------------
      for(int i=0; i<dof; i++)
	{
	  double sval[1];
	  sval[0] = temsat(i);
	  double lambda = totalmobility->Eval(sval);
	  permdata[0][i] = lambda*initperm(i);
	}   
      data->SetEllipticData(permdata);
    
      FEM *fem = new QuadraticFEM2D(mesh, data);
      fem->Assemble();
      fem->FinalizeMatrix();
      fem->Solve(sol, "multifrontal");

      //---------------------------------------------------------
      //  Post-Process
      //---------------------------------------------------------
      cout << "Computing Conservative Flux..." << endl;
      Postprocessing *fpp = new Postprocessing(dual, fem, data);
      fpp->ComputeConservativeFlux(truepressure, sol, flux, ppsol, ppflux, ppbdrflux);

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

      //---Time Marching-----------------------------------------
      TPFTimeMarchSolve(fdt, num_finetimesteps, adual, areas, appflux, appbdrflux, fluxfuncs, satold, sat);
      satold = sat;
      count ++;

      /*
      if(count%1==0)
	{
	  cout << "Two Phase Flow Simulation Time Marching: " << 100*count*dt/final_time << "%" << endl;
	  sprintf(filename, "saturation_%d.vtk", count/1);
	  out.open(filename);
	  amesh->ParaviewPrint(out, sat); 
	  out.close();

	  //sprintf(filename, "pressure_%d.vtk", count/10);
	  //out.open(filename);
	  //mesh->ParaviewPrint(out, sol); 
	  //out.close();  
	}
      */

      for(int j=0; j<N[1]; j++)
	{
	  for(int i=0; i<N[0]; i++)
	    {
	      Array<int> indl;
	      Array<int> indu;
	      fem->getElementDOF(2*(j*N[0]+i), indl);
	      fem->getElementDOF(2*(j*N[0]+i)+1, indu);
	      int rind = 2*j*anx + 2*i;

	      temsat(indl[0]) = sat(rind);
	      temsat(indl[1]) = sat(rind+2);
	      temsat(indl[2]) = sat(rind+2*anx);
	      temsat(indl[3]) = sat(rind+1);
	      temsat(indl[4]) = sat(rind+1+anx);
	      temsat(indl[5]) = sat(rind+anx);

	      temsat(indu[1]) = sat(rind+2+2*anx);
	      temsat(indu[3]) = sat(rind+2+anx);
	      temsat(indu[4]) = sat(rind+1+2*anx);
	    }
	}
      delete fpp;
      delete fem;
    }
  cout << "--------Time Marching Done------------" << endl << endl;

  /*
  sprintf(filename, "ex21srefqsat_%d", aN[0]);
  out.open(filename);
  for(int i=0; i<dof; i++)
    {
      out << sat(i) << endl; 
    }
  out.close(); 
  */

  
  //----L2 Error---------------
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
  sprintf(refsatfile, "ex22srefqsat_256");
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
	  ind[0] = j*anx + i;
	  ind[1] = ind[0] + 1;
	  ind[2] = ind[1] + anx;
	  ind[3] = ind[0] + anx;

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
      if(usat(k)>1) cout << usat(k) << endl;
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
  for(int i=0; i<appflux.Size(); i++) delete []appflux[i];
  
  delete temfem;     
  delete mesh;
  delete []permdata[0];
  delete []forcedata[0];
  delete totalmobility;
  delete fluxfuncs;
  delete timesat;
  delete data; 
  delete adual;
  delete amesh;
  delete dual;
  delete tpressure;
  delete force;
  delete permfunc;
  for(int k=0; k<nbdryval.Size(); k++) delete []nbdryval[k];
  for(int k=0; k<dbdryval.Size(); k++) delete []dbdryval[k];
  cout << endl << "fin" << endl;
}



/*

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
  ofstream out; char filename[256];

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

  cout << endl << "Mesh: " << "\t" << N[0] << "\t" << N[1] << endl;

  int nv = mesh->GetNV(); //number of vertices
  int nbe = mesh->GetNBE();
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

  cout << endl << "Final Time: " << final_time << endl;
  cout << "Course and Fine Time Steps: " << num_timesteps << "\t" << num_finetimesteps << endl << endl;

  ReadIdentifier(input_file, buffer, bufflen);
  Function *totalmobility = ReadFunction(input_file);

  ReadIdentifier(input_file, buffer, bufflen);
  Function *fluxfuncs = ReadFunction(input_file);

  // Close the main input file.
  input_file.close();

  double dt = final_time/(double)num_timesteps;
  double fdt = dt/(double)num_finetimesteps;
  Vector satold(dof);
  Vector sat(dof);
  Vector sol(dof);
  sol = 0.0;
  satold = 0.0;
  for(int j=0; j<=N[1]; j++){ satold(j*(N[0]+1)) = 1.0; }
  for(int j=0; j<nbe/4; j++)
    {
      int temp = mesh->GetBdrElementEdgeIndex(j+3*nbe/4);
      satold(temp+nv) = 1.0;
    }
  sat = satold;
  sol = sat;


  Vector initperm(dof);
  for(int i=0; i<dof; i++){ initperm(i) =  permdata[0][i]; }

  Vector perm(nv);   
  for(int i=0; i<nv; i++){perm(i) = log(initperm(i)); }   
  sprintf(filename, "perm.vtk");
  out.open(filename);
  mesh->ParaviewPrint(out, perm); 
  out.close();  
  
  DualMesh *dual = new DualMesh(mesh, 2); 
  Array<double> areas(dof);
  for(int i=0; i<dof; i++) areas[i] = 0.0;
      
  FEM *temfem = new QuadraticFEM2D(mesh, data);
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<double> locareas(1);  
      Array<int> ind;
      temfem->getElementDOF(i, ind);
	  
      dual->GetElemDualAreas(i, locareas);
	  
      for(int j=0; j<3; j++) areas[ind[j]] += locareas[0];
      for(int j=3; j<6; j++) areas[ind[j]] += 3.0*locareas[0];
    } 

  Array<double *> flux(9);
  Array<double *> ppsol(6);
  Array<double *> ppflux(9);
  for(int i=0; i<flux.Size(); i++) flux[i] = new double[mesh->GetNE()];
  for(int i=0; i<ppsol.Size(); i++) ppsol[i] = new double[mesh->GetNE()];
  for(int i=0; i<ppflux.Size(); i++) ppflux[i] = new double[mesh->GetNE()];

  Array<double> ppbdrflux(dof);
  
  int count = 1;
  while(count <= num_timesteps)
    {
      cout << "Two Phase Flow Simulation Time Marching: " << 100*count*dt/final_time << "%" << endl;

      //---------------------------------------------------------
      //  Print
      //---------------------------------------------------------
      sprintf(filename, "saturation_%d.vtk", count);
      out.open(filename);
      mesh->ParaviewPrint(out, satold); 
      out.close();

      sprintf(filename, "pressure_%d.vtk", count);
      out.open(filename);
      mesh->ParaviewPrint(out, sol); 
      out.close();  

      //---------------------------------------------------------
      //  QuadraticFEM for Pressure
      //--------------------------------------------------------- 
      for(int i=0; i<dof; i++)
	{
	  double sval[1];
	  sval[0] = sat(i);
	  double lambda = totalmobility->Eval(sval);
	  permdata[0][i] = lambda*initperm(i);
	}   
      data->SetEllipticData(permdata);
    
      FEM *fem = new QuadraticFEM2D(mesh, data);
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

      count ++;
      delete fpp;
      delete fem;
    }
  cout << "--------Time Marching Done------------" << endl << endl;

  int nx = 128;
  int ny = 128;
  int rx = nx/N[0];
  int ry = ny/N[1];

  int rdof = (nx+1)*(ny+1);
  Vector rsat(rdof);
  Vector usat(rdof);
  rsat = 0.0;
  usat = 0.0;

  char refsatfile[256];
  sprintf(refsatfile, "refsat_128");
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
	  Array<int> indl;
	  Array<int> indu;
	  temfem->getElementDOF(2*(j*N[0]+i), indl);
	  temfem->getElementDOF(2*(j*N[0]+i)+1, indu);

	  rind0 = j*ry*(nx+1) + i*rx;
	  for(int jj=0; jj<ry+1; jj++)
	    {
	      for(int ii=0; ii<rx+1; ii++)
		{
		  rind = rind0 + jj*(nx+1) + ii;
		  usat(rind) = sat(indl[0]) * ( 1.0 - ii/double(rx) ) * ( 1.0 - 2.0*ii/double(rx) ) * ( 1.0 - jj/double(ry) ) * ( 1.0 - 2.0*jj/double(ry) ) 
		    + sat(indl[1]) * ( - ii/double(rx) ) * ( 1.0 - 2.0*ii/double(rx) ) * ( 1.0 - jj/double(ry) ) * ( 1.0 - 2.0*jj/double(ry) ) 
		    + sat(indu[1]) * ( - ii/double(rx) ) * ( 1.0 - 2.0*ii/double(rx) ) * ( - jj/double(ry) ) * ( 1.0 - 2.0*jj/double(ry) ) 
		    + sat(indl[2]) * ( 1.0 - ii/double(rx) ) * ( 1.0 - 2.0*ii/double(rx) ) * ( - jj/double(ry) ) * ( 1.0 - 2.0*jj/double(ry) ) 
		    + sat(indl[3]) * ( 1.0 - ii/double(rx) ) * ( 4.0*ii/double(rx) ) * ( 1.0 - jj/double(ry) ) * ( 1.0 - 2.0*jj/double(ry) ) 
		    + sat(indu[3]) * ( - ii/double(rx) ) * ( 1.0 - 2.0*ii/double(rx) ) * ( 1.0 - jj/double(ry) ) * ( 4.0*jj/double(ry) ) 
		    + sat(indu[4]) * ( 1.0 - ii/double(rx) ) * ( 4.0*ii/double(rx) ) * ( - jj/double(ry) ) * ( 1.0 - 2.0*jj/double(ry) ) 
		    + sat(indl[5]) * ( 1.0 - ii/double(rx) ) * ( 1.0 - 2.0*ii/double(rx) ) * ( 1.0 - jj/double(ry) ) * ( 4.0*jj/double(ry) ) 
		    + sat(indl[4]) * ( 1.0 - ii/double(rx) ) * ( 4.0*ii/double(rx) ) * ( 1.0 - jj/double(ry) ) * ( 4.0*jj/double(ry) );
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
  
  delete temfem;     
  delete mesh;
  delete []permdata[0];
  delete []forcedata[0];
  delete totalmobility;
  delete fluxfuncs;

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
