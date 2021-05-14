#include "../linalg/linalg_header.h"
#include "../fem/fem_header.h"
#include "../mesh/mesh_header.h"
#include "../msrobin/msrobin_header.h"
#include <iostream>
#include <fstream>
#include <cmath>


using namespace std;

//======================================================

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

  int Nx, Ny, nx, ny, mtype;
  input_file.getline( buffer, bufflen, ' ' );
  input_file >> Nx >> Ny >> nx >> ny >> mtype;

  //build Global Domain
  Array<int> N(2);
  N[0] = Nx*nx; N[1] = Ny*ny;

  Mesh *mesh;
  if(mtype==Element::TRIANGLE){mesh = new TMesh<2>(N, L, Element::TRIANGLE); }
  if(mtype==Element::QUADRILATERAL){ mesh = new TMesh<2>(N, L, Element::QUADRILATERAL); }

  //print
  ofstream out;
  out.open("mesh.out");
  mesh->Print(out, mtype);
  out.close();
  
  mesh->ParaviewPrint(mtype);
  

  //Build Subdomains
  int Nsd = 3*Nx*Ny/4;
  int dim = 2;
  double Hx = 1.0/((double) Nx);
  double Hy = 1.0/((double) Ny);
  cout << Nsd << endl;
  Array<msd*> subdomains(Nsd);
   
  int cind = 0;//coarse index
  for(int j=0; j<Ny; j++)
    {
      for(int i=0; i<Nx; i++)
	{ 
	  if(i<Nx/2 || j<Ny/2)
	    {
	      int type = (cind%5==0) ? 2 : 3;
	      Array<int> n(2);
	      n[0] = nx; n[1] = ny;
	      Array<double> l(4);
	      l[0] = i*Hx; l[1] = (i+1)*Hx;
	      l[2] = j*Hy; l[3] = (j+1)*Hy;
	      
	      subdomains[cind] = new msd(n, l, type);
	      cind++;
	    }
	}
    }
  
  int startindex=0;
  for(int k=0; k<Nsd; k++)
    {
      subdomains[k]->setgbvindex(startindex);
      startindex += subdomains[k]->gettotvert();
    }
  
  cout << "subdomains built" << endl;
  //---------------------------------------------------------
  // Build Data
  //--------------------------------------------------------- 

  int nnn = mesh->GetNV(); //number of vertices

  //Gather Permeability Data
  Array<double *> ellipticdata(1);
  ellipticdata[0] = new double[nnn];

  char permflag[bufflen]; 
  input_file.getline( buffer, bufflen, ' ' );
  ReadIdentifier(input_file, buffer, bufflen);
  input_file >> permflag;

  if (!strcmp(permflag, "equation"))
    {
      Function *permfunc;
      ReadIdentifier(input_file, buffer, bufflen);
      permfunc = ReadFunction(input_file);

      for(int k=0; k<nnn; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  ellipticdata[0][k] = permfunc->Eval(coord);
	}
      
      cind = 0;//coarse index
      for(int j=0; j<Ny; j++)
	{
	  for(int i=0; i<Nx; i++)
	    { 
	      if(i<Nx/2 || j<Ny/2)
		{
		  subdomains[cind]->setlocellipticdata(permfunc);	      
		  cind++;
		}
	    }
	}
      
      delete permfunc;
    }

  //Gather Force Data
  Array<double *> forcedata(1);
  forcedata[0] = new double[nnn];

  ReadIdentifier(input_file, buffer, bufflen);
  Function *force = ReadFunction(input_file);
  for(int k=0; k<nnn; k++)
    {
      double* coord = mesh->GetVertex(k);
      forcedata[0][k] = force->Eval(coord);
    }

  cind = 0;//coarse index
  for(int j=0; j<Ny; j++)
    {
      for(int i=0; i<Nx; i++)
	{ 
	  if(i<Nx/2 || j<Ny/2)
	    {
	      subdomains[cind]->setlocforcedata(force);	      
	      cind++;
	    }
	}
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
      int nn = mesh->GetNBdryElem(b)+1;
      nbdryval[b] = new double[nn];
      rbdryval[b] = new double[nn];
      rcoeff[b] = new double[nn];
      dbdryval[b] = new double[nn];

      *nbdryval[b] = 0.0;
      *rbdryval[b] = 0.0;
      *rcoeff[b] = 0.0;
      *dbdryval[b] = 0.0;
      
      input_file.getline( buffer, bufflen, ' ' );
      input_file >> bdryflag;
      
      if (!strcmp(bdryflag, "neumann"))
	{
	  BdrNeumann[b] = true;
	  BdrRobin[b] = false;
	  BdrDirichlet[b] = false;

	  ReadIdentifier(input_file, buffer, bufflen);
	  Function *bdryval = ReadFunction(input_file);
	  for(int k=0; k<nn; k++)
	    {
	      int indie = mesh->GetBdryGindex(b, k);
	      double* coord = mesh->GetVertex(indie);

	      nbdryval[b][k] = bdryval->Eval(coord);
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
	  for(int k=0; k<nn; k++)
	    {
	      int indie = mesh->GetBdryGindex(b, k);
	      double* coord = mesh->GetVertex(indie);

	      dbdryval[b][k] = bdryval->Eval(coord);
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
	  for(int k=0; k<nn; k++)
	    {
	      int indie = mesh->GetBdryGindex(b, k);
	      double* coord = mesh->GetVertex(indie);

	      rbdryval[b][k] = bdryval->Eval(coord);
	      rcoeff[b][k] = coeff->Eval(coord);
	    }

	  delete bdryval;
	  delete coeff;
	}
      
    }

  double gamma = -1.0;
  cind = 0;//coarse index
  for(int j=0; j<Ny; j++)
    {
      for(int i=0; i<Nx; i++)
	{ 
	  if(i<Nx/2 || j<Ny/2)
	    {
	      int numbdrs = subdomains[cind]->getnbdrs();	  
	      Array<double *> locnbdryval(numbdrs);
	      Array<double *> locrbdryval(numbdrs);
	      Array<double *> locrcoeff(numbdrs);
	      Array<double *> locdbdryval(numbdrs);
	      Array<bool> locBdrNeumann(numbdrs);
	      Array<bool> locBdrRobin(numbdrs);
	      Array<bool> locBdrDirichlet(numbdrs);
	      Array<bool> locGBdry(numbdrs);
	  
	      for(int b=0; b<numbdrs; b++)
		{
		  int numbv = subdomains[cind]->getnbdryelem(b)+1;
		  locnbdryval[b] = new double[numbv];
		  locrbdryval[b] = new double[numbv];
		  locrcoeff[b] = new double[numbv];
		  locdbdryval[b] = new double[numbv];

		  if((j==0 && b==0) || (j==Ny-1 && i<Nx/2 && b==2) || (j==Ny/2-1 && i>=Nx/2 && b==2) )
		    {
		      locBdrRobin[b] = BdrRobin[b];
		      locBdrNeumann[b] = BdrNeumann[b];
		      locBdrDirichlet[b] = BdrDirichlet[b];
		      locGBdry[b] = true;
		  
		      for(int k=0; k<numbv; k++)
			{
			  int dd = i*nx+k;
			  locrcoeff[b][k] = rcoeff[b][dd];
			  locrbdryval[b][k] = rbdryval[b][dd];
			  locdbdryval[b][k] = dbdryval[b][dd];
			  locnbdryval[b][k] = nbdryval[b][dd];
			}
		    }
		  else if((i==0 && b==3) || (i==Nx-1 && j<Ny/2 &&  b==1) || (i==Nx/2-1 && j>=Ny/2 && b==1))
		    {
		      locBdrRobin[b] = BdrRobin[b];
		      locBdrNeumann[b] = BdrNeumann[b];
		      locBdrDirichlet[b] = BdrDirichlet[b];
		      locGBdry[b] = true;
		  
		      for(int k=0; k<numbv; k++)
			{
			  int dd = j*ny+k;
			  locrcoeff[b][k] = rcoeff[b][dd];
			  locrbdryval[b][k] = rbdryval[b][dd];
			  locdbdryval[b][k] = dbdryval[b][dd];
			  locnbdryval[b][k] = nbdryval[b][dd];
			}

		    }
		  else
		    {
		      locBdrRobin[b] = true;
		      locBdrNeumann[b] = false;
		      locBdrDirichlet[b] = false;
		      locGBdry[b] = false;
		  
		      for(int k=0; k<numbv; k++)
			{
			  locrcoeff[b][k] = gamma;
			  locrbdryval[b][k] = 0.0;
			  locdbdryval[b][k] = 0.0;
			  locnbdryval[b][k] = 0.0;
			}
		    }
		}

	      subdomains[cind]->setlocbndrydata(locdbdryval, locnbdryval, 
						locrbdryval, locrcoeff, 
						locBdrDirichlet, locBdrNeumann,
						locBdrRobin, locGBdry);
	  
	      for(int k=0; k<numbdrs; k++)
		{
		  delete []locnbdryval[k];
		  delete []locdbdryval[k];
		  delete []locrbdryval[k];
		  delete []locrcoeff[k];
		}
	   

	      cind++;
	    }
	}
    }
  
  Data *data;
  data = new Data();
  data->SetEllipticData(ellipticdata);
  data->SetForceData(forcedata);
  data->SetNeumannData(nbdryval, BdrNeumann);
  data->SetRobinData(rcoeff, rbdryval, BdrRobin);
  data->SetDirichletData(dbdryval, BdrDirichlet);
  cout << "Data Set" << endl;

  //---------------------------------------------------------
  // Compute Multiscale Basis Functions
  //--------------------------------------------------------- 
  
  Array<int> ns(nbdry);
  input_file.getline( buffer, bufflen, ' ' );
  input_file >> ns[0] >> ns[1] >> ns[2] >> ns[3];

  cind=0;
  for(int j=0; j<Ny; j++)
    {
      for(int i=0; i<Nx; i++)
	{
	  if(i<Nx/2 || j<Ny/2)
	    {
	      subdomains[cind]->buildmsbasisfcts(ns);

	      if(i==0 && j==0){exit(2); }
	      cind++;
	    }
	}
    }

   
  //print out the mesh and basis functions  
  
  out.open("mesh.out");
  subdomains[0]->printmesh(out);
  out.close();

  int tbf = subdomains[0]->getnumbasis();
  cout << tbf << endl;
  for(int j=0; j<tbf; j++)
    {
      char fname[256];
      sprintf(fname, "irbasis_%d.out", j);
      out.open(fname);
      subdomains[0]->printbasis(out, j);
      out.close();
    }
  

  //---------------------------------------------------------
  // Set Global Index of  Multiscale Basis Functions
  //--------------------------------------------------------- 

  int gbindex=0;
  for(int k=0; k<Nsd; k++)
    {
      int numsegl = subdomains[k]->getnumseg(3);
      int numsegr = subdomains[k]->getnumseg(1);
      if(numsegl != 0)
	{
	  subdomains[k]->setgbindex(3, gbindex);
	  gbindex += numsegl+1; 
	}
      if(numsegr != 0)
	{
	  subdomains[k]->setgbindex(1, gbindex);
	  gbindex += numsegr+1; 
	}
    }

  for(int i=0; i<Nx; i++)
    {
      for(int j=0; j<Ny; j++)
	{
	  if(i<Nx/2 || j<Ny/2)
	    {
	      int Gi = (j>=Ny/2) ? Nx*Ny/2 + (j-Ny/2)*Nx/2 + i : j*Nx + i;
	      int numsegb = subdomains[Gi]->getnumseg(0);
	      int numsegt = subdomains[Gi]->getnumseg(2);
	      if(numsegb != 0)
		{
		  subdomains[Gi]->setgbindex(0, gbindex);
		  gbindex += numsegb+1; 
		}
	      if(numsegt != 0)
		{
		  subdomains[Gi]->setgbindex(2, gbindex);
		  gbindex += numsegt+1; 
		}
	    }
	}
    }

  //---------------------------------------------------------
  // Set subdomain neighbor info
  //--------------------------------------------------------- 
  
      for(int j=0; j<Ny; j++)
	{
	  for(int i=0; i<Nx; i++)
	    {
	      if(j<Ny/2)
		{
		  int k = j*Nx + i;
		  int numbdrs = subdomains[k]->getnbdrs();
		  Array<int*> nbinfo(2);  
		  for(int j=0; j<2; j++){ nbinfo[j] = new int[numbdrs]; }
	      
		  for(int b=0; b<numbdrs; b++){ nbinfo[0][b] = -1; nbinfo[1][b] = -1; }
		  if(subdomains[k]->getGBdry(0)==false){ nbinfo[0][0] = k-Nx; nbinfo[1][0]=2; }
		  if(subdomains[k]->getGBdry(1)==false){ nbinfo[0][1] = k+1;  nbinfo[1][1]=3; }
		  if(subdomains[k]->getGBdry(2)==false){ nbinfo[0][2] = k+Nx; nbinfo[1][2]=0; }
		  if(subdomains[k]->getGBdry(3)==false){ nbinfo[0][3] = k-1;  nbinfo[1][3]=1; }
	      
		  subdomains[k]->setnbinfo(nbinfo);
	      
		  for(int j=0; j<2; j++){ delete []nbinfo[j]; }  	  
		}
	      if(i<Nx/2 && j>=Ny/2)
		{
		  int k = Nx*Ny/2 + (j-Ny/2)*Nx/2 + i;
		  int numbdrs = subdomains[k]->getnbdrs();
		  Array<int*> nbinfo(2);  
		  for(int j=0; j<2; j++){ nbinfo[j] = new int[numbdrs]; }
	      
		  int Nxx = (j==Ny/2) ? Nx : Nx/2;
		  for(int b=0; b<numbdrs; b++){ nbinfo[0][b] = -1; nbinfo[1][b] = -1; }
		  if(subdomains[k]->getGBdry(0)==false){ nbinfo[0][0] = k-Nxx; nbinfo[1][0]=2; }
		  if(subdomains[k]->getGBdry(1)==false){ nbinfo[0][1] = k+1;    nbinfo[1][1]=3; }
		  if(subdomains[k]->getGBdry(2)==false){ nbinfo[0][2] = k+Nx/2; nbinfo[1][2]=0; }
		  if(subdomains[k]->getGBdry(3)==false){ nbinfo[0][3] = k-1;    nbinfo[1][3]=1; }
	      
		  subdomains[k]->setnbinfo(nbinfo);
	      
		  for(int j=0; j<2; j++){ delete []nbinfo[j]; }  

		}
		}
	    }

  //---------------------------------------------------------
  //  Read and create exact solution if known.
  //--------------------------------------------------------- 

  Array<Function *> exactderfunc;
  Function *exactsolfunc;
  ReadIdentifier(input_file, buffer, bufflen);
  if (!strcmp(buffer, "exactsolknown"))
    {
      ReadIdentifier(input_file, buffer, bufflen);
      exactsolfunc = ReadFunction(input_file);
      exactderfunc.SetSize(2);
      for (int i=0; i<2; i++)
	{
	  ReadIdentifier(input_file, buffer, bufflen);
	  exactderfunc[i] = ReadFunction(input_file);
	}
    }

  // Close the main input file.
  input_file.close();

  cout << "Exact Solution Set" << endl;

  //---------------------------------------------------------
  //  Build MsRM
  //--------------------------------------------------------- 
  directMsRM msolver(subdomains, gbindex, Nx, Ny);
  msolver.paraviewprintmesh();
  msolver.solve();
  // Vector msol;
  //msolver.buildsol(msol, Nx, Ny, nx, ny);

  //---------------------------------------------------------
  //  Build FEM 
  //--------------------------------------------------------- 
  FEM *fem;
  if(mtype==Element::QUADRILATERAL){ fem = new BilinearFEM2D(mesh, data); }
  if(mtype==Element::TRIANGLE){ fem = new LinearFEM2D(mesh, data); }
  fem->Assemble();
  cout << "FEM Built" << endl;

  //---------------------------------------------------------
  //  Solve
  //---------------------------------------------------------
  Vector sol(nnn);
  fem->Solve(sol, "multifrontal");
  cout << "FEM Problem Solved" << endl;  

  //---------------------------------------------------------
  //  Print
  //---------------------------------------------------------
  
  out.open("sol.out");
  mesh->Print(out, sol);
  out.close();
  
  mesh->ParaviewPrint(sol, mtype);
  msolver.paraviewprintsol();
  msolver.paraviewprintellipticdata();

  cout << endl;
  cout << "msrobin: " << endl;
  cout << "L2error: " << msolver.computel2error(exactsolfunc) << endl;
  cout << "H1error: " << msolver.computeh1error(exactderfunc) << endl;

  cout << endl << "fem: " << endl;
  cout << "L2error: " << fem->ComputeL2Error(exactsolfunc, sol) << endl;
  cout << "H1error: " << fem->ComputeH1Error(exactderfunc, sol) << endl;

  //---------------------------------------------------------
  //  Freedom!!!!
  //---------------------------------------------------------
  delete mesh;
  delete []ellipticdata[0];
  delete []forcedata[0];
  delete data;
  delete fem;
  delete exactsolfunc;

  for (int k=0; k<exactderfunc.Size(); k++)
    delete exactderfunc[k]; 

  for(int k=0; k<nbdryval.Size(); k++)
    delete []nbdryval[k];

  for(int k=0; k<dbdryval.Size(); k++)
    delete []dbdryval[k];

  for(int k=0; k<rbdryval.Size(); k++)
    delete []rbdryval[k];

  for(int k=0; k<rcoeff.Size(); k++)
    delete []rcoeff[k];

  for(int k=0; k<subdomains.Size(); k++)
    delete subdomains[k];

  cout << endl << "fin" << endl;
}

