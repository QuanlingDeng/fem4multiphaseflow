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
  if(mtype==Element::TRIANGLE){mesh = new TMesh<2>(N, L, Element::TRIANGLE); }//triangle = 2
  if(mtype==Element::QUADRILATERAL){ mesh = new TMesh<2>(N, L, Element::QUADRILATERAL); }//quadrilateral=3

  //print
  ofstream out;
  out.open("mesh.out");
  mesh->Print(out, mtype);
  out.close();
  
  mesh->ParaviewPrint(mtype);
  

  //Build Subdomains
  int Nsd = Nx*Ny;
  int dim = 2;
  double Hx = (L[1] - L[0])/((double) Nx);
  double Hy = (L[3] - L[2])/((double) Ny);
  Array<msd*> subdomains(Nsd);
  
  int nbdrcon = 2; 
  int cind = 0;//coarse index
  for(int j=0; j<Ny; j++)
    {
      for(int i=0; i<Nx; i++)
	{ 
	  int type = mtype; //(cind%3==0) ? 2 : 3;
	  Array<int> n(2);
	  n[0] = nx; n[1] = ny;
	  Array<double> l(4);
	  l[0] = i*Hx; l[1] = (i+1)*Hx;
	  l[2] = j*Hy; l[3] = (j+1)*Hy;
	  
	  subdomains[cind] = new msd(n, l, type, nbdrcon);
	  cind++;
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
  int dof = nbdrcon*nnn; //number of basis functions

  //Gather Permeability Data
  Array<double *> ellipticdata(nbdrcon);
  ellipticdata[0] = new double[nnn];
  ellipticdata[1] = new double[nnn];

  char elmodflag[bufflen]; 
  input_file.getline( buffer, bufflen, ' ' );
  ReadIdentifier(input_file, buffer, bufflen);
  input_file >> elmodflag;

  if (!strcmp(elmodflag, "equation"))
    {
      Function *modfunc;
      ReadIdentifier(input_file, buffer, bufflen);
      modfunc = ReadFunction(input_file);
      
      double prat = 0.25;
      for(int k=0; k<nnn; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  double data = modfunc->Eval(coord);
	  ellipticdata[0][k] = data*prat/((1.0+prat)*(1.0-2.0*prat));
	  ellipticdata[1][k] = data/(1+prat);
	}  
    
      cind = 0;//coarse index
      for(int j=0; j<Ny; j++)
	{
	  for(int i=0; i<Nx; i++)
	    { 
	      subdomains[cind]->setmaterialdata(modfunc,prat);	      
	      cind++;
	    }
	}
      delete modfunc;
    }
  else
    {
      char modfile[bufflen];
      double sigma, prat;
      ifstream modin;
      input_file >> modfile >> sigma >> prat;
      modin.open(modfile);
      if(!modin) { cout << modfile << " does not exist.\n";  exit(1); }
      double val, data;
      for (int k=0; k<nnn; k++)
	{
	  modin >> val; 
	  data = exp(sigma * val);
	  ellipticdata[0][k] = data*prat/((1.0+prat)*(1.0-2.0*prat));
	  ellipticdata[1][k] = data/(1+prat);
	}
      modin.close();

      for (int j=0; j<Ny; j++)
	{
	  for (int i=0; i<Nx; i++)
	    {
	      int Gind = j*Nx + i;
	      int nv = subdomains[Gind]->gettotvert();
	      Array<double *> tempellipticdata(nbdrcon);
	      tempellipticdata[0] = new double[nv];
	      tempellipticdata[1] = new double[nv];
	      for (int k=0; k<=ny; k++)
		{
		  for (int l=0; l<=nx; l++)
		    {
		      int FGind = l + i*nx + (k + j*ny)*(Nx*nx + 1);
		      int ind =  k*(nx+1) + l;
		      tempellipticdata[0][ind] = ellipticdata[0][FGind];
		      tempellipticdata[1][ind] = ellipticdata[1][FGind];
		    }
		}
	      subdomains[Gind]->setmaterialdata(tempellipticdata);

	      for(int k=0; k<tempellipticdata.Size(); k++)
		{
		  delete []tempellipticdata[k];
		}
	    }
	}

    }


  //Gather Force Data
  Array<double *> forcedata(nbdrcon);
  forcedata[0] = new double[nnn];
  forcedata[1] = new double[nnn];


  ReadIdentifier(input_file, buffer, bufflen);
  Function *forcex = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *forcey = ReadFunction(input_file);

  for(int k=0; k<nnn; k++)
    {
      double* coord = mesh->GetVertex(k);
      forcedata[0][k] = forcex->Eval(coord);
      forcedata[1][k] = forcey->Eval(coord);
    }


  cind = 0;//coarse index
  for(int j=0; j<Ny; j++)
    {
      for(int i=0; i<Nx; i++)
	{ 
	  subdomains[cind]->setlocforcedata(forcex,forcey);	      
	  cind++;
	}
    }
  delete forcex;
  delete forcey;

  
  double deltaforce;
  input_file.getline( buffer, bufflen, ' ' );
  input_file >> deltaforce;
  
  Array<double> nodevalue(1);
  Array<int> nodeindex(1);
  
  nodevalue[0] = deltaforce;
  nodeindex[0] = nnn+N[1]*(N[0]+1)/2+N[0]/2;
  
  cind = (Ny-1)*(Nx+1)/2+Nx/2;
  cout << cind << endl;
  Array<int> locnodeindex(1);
  locnodeindex[0] = subdomains[cind]->gettotvert() + ny*(nx+1)/2+nx/2;
  subdomains[cind]->setpointsourcedata(nodevalue, locnodeindex);

  cout  << locnodeindex[0] << " " << nodeindex[0] << endl;

  //Gather Boundary Data
  Function *gammafunc;
  ReadIdentifier(input_file, buffer, bufflen);
  gammafunc = ReadFunction(input_file);

  int nbredges = mesh->GetNBdrs();
  int nbdry = nbdrcon*nbredges;
  Array<double *> nbdryval(nbdry);
  Array<double *> rbdryval(nbdry);
  Array<double *> rcoeff(nbdry);
  Array<double *> dbdryval(nbdry);
  Array<bool> BdrNeumann(nbdry);
  Array<bool> BdrRobin(nbdry);
  Array<bool> BdrDirichlet(nbdry);

  char bdryflag[bufflen]; 
  for(int b=0; b<nbredges; b++)
    {
      int nv = mesh->GetNBdryElem(b)+1;
      for(int i=0; i<nbdrcon; i++)
	{
	  int bb = b+i*nbredges;
	  nbdryval[bb] = new double[nv];
	  rbdryval[bb] = new double[nv];
	  rcoeff[bb] = new double[nv];
	  dbdryval[bb] = new double[nv];
	  
	  input_file.getline( buffer, bufflen, ' ' );
	  input_file >> bdryflag;

	  
	  if (!strcmp(bdryflag, "neumann"))
	    {
	      BdrNeumann[bb] = true;
	      BdrRobin[bb] = false;
	      BdrDirichlet[bb] = false;
	      
	      ReadIdentifier(input_file, buffer, bufflen);
	      Function *bdryval = ReadFunction(input_file);
	      for(int k=0; k<nv; k++)
		{
		  int indie = mesh->GetBdryGindex(b, k);
		  double* coord = mesh->GetVertex(indie);
		  nbdryval[bb][k] = bdryval->Eval(coord);
		  
		  rbdryval[bb][k] = 0.0;
		  rcoeff[bb][k] = 0.0;
		  dbdryval[bb][k] = 0.0;
		}

	      delete bdryval;
	    }
	  
	  if (!strcmp(bdryflag, "dirichlet"))
	    {
	      BdrDirichlet[bb] = true;
	      BdrNeumann[bb] = false;
	      BdrRobin[bb] = false;

	      ReadIdentifier(input_file, buffer, bufflen);
	      Function *bdryval = ReadFunction(input_file);
	      for(int k=0; k<nv; k++)
		{
		  int indie = mesh->GetBdryGindex(b, k);
		  double* coord = mesh->GetVertex(indie);

		  dbdryval[bb][k] = bdryval->Eval(coord);

		  rbdryval[bb][k] = 0.0;
		  rcoeff[bb][k] = 0.0;
		  nbdryval[bb][k] = 0.0;
		}
	      delete bdryval;
	    }
      
	  if (!strcmp(bdryflag, "robin"))
	    {
	      BdrRobin[bb] = true;
	      BdrNeumann[bb] = false;
	      BdrDirichlet[bb] = false;

	      ReadIdentifier(input_file, buffer, bufflen);
	      Function *bdryval = ReadFunction(input_file);
	      ReadIdentifier(input_file, buffer, bufflen);
	      Function *coeff = ReadFunction(input_file);
	      for(int k=0; k<nv; k++)
		{
		  int indie = mesh->GetBdryGindex(b, k);
		  double* coord = mesh->GetVertex(indie);

		  rbdryval[bb][k] = bdryval->Eval(coord);
		  rcoeff[bb][k] = coeff->Eval(coord);

		  dbdryval[bb][k] = 0.0;
		  nbdryval[bb][k] = 0.0;
		}
	      delete bdryval;
	      delete coeff;
	    }
      
	}
    }


  cind = 0;//coarse index
  for(int j=0; j<Ny; j++)
    {
      for(int i=0; i<Nx; i++)
	{ 
	  int numbdrs = nbdrcon*subdomains[cind]->getnbdrs();	  
	  Array<double *> locnbdryval(numbdrs);
	  Array<double *> locrbdryval(numbdrs);
	  Array<double *> locrcoeff(numbdrs);
	  Array<double *> locdbdryval(numbdrs);
	  Array<bool> locBdrNeumann(numbdrs);
	  Array<bool> locBdrRobin(numbdrs);
	  Array<bool> locBdrDirichlet(numbdrs);
	  Array<bool> locGBdry(numbdrs);
	  
	  for(int u=0; u<nbdrcon; u++)
	    {
	      for(int b=0; b<subdomains[cind]->getnbdrs(); b++)
		{		  
		  int bb = b + u*subdomains[cind]->getnbdrs();
		  int numbv = subdomains[cind]->getnbdryelem(b)+1;
		  locnbdryval[bb] = new double[numbv];
		  locrbdryval[bb] = new double[numbv];
		  locrcoeff[bb] = new double[numbv];
		  locdbdryval[bb] = new double[numbv];
		  
		  if((j==0 && b==0) || (j==Ny-1 && b==2))
		    {
		      locBdrRobin[bb] = BdrRobin[bb];
		      locBdrNeumann[bb] = BdrNeumann[bb];
		      locBdrDirichlet[bb] = BdrDirichlet[bb];
		      locGBdry[bb] = true;
		      
		      for(int k=0; k<numbv; k++)
			{
			  int dd = i*nx+k;
			  locrcoeff[bb][k] = rcoeff[bb][dd];
			  locrbdryval[bb][k] = rbdryval[bb][dd];
			  locdbdryval[bb][k] = dbdryval[bb][dd];
			  locnbdryval[bb][k] = nbdryval[bb][dd];
			}
		    } 
		  else if((i==0 && b==3) || (i==Nx-1 && b==1))
		    {
		      locBdrRobin[bb] = BdrRobin[bb];
		      locBdrNeumann[bb] = BdrNeumann[bb];
		      locBdrDirichlet[bb] = BdrDirichlet[bb];
		      locGBdry[bb] = true;
		      
		      for(int k=0; k<numbv; k++)
			{
			  int dd = j*ny+k;
			  locrcoeff[bb][k] = rcoeff[bb][dd];
			  locrbdryval[bb][k] = rbdryval[bb][dd];
			  locdbdryval[bb][k] = dbdryval[bb][dd];
			  locnbdryval[bb][k] = nbdryval[bb][dd];
			}

		    }
		  else
		    {
		      locBdrRobin[bb] = true;
		      locBdrNeumann[bb] = false;
		      locBdrDirichlet[bb] = false;
		      locGBdry[bb] = false;
		  
		      for(int k=0; k<numbv; k++)
			{
			  int indie = mesh->GetBdryGindex(b, k);
			  double* coord = mesh->GetVertex(indie);

			  locrcoeff[bb][k] = gammafunc->Eval(coord);
			  locrbdryval[bb][k] = 0.0;
			  locdbdryval[bb][k] = 0.0;
			  locnbdryval[bb][k] = 0.0;
			}
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
  
  Data *data;
  data = new Data();
  data->SetNumberofBoundaryConditions(nbdrcon);
  data->SetEllipticData(ellipticdata);
  data->SetForceData(forcedata);
  data->SetNeumannData(nbdryval, BdrNeumann);
  data->SetRobinData(rcoeff, rbdryval, BdrRobin);
  data->SetDirichletData(dbdryval, BdrDirichlet);
  data->SetPointSourceData(nodevalue, nodeindex);
  cout << "Data Set" << endl;
  //---------------------------------------------------------
  // Compute Multiscale Basis Functions
  //--------------------------------------------------------- 
  
  Array<int> ns(nbdrcon*nbdry);
  input_file.getline( buffer, bufflen, ' ' );
  input_file >> ns[0] >> ns[1] >> ns[2] >> ns[3] >> ns[4] >> ns[5] >> ns[6] >> ns[7];

  cind=0;
  for(int j=0; j<Ny; j++)
    {
      for(int i=0; i<Nx; i++)
	{
	  subdomains[cind]->buildmsbasisfcts(ns);
	  cind++;
	}
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
	  int Gi = j*Nx + i;
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

  for(int k=0; k<Nsd; k++)
    {
      int numsegl = subdomains[k]->getnumseg(3+nbredges);
      int numsegr = subdomains[k]->getnumseg(1+nbredges);
      if(numsegl != 0)
	{
	  subdomains[k]->setgbindex(3+nbredges, gbindex);
	  gbindex += numsegl+1; 
	}
      if(numsegr != 0)
	{
	  subdomains[k]->setgbindex(1+nbredges, gbindex);
	  gbindex += numsegr+1; 
	}
    }

  for(int i=0; i<Nx; i++)
    {
      for(int j=0; j<Ny; j++)
	{
	  int Gi = j*Nx + i;
	  int numsegb = subdomains[Gi]->getnumseg(nbredges);
	  int numsegt = subdomains[Gi]->getnumseg(2+nbredges);
	  if(numsegb != 0)
	    {
	      subdomains[Gi]->setgbindex(nbredges, gbindex);
	      gbindex += numsegb+1; 
	    }
	  if(numsegt != 0)
	    {
	      subdomains[Gi]->setgbindex(nbredges+2, gbindex);
	      gbindex += numsegt+1; 
	    }
	}
    }

  //---------------------------------------------------------
  // Set subdomain neighbor info
  //--------------------------------------------------------- 

  for(int k=0; k<Nsd; k++)
    {
      int numbdrs = subdomains[k]->getnbdrs();
      Array<int*> nbinfo(2);  
      for(int j=0; j<2; j++){ nbinfo[j] = new int[nbdrcon*numbdrs]; }
      for(int b=0; b<nbdrcon*numbdrs; b++){ nbinfo[0][b] = -1; nbinfo[1][b] = -1; }
      
      for(int i=0; i<nbdrcon; i++)
	{
	  int count = i*numbdrs;
	  if(subdomains[k]->getGBdry(0+count)==false){ nbinfo[0][0+count] = k-Nx; nbinfo[1][0+count]=2+count; }
	  if(subdomains[k]->getGBdry(1+count)==false){ nbinfo[0][1+count] = k+1;  nbinfo[1][1+count]=3+count; }
	  if(subdomains[k]->getGBdry(2+count)==false){ nbinfo[0][2+count] = k+Nx; nbinfo[1][2+count]=0+count; }
	  if(subdomains[k]->getGBdry(3+count)==false){ nbinfo[0][3+count] = k-1;  nbinfo[1][3+count]=1+count; }
	}

      //for(int b=0; b<nbdrcon*numbdrs; b++){ cout << k << " " << nbinfo[0][b] << " " <<  nbinfo[1][b] << endl; }
      subdomains[k]->setnbinfo(nbinfo);

      for(int j=0; j<2; j++){ delete []nbinfo[j]; }  	  
    }


  //---------------------------------------------------------
  //  Read and create exact solution if known.
  //--------------------------------------------------------- 
  Function *exactsolfuncx; 
  Function *exactsolfuncy;
  bool solknown = false;
  ReadIdentifier(input_file, buffer, bufflen);
  if (!strcmp(buffer, "exactsolknown"))
    {
      ReadIdentifier(input_file, buffer, bufflen);
      exactsolfuncx = ReadFunction(input_file);
      ReadIdentifier(input_file, buffer, bufflen);
      exactsolfuncy = ReadFunction(input_file);
      solknown = true;
    }
  // Close the main input file.
  input_file.close();

  //---------------------------------------------------------
  //  Build MsRM
  //--------------------------------------------------------- 
  directMsRM msolver(subdomains, gbindex, nbdrcon);

  msolver.solve();
  Vector msol;
  msolver.buildsol(msol, Nx, Ny, nx, ny);
  msolver.paraviewprintmesh();
  msolver.paraviewprintdeformedmesh();

  //---------------------------------------------------------
  //  Build FEM 
  //--------------------------------------------------------- 
  FEM *fem;
  if(mtype==Element::QUADRILATERAL && nbdrcon==1){ fem = new BilinearFEM2D(mesh, data); }
  if(mtype==Element::TRIANGLE && nbdrcon==1){ fem = new LinearFEM2D(mesh, data); }
  if(mtype==Element::QUADRILATERAL && nbdrcon>1){ fem = new BilinearFEMELAST2D(mesh, data); }
  if(mtype==Element::TRIANGLE && nbdrcon>1){ fem = new LinearFEMELAST2D(mesh, data); }
  fem->Assemble();
  fem->FinalizeMatrix();
  cout << "FEM Built" << endl;

  //---------------------------------------------------------
  //  Solve
  //---------------------------------------------------------
  Vector sol(dof);
  fem->Solve(sol, "multifrontal");
  cout << "FEM Problem Solved" << endl;  
  ofstream output;
  output.open("avmsdefmesh.vtk"); 
  mesh->ParaviewPrintDeformedMesh(output, msol, mtype);  
  output.close(); 

  output.open("refdefmesh.vtk"); 
  mesh->ParaviewPrintDeformedMesh(output, sol, mtype);  
  output.close(); 

  msolver.paraviewprintellipticdata();

  //---------------------------------------------------------
  //  Error Estimates
  //---------------------------------------------------------
  if (solknown==true)
    {    
      cout << endl;
      cout << "MsRobin: " << endl;
      double error = msolver.computel2error(exactsolfuncx, exactsolfuncy);

      Vector u1(nnn);
      Vector u2(nnn); 
      cout << "Averaged MsRobin" << endl;
      for(int i=0; i<nnn; i++){u1(i) = fabs(msol(i)-sol(i)); u2(i) = fabs(msol(i+nnn)-sol(i+nnn));}
      cout << "L2error u1: " << fem->ComputeL2Error(exactsolfuncx, u1) << endl;
      cout << "L2error u2: " << fem->ComputeL2Error(exactsolfuncy, u2) << endl;
      cout << endl;

      cout << "Trad: " << endl;
      for(int i=0; i<nnn; i++){u1(i) = sol(i); u2(i) = sol(i+nnn);}
      cout << "L2error u1: " << fem->ComputeL2Error(exactsolfuncx, u1) << endl;
      cout << "L2error u2: " << fem->ComputeL2Error(exactsolfuncy, u2) << endl;

      delete exactsolfuncx;
      delete exactsolfuncy;
    }

  //---------------------------------------------------------
  //  Freedom!!!!
  //---------------------------------------------------------
  delete mesh;
  delete data;
  delete fem;

  for(int k=0; k<ellipticdata.Size(); k++)
    delete []ellipticdata[k];

  for(int k=0; k<forcedata.Size(); k++)
    delete []forcedata[k];
 
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

