#include <iostream>
#include <fstream>
#include <cmath>
#include "../linalg/linalg_header.h"
#include "../mesh/mesh_header.h"
#include "../linalg/linalg_header.h"
#include "../fem/fem_header.h"

using namespace std;


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


 cout << N[0] << " " << N[1] << " " << meshtype << endl;

  double hx = 1.0/N[0];
  double hy = 1.0/N[1];

  cout << meshtype << " " << Element::TRIANGLE << " " << Element::QUADRILATERAL << endl;
  
  //build
  Mesh *mesh;
  if(meshtype==Element::TRIANGLE){mesh = new TMesh<2>(N, L, Element::TRIANGLE);}
  if(meshtype==Element::QUADRILATERAL){ mesh = new TMesh<2>(N, L, Element::QUADRILATERAL); }

  cout << endl << "Mesh Constructed" << endl;

  //print
  ofstream out;
  out.open("mesh.out");
  mesh->Print(out, meshtype);
  out.close();
  mesh->ParaviewPrint(meshtype);
   
  //---------------------------------------------------------
  // Build Data
  //--------------------------------------------------------- 
  int nbdrcon = 2;
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
      
      for(int k=0; k<nnn; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  double data = modfunc->Eval(coord);
	  double prat = 0.3;
	  ellipticdata[0][k] = data*prat/((1.0+prat)*(1.0-2.0*prat));
	  ellipticdata[1][k] = data/(1+prat);
	}  
    
      delete modfunc;
    }
  else if (!strcmp(elmodflag, "lame"))
    {
      Function *lambda;
      Function *mu;
      ReadIdentifier(input_file, buffer, bufflen);
      lambda = ReadFunction(input_file);
      ReadIdentifier(input_file, buffer, bufflen);
      mu = ReadFunction(input_file);
      for(int k=0; k<nnn; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  ellipticdata[0][k] = lambda->Eval(coord);
	  ellipticdata[1][k] = mu->Eval(coord);
	}  
    
      delete lambda;
      delete mu;
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
    }
 
  cout << "end" << endl;
  Array<double *> forcedata(nbdrcon);
  forcedata[0] = new double[nnn];
  forcedata[1] = new double[nnn];

  for(int i=0; i<forcedata.Size(); i++)
    {
      ReadIdentifier(input_file, buffer, bufflen);
      Function *force = ReadFunction(input_file);
      for(int k=0; k<nnn; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  forcedata[i][k] = force->Eval(coord);
	  //cout <<  force->Eval(coord) << endl;
	}
      delete force;
    }


  //Gather Boundary Data
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
	  int bb = b + i*nbredges;
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

  Data *data;
  data = new Data();
  data->SetNumberofBoundaryConditions(nbdrcon);
  data->SetEllipticData(ellipticdata);
  data->SetForceData(forcedata);
  data->SetNeumannData(nbdryval, BdrNeumann);
  data->SetRobinData(rcoeff, rbdryval, BdrRobin);
  data->SetDirichletData(dbdryval, BdrDirichlet);
  cout << "Data Set" << endl;
  
  double deltaforce;
  input_file.getline( buffer, bufflen, ' ' );
  input_file >> deltaforce;

  Array<double> nodevalue(1);
  Array<int> nodeindex(1);
  
  nodevalue[0] = deltaforce;
  nodeindex[0] = nnn+N[1]*(N[0]+1)/2+N[0]/2;
  data->SetPointSourceData(nodevalue, nodeindex);
  
  //---------------------------------------------------------
  //  Read and create exact solution if known.
  //--------------------------------------------------------- 
  /*  
Function *exactsolfuncx; 
  Function *exactsolfuncy;
  ReadIdentifier(input_file, buffer, bufflen);

  if (!strcmp(buffer, "exactsolknown"))
    {
      ReadIdentifier(input_file, buffer, bufflen);
      exactsolfuncx = ReadFunction(input_file);
      ReadIdentifier(input_file, buffer, bufflen);
      exactsolfuncy = ReadFunction(input_file);
    }
*/
  // Close the main input file.
  input_file.close();

  //---------------------------------------------------------
  //  Build FEM 
  //---------------------------------------------------------
  //BilinearFEMELAST2D *fem = new BilinearFEMELAST2D(mesh, data);
  //LinearFEMELAST2D *fem = new LinearFEMELAST2D(mesh, data);

     
  FEM *fem;
  if(meshtype==Element::QUADRILATERAL){ fem = new BilinearFEMELAST2D(mesh, data); }
  if(meshtype==Element::TRIANGLE){ fem = new LinearFEMELAST2D(mesh, data); }
    
  fem->Assemble();
  fem->FinalizeMatrix();
  cout << "FEM Built" << endl;


  //---------------------------------------------------------
  //  Solve
  //---------------------------------------------------------
  Vector sol(dof);
  sol = 0;

  cout << "solving..." << endl;
  fem->Solve(sol, "multifrontal");
  cout << "solved" << endl;
  
  
  Vector u1(nnn);
  Vector u2(nnn);
  for(int i=0; i<nnn; i++){u1(i) = sol(i); u2(i) = sol(i+nnn);}
  //cout << "L2error u1: " << fem->ComputeL2Error(exactsolfuncx, u1) << endl;
  //cout << "L2error u2: " << fem->ComputeL2Error(exactsolfuncy, u2) << endl;
  
  //---------------------------------------------------------
  //  Print
  //---------------------------------------------------------
  ofstream output;
  output.open("defmesh.vtk"); 
  mesh->ParaviewPrintDeformedMesh(output, sol, meshtype);  
  output.close(); 

  output.open("u1.vtk"); 
  mesh->ParaviewPrint(output, u1, meshtype);  
  output.close(); 

  output.open("u2.vtk"); 
  mesh->ParaviewPrint(output, u2, meshtype);  
  output.close(); 
  //---------------------------------------------------------
  //  Freedom!!!!
  //---------------------------------------------------------
  delete mesh;
  delete data;
  delete fem;
  // delete exactsolfuncx;
  //delete exactsolfuncy;

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

  cout << endl << "fin" << endl;

}


 
