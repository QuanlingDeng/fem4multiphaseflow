#include <iostream>
#include <fstream>
#include <cmath>
#include "../linalg/linalg_header.h"
#include "../mesh/mesh_header.h"
#include "../linalg/linalg_header.h"
#include "../fem/fem_header.h"

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

  Array<int> N(2);
  int meshtype;
  input_file.getline( buffer, bufflen, ' ' );
  input_file >> N[0] >> N[1] >> meshtype;

  double hx = 1.0/N[0];
  double hy = 1.0/N[1];
  cout <<  Element::TRIANGLE << endl;
  //build
  Mesh *mesh;
  if(meshtype==Element::TRIANGLE){mesh = new TMesh<2>(N, L, Element::TRIANGLE); cout << "meshbuilding" << endl; }
  if(meshtype==Element::QUADRILATERAL){ mesh = new TMesh<2>(N, L, Element::QUADRILATERAL); }
  /*
  cout << endl << "Mesh Constructed" << endl << endl;
  DualMesh *dual = new DualMesh(mesh);

  delete dual;
  */
  //print
  ofstream out;
  out.open("mesh.out");
  mesh->Print(out, meshtype);
  out.close();
  mesh->ParaviewPrint(meshtype);
   
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

  Data *data;
  data = new Data();
  data->SetEllipticData(ellipticdata);
  data->SetForceData(forcedata);
  data->SetNeumannData(nbdryval, BdrNeumann);
  data->SetRobinData(rcoeff, rbdryval, BdrRobin);
  data->SetDirichletData(dbdryval, BdrDirichlet);
  /*
  Array<double> nodevalue(N[1]+1);
  Array<int> nodeindex(N[1]+1);
  for(int i=0; i<=N[1]; i++)
    {
      nodevalue[i] = 1.0;
      nodeindex[i] =  i*(N[0]+1)+N[0]/2;
    }
  */
  Array<double> nodevalue(1);
  Array<int> nodeindex(1);
  nodevalue[0] = 1.0;
  nodeindex[0] = N[1]*(N[0]+1)/2+N[0]/2;
  

  data->SetPointSourceData(nodevalue, nodeindex);
  cout << "Data Set" << endl;
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
  //  Build FEM 
  //--------------------------------------------------------- 
  FEM *fem;
  if(meshtype==Element::QUADRILATERAL){ fem = new BilinearFEM2D(mesh, data); }
  if(meshtype==Element::TRIANGLE){ fem = new LinearFEM2D(mesh, data); }
  fem->Assemble();
  fem->FinalizeMatrix();
  cout << "FEM Built" << endl;

  //---------------------------------------------------------
  //  Solve
  //---------------------------------------------------------
  Vector sol(nnn);
  fem->Solve(sol, "multifrontal");
  cout << "FEM Problem Solved" << endl;  

  DualMesh *dual = new DualMesh(mesh);
  dual->PrintDualMesh();

  //---------------------------------------------------------
  //  Print
  //---------------------------------------------------------
  /* 
 out.open("sol.out");
  mesh->Print(out, sol);
  out.close();
  
  ofstream output;
  output.open("sol.vtk"); 
  mesh->ParaviewPrintDeformedMesh(output, sol, meshtype);  
  output.close(); 

  mesh->ParaviewPrint(sol, meshtype);

  cout << "L2error: " << fem->ComputeL2Error(exactsolfunc, sol) << endl;
  cout << "H1error: " << fem->ComputeH1Error(exactderfunc, sol) << endl;
  */
  //---------------------------------------------------------
  //  Freedom!!!!
  //---------------------------------------------------------
  delete mesh;
  delete []ellipticdata[0];
  delete []forcedata[0];
  delete data;
  delete fem;
  delete dual;
  //delete exactsolfunc;

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

  cout << endl << "fin" << endl;
}

