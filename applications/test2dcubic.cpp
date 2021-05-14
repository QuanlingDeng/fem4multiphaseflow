#include <cmath>
#include "../fem/fem_header.h"
#include "../timedependent/timedependent_header.h"

using namespace std;

int main(int argc, const char * argv[])
{
  //-----------------------------
  // Input File
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

  //-----------------------------
  // Build the Mesh
  //-----------------------------
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
  if(meshtype==Element::TRIANGLE) { cout << "meshtype: Triangle" << endl; }
  if(meshtype==Element::QUADRILATERAL) { cout << "meshtype: Quadrilateral" << endl; }
  cout << endl;

  //build mesh
  Mesh *mesh;
  if(meshtype==Element::TRIANGLE) { mesh = new TMesh<2>(N, L, Element::TRIANGLE); }
  if(meshtype==Element::QUADRILATERAL) { mesh = new TMesh<2>(N, L, Element::QUADRILATERAL); }
  cout << "Mesh Built" << endl;

  //print mesh
  //mesh->ParaviewPrint(meshtype);
  int dof = mesh->GetNV() + 2 * mesh->GetNEdges() + mesh->GetNE();
  
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

  Array<int> ind;
  for(int i=0; i<mesh->GetNE(); i++)
    {
      mesh->GetElementVertices(i, ind);
      
      Array<int> edges;
      mesh->GetElementEdges(i, edges);
      if(i%2==0)
	{
	  ind.Append( mesh->GetNV() + 2*edges[0] );
	  ind.Append( mesh->GetNV() + 2*edges[0] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[1] );
	  ind.Append( mesh->GetNV() + 2*edges[1] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[2] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[2] );
	}
      else
	{
	  ind.Append( mesh->GetNV() + 2*edges[0] );
	  ind.Append( mesh->GetNV() + 2*edges[0] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[1] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[1] );
	  ind.Append( mesh->GetNV() + 2*edges[2] + 1);
	  ind.Append( mesh->GetNV() + 2*edges[2] );
	}
      ind.Append( mesh->GetNV() + 2 * mesh->GetNEdges() + i );
      
      Array<double> x(10);
      Array<double> y(10);
      for(int k=0; k<3; k++)
	{
	  x[k] = mesh->GetVertex(ind[k], 0);
	  y[k] = mesh->GetVertex(ind[k], 1);	  
	}

      x[3] = ( 2.0*x[0] + x[1] ) / 3.0;
      y[3] = ( 2.0*y[0] + y[1] ) / 3.0;
      
      x[4] = ( 2.0*x[1] + x[0] ) / 3.0;
      y[4] = ( 2.0*y[1] + y[0] ) / 3.0;
      
      x[5] = ( 2.0*x[1] + x[2] ) / 3.0;
      y[5] = ( 2.0*y[1] + y[2] ) / 3.0;
     
      x[6] = ( 2.0*x[2] + x[1] ) / 3.0;
      y[6] = ( 2.0*y[2] + y[1] ) / 3.0;
      
      x[7] = ( 2.0*x[2] + x[0] ) / 3.0;
      y[7] = ( 2.0*y[2] + y[0] ) / 3.0;
      
      x[8] = ( 2.0*x[0] + x[2] ) / 3.0;
      y[8] = ( 2.0*y[0] + y[2] ) / 3.0;
      
      x[9] = ( x[0] + x[1] + x[2] ) / 3.0;
      y[9] = ( y[0] + y[1] + y[2] ) / 3.0;
      
      for(int k=0; k<10; k++)
	{
	  double coord[2];
	  coord[0] = x[k];
	  coord[1] = y[k];
	  permdata[0][ind[k]] = permfunc->Eval(coord);
	  forcedata[0][ind[k]] = force->Eval(coord);	  
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
      nbdryval[b] = new double[3*nvb-2];
      dbdryval[b] = new double[3*nvb-2];
      
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
	  coord[0] = (2.0*coordl[0] + coordr[0]) / 3.0;
	  coord[1] = (2.0*coordl[1] + coordr[1]) / 3.0;
	  	  
	  if(BdrDirichlet[b]==true){ dbdryval[b][2*k+nvb] = boundaryfuncs[b]->Eval(coord); }
	  else if(BdrNeumann[b]==true){ nbdryval[b][2*k+nvb] = boundaryfuncs[b]->Eval(coord); }

	  coord[0] = (coordl[0] + 2.0*coordr[0]) / 3.0;
	  coord[1] = (coordl[1] + 2.0*coordr[1]) / 3.0;
	  	  
	  if(BdrDirichlet[b]==true){ dbdryval[b][2*k+1+nvb] = boundaryfuncs[b]->Eval(coord); }
	  else if(BdrNeumann[b]==true){ nbdryval[b][2*k+1+nvb] = boundaryfuncs[b]->Eval(coord); }
	} 
    }

  //set data
  Data *data;
  data = new Data();
  //data->SetEllipticData(permdata);
  data->SetEllipticFunction(permfunc);
  //data->SetForceData(forcedata);
  data->SetForceFunction(force);
  data->SetDirichletData(dbdryval, BdrDirichlet);
  data->SetNeumannData(nbdryval, BdrNeumann);
  cout << "Data Set" << endl;

  //---------------------------------------------------------
  //  Set constructor and solve the problem
  //---------------------------------------------------------
  FEM *fem;
  if(meshtype==Element::QUADRILATERAL)
    { cout << "Not implemented" << endl; exit(1); }
  if(meshtype==Element::TRIANGLE){ fem = new CubicFEM2D(mesh, data); }
  fem->Assemble(); 
  fem->FinalizeMatrix();

  Vector sol(dof);
  sol = 0.0; 
  //fem->Solve(sol, "pcg", 10000, 0, 1.0e-20, 1.0e-40);
  fem->Solve(sol, "multifrontal", 10000, 1, 1.0e-12, 1.0e-24);
  //sol.Print();
  
  //if(meshtype==Element::TRIANGLE){ mesh->ParaviewPrint(sol, Element::TRIANGLE); }

  //---------------------------------------------------------
  //  Set Exact Solution
  //---------------------------------------------------------
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solu = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solux = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *soluy = ReadFunction(input_file);

  Array<Function *> dersol(2);
  dersol[0] = solux;
  dersol[1] = soluy;
  
  //cout << "L2 error: " << fem->ComputeL2Error(solu, sol) << endl;
  //cout << "The H1 SemiNorm error: " << fem->ComputeH1Error(dersol, sol) << endl;

  //---------------------------------------------------------
  //  Post-Proccessing
  //---------------------------------------------------------
  cout << "Post-processing..." << endl;
  DualMesh *dualmesh = new DualMesh(mesh, 3);

  Postprocessing *pp = new Postprocessing(dualmesh, fem, data);

  Array<double *> flux(18);
  Array<double *> ppsol(10);
  Array<double *> ppflux(18);
  for(int i=0; i<flux.Size(); i++) flux[i] = new double[mesh->GetNE()];
  for(int i=0; i<ppsol.Size(); i++) ppsol[i] = new double[mesh->GetNE()];
  for(int i=0; i<ppflux.Size(); i++) ppflux[i] = new double[mesh->GetNE()];

  Array<double> ppbdrflux(dof);
  
  pp->ComputeConservativeFlux(sol, flux, ppsol, ppflux, ppbdrflux);

  Array<double> l2errs(3);
  Array<double> h1errs(3);
  fem->ComputePPL2Error(l2errs, solu, sol, ppsol);
  fem->ComputePPH1Error(h1errs, dersol, sol, ppsol);

  cout << "femL2err: " << l2errs[0] << endl;
  cout << "femH1err: " << h1errs[0] << endl;

  cout << endl;
  cout << "ppL2err: " << l2errs[1] << endl;
  cout << "ppH1err: " << h1errs[1] << endl;

  cout << endl;
  cout << "diffL2err: " << l2errs[2] << endl;
  cout << "diffH1err: " << h1errs[2] << endl;
    
  cout << "Post-processing done" << endl;
  
  //------ Memory Free ----------------
  delete pp;
  delete dualmesh;
  for(int i=0; i<flux.Size(); i++) delete []flux[i];
  for(int i=0; i<ppflux.Size(); i++) delete []ppflux[i];
  for(int i=0; i<ppsol.Size(); i++) delete []ppsol[i];
  delete permfunc;
  delete force;
  delete solu;
  delete solux;
  delete soluy;
  delete data;
  delete mesh;
  delete []permdata[0];
  delete []forcedata[0];
  for(int k=0; k<dbdryval.Size(); k++) delete []dbdryval[k];
  for(int k=0; k<nbdryval.Size(); k++) delete []nbdryval[k];
  for (int k=0; k<boundaryfuncs.Size(); k++) delete boundaryfuncs[k]; 

  return 0;
}
