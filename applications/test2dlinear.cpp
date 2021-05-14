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
  mesh->ParaviewPrint(meshtype);
  int NV = mesh->GetNV();

  cout << NV << endl;

  //---------------------------------------------------------
  // Gather Data for Diffusion Equation
  //---------------------------------------------------------
  //Set Permeability Data
  Array<double *> permdata(1);
  permdata[0] = new double[NV];

  char permflag[bufflen]; 
  ReadIdentifier(input_file, buffer, bufflen);
  input_file >> permflag;
  cout << permflag << endl;
  
  Function *permfunc;
  if (!strcmp(permflag, "equation"))
    {
      ReadIdentifier(input_file, buffer, bufflen);
      permfunc = ReadFunction(input_file);

      for(int k=0; k<NV; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  permdata[0][k] = permfunc->Eval(coord);
	}      
    }

  //Gather Force Data
  ReadIdentifier(input_file, buffer, bufflen);
  Function *force = ReadFunction(input_file);
  Array<double *> forcedata(1);
  forcedata[0] = new double[NV];

  for(int k=0; k<NV; k++)
    {
      double* coord = mesh->GetVertex(k);
      forcedata[0][k] = force->Eval(coord);
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
      nbdryval[b] = new double[nvb];
      dbdryval[b] = new double[nvb];
      
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

  //---------------------------------------------------------
  //  Set constructor and solve the problem
  //---------------------------------------------------------
  FEM *fem;
  if(meshtype==Element::QUADRILATERAL){ fem = new BilinearFEM2D(mesh, data); }
  if(meshtype==Element::TRIANGLE){ fem = new LinearFEM2D(mesh, data); }
  fem->Assemble(); 
  fem->FinalizeMatrix();

  Vector sol(NV);
  sol = 0.0;
  fem->Solve(sol, "multifrontal", 10000, 0, 1.0e-20, 1.0e-40);
  //sol.Print();
  if(meshtype==Element::QUADRILATERAL){ mesh->ParaviewPrint(sol, Element::QUADRILATERAL); }
  if(meshtype==Element::TRIANGLE){ mesh->ParaviewPrint(sol, Element::TRIANGLE); }

  //---------------------------------------------------------
  //  Set Exact Solution
  //---------------------------------------------------------
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solu = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solux = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *soluy = ReadFunction(input_file);

  Vector truesol(NV);
  for(int k=0; k<NV; k++)
    {
      double* coord = mesh->GetVertex(k);
      truesol(k) = solu->Eval(coord);
    }

  Array<Function *> dersol(2);
  dersol[0] = solux;
  dersol[1] = soluy;
  
  //cout << "L2 error: " << fem->ComputeL2Error(solu, sol) << endl;
  //cout << "The H1 SemiNorm error: " << fem->ComputeH1Error(dersol, sol) << endl;

  //---------------------------------------------------------
  //  Post-Proccessing
  //---------------------------------------------------------
  cout << "Post-processing..." << endl;
  DualMesh *dualmesh = new DualMesh(mesh);

  Postprocessing *pp = new Postprocessing(dualmesh, fem, data);

  int nl = (meshtype==Element::TRIANGLE) ? 3 : 4;
  Array<double *> flux(nl);
  Array<double *> ppsol(nl);
  Array<double *> ppflux(nl);
  for(int i=0; i<nl; i++) flux[i] = new double[mesh->GetNE()];
  for(int i=0; i<nl; i++) ppsol[i] = new double[mesh->GetNE()];
  for(int i=0; i<nl; i++) ppflux[i] = new double[mesh->GetNE()];

  Array<double> ppbdrflux(mesh->GetNV());
  
  pp->ComputeConservativeFlux(truesol, sol, flux, ppsol, ppflux, ppbdrflux);

  Array<double> l2errs(3);
  Array<double> h1errs(3);
  fem->ComputePPL2Error(l2errs, solu, sol, ppsol);
  fem->ComputePPH1Error(h1errs, dersol, sol, ppsol);

  cout << endl;
  //cout << "femL2err: " << l2errs[0] << endl;
  cout << "femH1err: " << h1errs[0] << endl;

  //cout << endl;
  //cout << "ppL2err: " << l2errs[1] << endl;
  cout << "ppH1err: " << h1errs[1] << endl;

  //cout << endl;
  //cout << "diffL2err: " << l2errs[2] << endl;
  //cout << "diffH1err: " << h1errs[2] << endl << endl;
    
  cout << "Post-processing done" << endl;
  
  //------ Memory Free ----------------
  delete pp;
  delete dualmesh;
  for(int i=0; i<flux.Size(); i++) delete []flux[i];
  for(int i=0; i<ppflux.Size(); i++) delete []ppflux[i];
  for(int i=0; i<ppsol.Size(); i++) delete []ppsol[i];
  delete force;
  delete permfunc;
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
