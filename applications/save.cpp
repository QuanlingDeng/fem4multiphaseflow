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

  //--------------------------------------------------------------------------------
  // Set time data.
  //--------------------------------------------------------------------------------
  double finaltime;
  int numtimesteps;
  
  input_file.getline( buffer, bufflen, ' ' );
  input_file >> finaltime >> numtimesteps;
  double dt = finaltime/(numtimesteps);

  cout << endl << "Time Discretiazations: " << endl;
  cout << "Final Time: " << finaltime << endl;
  cout << "Number of Time Steps: " << numtimesteps << endl;
  cout << "Time Step Length: " << dt << endl; 
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

  double hx = (L[1]-L[0])/N[0];
  double hy = (L[3]-L[2])/N[1];
  double h = sqrt(hx*hx + hy*hy); 

  //print mesh
  mesh->ParaviewPrint(meshtype);
  int nv = mesh->GetNV(); //number of vertices

  //---------------------------------------------------------
  // Gather Data for Hydrodynamic Pressure Equation
  //---------------------------------------------------------
  //Set Permeability Data
  Array<double *> permdata(1);
  permdata[0] = new double[nv];

  char permflag[bufflen]; 
  ReadIdentifier(input_file, buffer, bufflen);
  input_file >> permflag;
  cout << permflag << endl;
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
  ReadIdentifier(input_file, buffer, bufflen);
  Function *force = ReadFunction(input_file);

  //Gather Boundary Data
  int nbdry = mesh->GetNBdrs();
  Array<double *> nbdryval(nbdry);
  Array<double *> dbdryval(nbdry);
  Array<bool> BdrNeumann(nbdry);
  Array<bool> BdrDirichlet(nbdry);
  Array<Function *> pressureboundaryfuncs(nbdry);

  char bdryflag[bufflen]; 
  for(int b=0; b<nbdry; b++)
    {
      int nvb = mesh->GetNBdryElem(b)+1;
      nbdryval[b] = new double[nvb];
      dbdryval[b] = new double[nvb];
      
      input_file.getline( buffer, bufflen, ' ' );
      input_file >> bdryflag;
      ReadIdentifier(input_file, buffer, bufflen);
      pressureboundaryfuncs[b] = ReadFunction(input_file);

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
	  nbdryval[b][k] = 0.0;
	  dbdryval[b][k] = 0.0;
	}  
    }
  
  //set data for pressure equation
  Data *pressuredata;
  pressuredata = new Data();
  pressuredata->SetEllipticData(permdata);
  cout << "Pressure Data Set" << endl;


  //---------------------------------------------------------
  //  Gather Data for Poroelasticity Problem
  //--------------------------------------------------------- 

  int nnn = 2*nv; //number of dof
  //Gather elasticitiy Data
  Array<double *> lamecoeff(2);
  lamecoeff[0] = new double[nv];
  lamecoeff[1] = new double[nv];
  Vector elmod(nv);

  char elmodflag[bufflen]; 
  ReadIdentifier(input_file, buffer, bufflen);
  input_file >> elmodflag;
  double prat = 0.25;

  if (!strcmp(elmodflag, "equation"))
    {
      Function *modfunc;
      ReadIdentifier(input_file, buffer, bufflen);
      modfunc = ReadFunction(input_file);
      
      for(int k=0; k<nv; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  double data = modfunc->Eval(coord);
	  elmod(k) = log(data);
	  lamecoeff[0][k] = data*prat/((1.0+prat)*(1.0-2.0*prat));
	  lamecoeff[1][k] = data/(1+prat);
	}
      delete modfunc;
    }
  if(!strcmp(elmodflag, "data"))
    {
      char modfile[bufflen];
      double sigma;
      double a;
      ifstream modin;
      input_file >> modfile >> sigma >> a;
      modin.open(modfile);
      if(!modin) { cout << modfile << " does not exist.\n";  exit(1); }

      double val;
      double data;
      double datamin = a*1000.0;
      double datamax = 0.0;
      for (int k=0; k<nv; k++)
	{
	  modin >> val;
	  data = a*exp(sigma * val);
	  elmod(k) = log(data);
	  lamecoeff[0][k] = data*prat/((1.0+prat)*(1.0-2.0*prat));
	  lamecoeff[1][k] = data/(1+prat);
	  if(data > datamax){datamax = data; }
	  if(data < datamin){datamin = data; }
	}
      cout << "E Range: [" << datamin << ", " << datamax << "]" << endl;
      modin.close();
    }


  char filename[256];
  ofstream output;
  sprintf(filename, "elmod.vtk");
  output.open(filename);
  mesh->ParaviewPrint(output, elmod);
  output.close();

  cout << elmodflag << endl;
  ReadIdentifier(input_file, buffer, bufflen);
  Function *forcex = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *forcey = ReadFunction(input_file);

  //Set Boundary Data
  int nbredges = mesh->GetNBdrs();
  int nbdryel = 2*nbredges;
  Array<double *> elnbdryval(nbdryel);
  Array<double *> eldbdryval(nbdryel);
  Array<bool> elBdrNeumann(nbdryel);
  Array<bool> elBdrDirichlet(nbdryel);
  Array<Function *> elasticboundaryfuncs(nbdryel);

  for(int b=0; b<nbredges; b++)
    {
      int nvb = mesh->GetNBdryElem(b)+1;
      for(int i=0; i<2; i++)
	{
	  int bb = b + i*nbredges;
	  elnbdryval[bb] = new double[nvb];
	  eldbdryval[bb] = new double[nvb];
	  
	  input_file.getline( buffer, bufflen, ' ' );
	  input_file >> bdryflag;
	  ReadIdentifier(input_file, buffer, bufflen);
	  elasticboundaryfuncs[bb] = ReadFunction(input_file);
	  cout << bdryflag << endl;
	  if (!strcmp(bdryflag, "neumann"))
	    {
	      elBdrNeumann[bb] = true;
	      elBdrDirichlet[bb] = false; 
	    }
	  if (!strcmp(bdryflag, "dirichlet"))
	    {
	      elBdrDirichlet[bb] = true;
	      elBdrNeumann[bb] = false;
	    }
	  
	  for(int k=0; k<nvb; k++)
	    {
	      elnbdryval[bb][k] = 0.0;
	      eldbdryval[bb][k] = 0.0;
	    }  
	}
    }

  //set elasticity data
  Data *elasticdata;
  elasticdata = new Data();
  elasticdata->SetNumberofBoundaryConditions(2);
  elasticdata->SetEllipticData(lamecoeff);

  cout << "Elastic Data Set" << endl;   
  //---------------------------------------------------------
  //  Set Exact Solutions
  //---------------------------------------------------------
  
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solp = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solu1 = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solu2 = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solpx = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solpy = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solu1t = ReadFunction(input_file);
  ReadIdentifier(input_file, buffer, bufflen);
  Function *solu2t = ReadFunction(input_file);

  cout <<buffer << endl;
  //---------------------------------------------------------
  //  Time March 
  //---------------------------------------------------------
  Vector u1(nv);
  Vector u1old(nv);
  Vector u2(nv);
  Vector u2old(nv);
  Vector pressure(nv);
  Vector pressureold(nv);
  Vector displacement(2*nv);
  Vector displacementold(2*nv);

  Vector sol(3*nv);
  Vector solold(3*nv);
  u1old = 0.0; 
  u2old = 0.0;
  pressureold = 0.0;
  displacementold = 0.0;
  sol = 0.0;
  solold = 0.0;

  cout << endl;
  double currtime = 0.0;
  for (int n=1; n<=numtimesteps; n++) 
    { 
      DualMesh *dualmesh = new DualMesh(mesh);
      cout << "//********************************************" << endl;
      cout << "Current Time: " << endl << currtime/(24.0*60.0*60.0) << " Days, " << 100*currtime/finaltime << "%" << endl; 

      //---------------------------------------------------------
      //  Print
      //---------------------------------------------------------
      cout << "Printing..." << flush;  
      sprintf(filename, "pressure_%f.vtk", currtime);
      output.open(filename);
      mesh->ParaviewPrint(output, pressureold); 
      output.close();  

      sprintf(filename, "defmesh_%f.vtk", currtime);
      output.open(filename);
      mesh->ParaviewPrintDeformedMesh(output, displacementold);
      output.close();  

      currtime += dt;
      //---------------------------------------------------------
      //  Solve Hydrodynamics and Poroelastic System
      //---------------------------------------------------------
      SparseMatrix *globalmat;
      globalmat = new SparseMatrix(3*nv, 50);
      Vector globalrhs(3*nv);
      globalrhs = 0.0;

      cout << endl << "Building Global System:" << endl;
      cout << "Elasticity Block..." << flush;
      //Fill Ae to Global Matrix
      Array<double *> elfdata(2);
      elfdata[0] = new double[nv];
      elfdata[1] = new double[nv];

      for(int k=0; k<nv; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  double coords[3];
	  coords[0] = coord[0];
	  coords[1] = coord[1];
	  coords[2] = currtime;
	  elfdata[0][k] = forcex->Eval(coords);
	  elfdata[1][k] = forcey->Eval(coords);
	}
      elasticdata->SetForceData(elfdata);

      for(int b=0; b<mesh->GetNBdrs(); b++)
	{
	  int nvb = mesh->GetNBdryElem(b)+1;
	  for(int i=0; i<2; i++)
	    {
	      int bb = b + i*mesh->GetNBdrs();
	      for(int k=0; k<nvb; k++)
		{
		  int indie = mesh->GetBdryGindex(b, k);
		  double* coord = mesh->GetVertex(indie);
		  double coords[3];
		  coords[0] = coord[0];
		  coords[1] = coord[1];
		  coords[2] = currtime;
		  if(elBdrDirichlet[bb]==true){ eldbdryval[bb][k] = elasticboundaryfuncs[bb]->Eval(coords); }
		  else if(elBdrNeumann[bb]==true){ elnbdryval[bb][k] = elasticboundaryfuncs[bb]->Eval(coords); }
		}
	    }
	}
      elasticdata->SetDirichletData(eldbdryval, elBdrDirichlet);
      elasticdata->SetNeumannData(elnbdryval, elBdrNeumann);

      FEM *elasticfem;
      if(meshtype==Element::QUADRILATERAL){ elasticfem = new BilinearFEMELAST2D(mesh, elasticdata); }
      if(meshtype==Element::TRIANGLE){ elasticfem = new LinearFEMELAST2D(mesh, elasticdata); }
      elasticfem->AssembleMatrix();
      elasticfem->FinalizeMatrix();
      SparseMatrix elastmatrix = elasticfem->getMatrix();
      
      int *I = elastmatrix.GetI();
      int *J = elastmatrix.GetJ();
      for(int k=0; k<2*nv; k++)
	{
	  for(int i=I[k]; i<I[k+1]; i++)
	    {
	      globalmat->Elem(k, J[i]) = elastmatrix(k, J[i]);
	    }
	}

      Vector urhs(2*nv);
      elasticfem->getRHS(urhs);      
      for(int j=0; j<2*nv; j++)
	{
	  globalrhs(j) += urhs(j);
	}
      cout << "done" << endl;

      //Fill Ap to Global Matrix  
      cout << "Pressure Block..." << flush;
      Array<double *> forcedata(1);
      forcedata[0] = new double[nv];

      for(int k=0; k<nv; k++)
	{
	  double* coord = mesh->GetVertex(k);
	  double coords[3];
	  coords[0] = coord[0];
	  coords[1] = coord[1];
	  coords[2] = currtime;
	  forcedata[0][k] = force->Eval(coords);
	}
      pressuredata->SetForceData(forcedata);

      for(int b=0; b<nbdry; b++)
	{
	  int nvb = mesh->GetNBdryElem(b)+1;
	  for(int k=0; k<nvb; k++)
	    {
	      int indie = mesh->GetBdryGindex(b, k);
	      double* coord = mesh->GetVertex(indie);
	      double coords[3];
	      coords[0] = coord[0];
	      coords[1] = coord[1];
	      coords[2] = currtime;
	      if(BdrDirichlet[b]==true){ dbdryval[b][k] = pressureboundaryfuncs[b]->Eval(coords); }
	      else if(BdrNeumann[b]==true){ nbdryval[b][k] = pressureboundaryfuncs[b]->Eval(coords); }
	    }
	}     
      pressuredata->SetDirichletData(dbdryval, BdrDirichlet);
      pressuredata->SetNeumannData(nbdryval, BdrNeumann);
      
      FEM *pressurefem;
      if(meshtype==Element::QUADRILATERAL){ pressurefem = new BilinearFEM2D(mesh, pressuredata); }
      if(meshtype==Element::TRIANGLE){ pressurefem = new LinearFEM2D(mesh, pressuredata); }
      pressurefem->AssembleMatrix(); 
      pressurefem->FinalizeMatrix();
      SparseMatrix pressmatrix = pressurefem->getMatrix();
      
      I = pressmatrix.GetI();
      J = pressmatrix.GetJ();
      for(int k=0; k<nv; k++)
	{
	  for(int i=I[k]; i<I[k+1]; i++)
	    {
	      globalmat->Elem(k+2*nv, J[i]+2*nv) = pressmatrix(k, J[i]);
	    }
	}

      Vector prhs(nv);
      pressurefem->getRHS(prhs);      
      for(int j=0; j<nv; j++)
	{
	  globalrhs(j+2*nv) += prhs(j);
	}
      
      cout << "done" << endl;
      cout << "Stabilization Block..." << flush;
      Array<double *> corrector(1);
      corrector[0] = new double[nv];
      for(int i=0; i<nv; i++)
	{
	  corrector[0][i] = 0.25*h*h*(exp(elmod(i)))/(double(dt));
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

      //Add Corrector to RHS 
      Vector betap0(nv);
      betamatrix.Mult(pressureold, betap0);
      for(int k=0; k<nv; k++)
	{
	  globalrhs(k+2*nv) += betap0(k);
	}

      cout << "done" << endl;
      cout << "Off Diagonal Blocks..." << flush;      
      //Set off diagonal data
      Data *scalardatax;
      scalardatax = new Data();
      
      Array<double *> scalarcoeffx(2);
      scalarcoeffx[0] = new double[nv];
      scalarcoeffx[1] = new double[nv];
      for(int i=0; i<nv; i++)
	{
	  scalarcoeffx[0][i] = 1.0; 
	  scalarcoeffx[1][i] = 0.0; 
	}
      scalardatax->SetAdvectionData(scalarcoeffx);
      
      FEM *scalarfemx;
      if(meshtype==Element::QUADRILATERAL){ scalarfemx = new BilinearFEM2D(mesh, scalardatax); }
      if(meshtype==Element::TRIANGLE){ scalarfemx = new LinearFEM2D(mesh, scalardatax); }
      scalarfemx->AssembleMatrix();
      scalarfemx->FinalizeMatrix(); 
      SparseMatrix Q = scalarfemx->getMatrix();

      I = Q.GetI();
      J = Q.GetJ();
      for(int k=0; k<nv; k++)
	{
	  for(int i=I[k]; i<I[k+1]; i++)
	    {
	      globalmat->Elem(k, J[i]+2*nv) = Q(k, J[i]);
	      globalmat->Elem(k+2*nv, J[i]) = Q(k, J[i])/dt;
	    }
	}
      Data *scalardatay;
      scalardatay = new Data();
      
      Array<double *> scalarcoeffy(2);
      scalarcoeffy[0] = new double[nv];
      scalarcoeffy[1] = new double[nv];
      for(int i=0; i<nv; i++)
	{
	  scalarcoeffy[0][i] = 0.0; 
	  scalarcoeffy[1][i] = 1.0; 
	}
      scalardatay->SetAdvectionData(scalarcoeffy);

      FEM *scalarfemy;
      if(meshtype==Element::QUADRILATERAL){ scalarfemy = new BilinearFEM2D(mesh, scalardatay); }
      if(meshtype==Element::TRIANGLE){ scalarfemy = new LinearFEM2D(mesh, scalardatay); }
      scalarfemy->AssembleMatrix(); 
      scalarfemy->FinalizeMatrix(); 
      SparseMatrix R = scalarfemy->getMatrix();

      I = R.GetI();
      J = R.GetJ();
      for(int k=0; k<nv; k++)
	{
	  for(int i=I[k]; i<I[k+1]; i++)
	    {
	      globalmat->Elem(k+nv, J[i]+2*nv) = R(k, J[i]);
	      globalmat->Elem(k+2*nv, J[i]+nv) = R(k, J[i])/dt;
	    }
	}


      Vector Qu1(nv);
      Vector Ru2(nv);
      Q.Mult(u1old, Qu1);
      R.Mult(u2old, Ru2);
      
      for(int j=0; j<nv; j++)
	{
	  globalrhs(j+2*nv) += (Qu1(j)+Ru2(j))/dt;
	}

      cout << "done" << endl;
      cout << "Boundary Conditions..." << flush;  
  
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
      cout << "done" << endl << endl;
      cout << "Solving..." << flush;  

      globalmat->Finalize();
      int mem_fact = 100;
      MultiFrontalMatrixInverse *InvMat = new MultiFrontalMatrixInverse(*globalmat, mem_fact);
      InvMat->Mult(globalrhs, sol);
      delete InvMat;
      cout << "done" << endl;

      for(int i=0; i<nv; i++){u1(i) = sol(i); u2(i) = sol(i+nv); pressure(i) = sol(i+2*nv);}
      for(int i=0; i<2*nv; i++){displacement(i) = sol(i); }
      
      //---------------------------------------------------------
      //  Post-Proccessing
      //---------------------------------------------------------
      Array<double *> flux;
      Array<double *> ppsol;
      GetPPDarcy(dualmesh, pressuredata, sol, solold, pressurefem, scalarfemx, scalarfemy, betafem, dt, flux, ppsol);

      double pherr = ComputeDiscH1SemiNormError(mesh, solpx, solpy, ppsol, currtime);
      cout << "The H1 SemiNorm error for the post-processed pressure is: " << pherr << endl;

      double preherr = ComputeH1SemiNormError(mesh, solpx, solpy, pressure, currtime);
      cout << "The H1 SemiNorm error for the pressure is: " << preherr << endl;

      // Post-processing for displacements
      Vector ut(2*nv);
      for(int i=0; i<2*nv; i++) ut(i) = ( displacement(i) - displacementold(i) ) / dt;

      Array<double *> uhdotn(dualmesh->GetNumLocalDOF());
      for(int i=0; i<uhdotn.Size(); i++ ) uhdotn[i] = new double[dualmesh->GetNE()];

      Array<double> buhdotn(nv);
      Computeuhdotn(dualmesh, ut, uhdotn, buhdotn);

      Array<double> bdrflux(nv);
      for(int i=0; i<bdrflux.Size(); i++){ bdrflux[i] = 0.0; }
      for(int b=0; b<mesh->GetNBdrs(); b++)
	{
	  for(int i=0; i<=mesh->GetNBdryElem(b); i++)
	    {
	      int numlocaldof = dualmesh->GetNumLocalDOF();
	      int vindex = mesh->GetBdryGindex(b, i);
	      double internalut = 0.0;
	      for(int k=0; k<dualmesh->Getcvnumel(vindex); k++)
		{
		  int el = dualmesh->Getcvelindex(vindex, k);
		  int locvert = dualmesh->Getcvlocalindex(vindex,k);
		  internalut += uhdotn[locvert][el];

		  int locbdryind = (locvert+numlocaldof-1)%numlocaldof; 
		  internalut -= uhdotn[locbdryind][el];
		}
	       bdrflux[vindex] = -internalut; 
	    }
	}    
  
      Vector force(nv);
      for(int k=0; k<nv; k++){ force(k) = pressuredata->GetNodalForceCoeff(k); }
      Array<double> tempbdrflux;
      dualmesh->GetFluxonBdrCVEdges(flux, tempbdrflux, force);
      for(int i=0; i<nv; i++){ bdrflux[i] += tempbdrflux[i]-buhdotn[i]; }
      
      //Data for upwinding
      Array<FData*> fluxdata(2);
      fluxdata[0] = new FData(flux, bdrflux);
      fluxdata[1] = new FData(uhdotn, buhdotn);
      
      //---------------------------------------------------------
      //  Check Local Conservation
      //---------------------------------------------------------
      // cout << "Local Conservation Checking..." << flush;
      Vector cverr(nv);
      cverr = 0.0;
      for(int k=0; k<mesh->GetNE(); k++)
	{
	  Array<int> ind;
	  dualmesh->GetDOF(k, ind);

	  Array<double> areas;     
	  dualmesh->GetAreas(k, areas);
	  if(dualmesh->GetNumLocalDOF()==3)
	    {
	      cverr(ind[0]) += flux[0][k] - flux[2][k] + uhdotn[0][k] - uhdotn[2][k];
	      cverr(ind[1]) += flux[1][k] - flux[0][k] + uhdotn[1][k] - uhdotn[0][k];
	      cverr(ind[2]) += flux[2][k] - flux[1][k] + uhdotn[2][k] - uhdotn[1][k];
	    }
	  else if(dualmesh->GetNumLocalDOF()==4)
	    {
	      cverr(ind[0]) += flux[0][k] - flux[3][k] + uhdotn[0][k] - uhdotn[3][k];
	      cverr(ind[1]) += flux[1][k] - flux[0][k] + uhdotn[1][k] - uhdotn[0][k];
	      cverr(ind[2]) += flux[2][k] - flux[1][k] + uhdotn[2][k] - uhdotn[1][k];
	      cverr(ind[3]) += flux[3][k] - flux[2][k] + uhdotn[3][k] - uhdotn[2][k];
	    }

	  for(int j=0; j<ind.Size(); j++){ cverr(ind[j]) -= pressuredata->GetNodalForceCoeff(ind[j])*areas[j]; }
	}
    
      double maxcverror= 0.0;
      int maxindex;
      for(int k=0; k<nv; k++)
	{
	  cverr(k) +=  bdrflux[k] + buhdotn[k]; 
	  // cout << k << " " << cverr(k) << endl;
	  if(maxcverror<fabs(cverr(k))){ maxcverror = fabs(cverr(k)); maxindex = k;}
	} 
      cout << "Local conservation error: " << maxcverror << " " << "in cv " << maxindex << endl;

      double sumcverr = 0.0;
      for(int k=0; k<nv; k++)
	{
	  sumcverr += fabs(cverr(k));
	}
      cout << "Sum of Local Conservation Error: " << sumcverr << endl;

      //ofstream fileout("lce.out");
      //for(int k=0; k<nv; k++)
      //fileout << k << "\t" << cverr(k) << endl;
      //fileout.close(); 

      sprintf(filename, "localconservationerror_%f.vtk", currtime);
      output.open(filename);
      mesh->ParaviewPrint(output, cverr); 
      output.close();    

      for(int i=0; i<nv; i++){ ; }
  
      Vector u1t(nv);
      Vector u2t(nv);
      for(int k=0; k<nv; k++)
	{
	  u1t(k) = (u1(k)-u1old(k))/dt;
	  u2t(k) = (u2(k)-u2old(k))/dt;
	} 
      //Estimate errors
      cout << "L2error pressure: " << pressurefem->ComputeL2Error(solp, pressure, currtime) << endl;	  
      cout << "L2error u1: " << elasticfem->ComputeL2Error(solu1, u1, currtime) << endl;
      cout << "L2error u2: " << elasticfem->ComputeL2Error(solu2, u2, currtime) << endl;
      cout << "L2error u1t: " << pressurefem->ComputeL2Error(solu1t, u1t, currtime) << endl;
      cout << "L2error u2t: " << elasticfem->ComputeL2Error(solu2t, u2t, currtime) << endl;

   
      pressureold = pressure;
      displacementold = displacement;
      u1old = u1;
      u2old = u2;
      solold = sol;

      
      for(int i=0; i<uhdotn.Size(); i++){ delete []uhdotn[i]; } 
      for(int i=0; i<fluxdata.Size(); i++){ delete fluxdata[i]; }
      for(int i=0; i<ppsol.Size(); i++){ delete []ppsol[i]; }
      for(int i=0; i<flux.Size(); i++){ delete []flux[i]; }
      for(int k=0; k<elfdata.Size(); k++){ delete []elfdata[k]; } 
      for(int k=0; k<scalarcoeffx.Size(); k++){ delete []scalarcoeffx[k]; }
      for(int k=0; k<scalarcoeffy.Size(); k++){ delete []scalarcoeffy[k]; };
      delete []forcedata[0];	  
      delete []corrector[0];
      delete globalmat;
      delete betafem;
      delete betadata;
      delete scalardatax;
      delete scalardatay;
      delete pressurefem;
      delete elasticfem;
      delete scalarfemx;
      delete scalarfemy;
      delete dualmesh;
      cout << endl;
    }

  //---------------------------------------------------------
  //  Freedom!!!!
  //---------------------------------------------------------
  delete mesh;
  delete pressuredata;
  delete elasticdata;
  delete force;
  delete forcex;
  delete forcey;
  delete []permdata[0];

  delete solp;
  delete solu1;
  delete solu2;
  delete solpx;
  delete solpy;
  delete solu1t;
  delete solu2t;

  for(int k=0; k<lamecoeff.Size(); k++)
    delete []lamecoeff[k];

  for(int k=0; k<dbdryval.Size(); k++)
    delete []dbdryval[k];

  for(int k=0; k<nbdryval.Size(); k++)
    delete []nbdryval[k];

  for(int k=0; k<eldbdryval.Size(); k++)
    delete []eldbdryval[k];

  for(int k=0; k<elnbdryval.Size(); k++)
    delete []elnbdryval[k];
 
  for (int k=0; k<elasticboundaryfuncs.Size(); k++)
    delete elasticboundaryfuncs[k]; 

  for (int k=0; k<pressureboundaryfuncs.Size(); k++)
    delete pressureboundaryfuncs[k]; 


  cout << endl << "fin" << endl;
}



