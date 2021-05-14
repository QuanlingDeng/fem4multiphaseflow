#include "fem_header.h"

FEM::FEM(Mesh *_mesh, Data *_data)
{
  mesh = _mesh;
  data = _data;
}

//============================================================================

FEM::~FEM()
{
  delete mat;
}

//============================================================================ 

void FEM::Assemble()
{
  // initialize rhs to zero
  rhs = 0.0;

  //mat is already initialized to zero in the constructor so it is not a problem
  // using += without first zeroing the appropriate entries

  Array<int> ind(NumLocalDOF);

  if (data->EllipticPartExists())
    {
      for(int i=0; i<mesh->GetNE(); i++)
	{
	  // get element dofs for fill in
	  GetElementDOF(i, ind);

	  //setup locmat
	  Array<double *> locmat(NumLocalDOF);
	  for(int j=0; j<NumLocalDOF; j++) 
	    { 
	      locmat[j] = new double[NumLocalDOF]; 
	    }

	  //get local system and fill in global matrix
	  ComputeEllipticLocalSystem(i, ind, locmat);
	  
	  for(int j=0; j<NumLocalDOF; j++)
	    {
	      for(int k=0; k<NumLocalDOF; k++)
		{
		  mat->Elem(ind[j], ind[k]) += locmat[j][k];
		}
	    }

	  for(int j=0; j<NumLocalDOF; j++) { delete []locmat[j]; }
	}
    }
  
  if (data->AdvectionPartExists())
    for(int i=0; i<mesh->GetNE(); i++)
      {
	// get element dofs for fill in
	GetElementDOF(i, ind);

	//setup locmat
	Array<double *> locmat(NumLocalDOF);
	for(int j=0; j<NumLocalDOF; j++) { locmat[j] = new double[NumLocalDOF]; }

	//get local system and fill in global matrix
	ComputeAdvectionLocalSystem(i, ind, locmat);

	for(int j=0; j<NumLocalDOF; j++)
	  for(int k=0; k<NumLocalDOF; k++)
	    mat->Elem(ind[j], ind[k]) += locmat[j][k];

	for(int j=0; j<NumLocalDOF; j++) { delete []locmat[j]; }
      }

  if (data->ConservativeAdvectionPartExists())
    for(int i=0; i<mesh->GetNE(); i++)
      {
	// get element dofs for fill in
	GetElementDOF(i, ind);

	//setup locmat
	Array<double *> locmat(NumLocalDOF);
	for(int j=0; j<NumLocalDOF; j++) { locmat[j] = new double[NumLocalDOF]; }

	//get local system and fill in global matrix
	ComputeConservativeAdvectionLocalSystem(i, ind, locmat);

	for(int j=0; j<NumLocalDOF; j++)
	  for(int k=0; k<NumLocalDOF; k++)
	    mat->Elem(ind[j], ind[k]) += locmat[j][k];

	for(int j=0; j<NumLocalDOF; j++) { delete []locmat[j]; }
      }

  if (data->ReactionPartExists())
    for(int i=0; i<mesh->GetNE(); i++)
      {
	// get element dofs for fill in
	GetElementDOF(i, ind);

	//setup locmat
	Array<double *> locmat(NumLocalDOF);
	for(int j=0; j<NumLocalDOF; j++) 
	  { 
	    locmat[j] = new double[NumLocalDOF]; 
	  }

	//get local system and fill in global matrix
	ComputeReactionLocalSystem(i, ind, locmat);

	for(int j=0; j<NumLocalDOF; j++)
	  for(int k=0; k<NumLocalDOF; k++)
	    mat->Elem(ind[j], ind[k]) += locmat[j][k];
 
	for(int j=0; j<NumLocalDOF; j++) { delete []locmat[j]; }
      }

  if (data->ForcePartExists())
    for(int i=0; i<mesh->GetNE(); i++)
      {
	// get element dofs for fill in
	GetElementDOF(i, ind);

	//setup locrhs
	Array<double> locrhs(NumLocalDOF);
      
	// get local system and fill in global matrix
	ComputeForceLocalSystem(i, ind, locrhs);

	for(int k=0; k<NumLocalDOF; k++)
	  {
	    rhs(ind[k]) += locrhs[k];
	  }
      }

  if (data->NeumannPartExists())
    {
      for(int i=0; i<mesh->GetNBE(); i++)
	{
	  // get boundary element dofs for fill in
	  GetBdrElementDOF(i, ind);
	  
	  //setup locrhs
	  Array<double> locrhs(ind.Size());
      
	  // get local system and fill in global matrix
	  ComputeNeumannLocalSystem(i, ind, locrhs);
	  for(int k=0; k<ind.Size(); k++)
	    {
	      rhs(ind[k]) += locrhs[k];
	    }
	}
    }
  
  
  if (data->BdrNeumannExists())
    {
      int bdryindex = 0;
      for(int b=0; b<mesh->GetNBdrs(); b++)
	{
	  int numel = mesh->GetNBdryElem(b);
	  if(data->BdryNeumann(b))
	    {
	      for(int i=0; i<numel; i++)
		{
		  // get boundary element dofs for fill in
		  GetBdrElementDOF(bdryindex+i, ind);

		  Array<int> tempind(ind.Size());
		  for(int k=0; k<ind.Size(); k++)
		    {
		      tempind[k] = i + k;
		    }
		  
		  //setup locrhs
		  Array<double> locrhs(ind.Size());
		  
		  // get local system and fill in global matrix
		  ComputeNeumannLocalSystem(b, i, tempind, locrhs);
		  for(int k=0; k<ind.Size(); k++)
		    {
		      rhs(ind[k]) += locrhs[k];
		    }
		}
	    }
	  bdryindex += numel;
	}
    }

  if (data->BdrRobinExists())
    {
      int bdryindex = 0;
      for(int b=0; b<mesh->GetNBdrs(); b++)
	{
	  int numel = mesh->GetNBdryElem(b);
	  if(data->BdryRobin(b))
	    {
	      for(int i=0; i<numel; i++)
		{
		  // get boundary element dofs for fill in
		  GetBdrElementDOF(bdryindex+i, ind);
		  		  
		  //setup locrhs
		  Array<double> locrhs(ind.Size());

		  //setup locmat
		  Array<double *> locmat(ind.Size());
		  //setup locind
		  Array<int> tempind(ind.Size());
		  for(int k=0; k<ind.Size(); k++) 
		    { 
		      locmat[k] = new double[ind.Size()]; 
		      tempind[k] = i + k;
		    }

		  // get local system and fill in global matrix
		  ComputeRobinLocalSystem(i, b, tempind, locmat, locrhs);

		  for(int j=0; j<ind.Size(); j++)
		    {
		      rhs(ind[j]) += locrhs[j];
		      for(int k=0; k<ind.Size(); k++)
			{
			  mat->Elem(ind[j], ind[k]) += locmat[j][k];
			}
		    }
		  for(int j=0; j<ind.Size(); j++) { delete []locmat[j]; }
		  
		}
	    }
	  bdryindex += numel;
	}
    }

  if (data->RobinPartExists())
    for(int i=0; i<mesh->GetNBE(); i++)
      {
	//setup locrhs
	Array<double> locrhs;
	    
	//setup locmat 
	Array<double *> locmat;
	/*
	if(strcmp(Type, "Stokes2DCR"))
	  {
	    GetBdrElementVertices(i, ind);

	    // use vert_to_ele to find element containing the boundary element
	    int elem = -1;
	    for(int k=0; k<mesh->GetNV(); i++)
	      {
		for(int p=0; p<mesh->GetVertToEle(k,0)+1; p++)
		  {
		    if(ind[0] == mesh->GetVertToEle(k,p))
		      {
			elem = mesh->GetVertToEle(k,p);
			break;
		      }		 
		  } 
		  if(elem != -1)
		    break;
	      }

	    locmat.SetSize(7);
	    locrhs.SetSize(7);

	    for(int j=0; j<7; j++) 
	      { 
		locmat[j] = new double[7]; 
	      }

	    GetElementDOF(elem, ind);

	    // get local system and fill in global matrix
	    ComputeRobinLocalSystem(elem, ind, locmat, locrhs);
	    
	    for(int j=0; j<ind.Size(); j++)
	      for(int k=0; k<ind.Size(); k++)
		mat->Elem(ind[j], ind[k]) += locmat[j][k];
	    
	    for(int j=0; j<ind.Size(); j++) { delete []locmat[j]; }
	    
	    
	    for(int k=0; k<ind.Size(); k++)
	      {
		rhs(ind[k]) += locrhs[k];
	      }
	    
	  }
	else
	  {
	*/
	GetBdrElementDOF(i, ind);

	//setup locmat 
	locmat.SetSize(ind.Size());
	locrhs.SetSize(ind.Size());

	for(int j=0; j<ind.Size(); j++) 
	  { 
	    locmat[j] = new double[ind.Size()]; 
	  }

	// get local system and fill in global matrix
	ComputeRobinLocalSystem(i, ind, locmat, locrhs);
	    
	for(int j=0; j<ind.Size(); j++)
	  for(int k=0; k<ind.Size(); k++)
	    mat->Elem(ind[j], ind[k]) += locmat[j][k];
	    
	for(int j=0; j<ind.Size(); j++) { delete []locmat[j]; }
	    
	    
	for(int k=0; k<ind.Size(); k++)
	  {
	    rhs(ind[k]) += locrhs[k];
	  }
	    
	// }
      }


  if (data->StablePartExists())
    for(int i=0; i<mesh->GetNE(); i++)
      { 
	// get element dofs for fill in
	GetElementDOF(i, ind);

	//setup locmat
	Array<double *> locmat(NumLocalDOF);
	for(int j=0; j<NumLocalDOF; j++) { locmat[j] = new double[NumLocalDOF]; }

	//get local system and fill in global matrix
	ComputeStableLocalSystem(i, ind, locmat);

	for(int j=0; j<NumLocalDOF; j++)
	  for(int k=0; k<NumLocalDOF; k++)
	    mat->Elem(ind[j], ind[k]) += locmat[j][k];

	for(int j=0; j<NumLocalDOF; j++) { delete []locmat[j]; }
      }

  if (data->StableRHSPartExists())
    for(int i=0; i<mesh->GetNE(); i++)
      {
	// get element dofs for fill in
	GetElementDOF(i, ind);

	//setup locrhs
	Array<double> locrhs(NumLocalDOF);
      
	// get local system and fill in global rhs
	ComputeStableRHSLocalSystem(i, ind, locrhs);

	for(int k=0; k<NumLocalDOF; k++)
	  {
	    rhs(ind[k]) += locrhs[k];
	  }
      }

  
  if (data->PointSourcePartExists())
    {
      for(int i=0; i<data->GetNumberofPointSources(); i++)
	{
	  int index = data->GetPointSourceIndex(i);
	  double value = data->GetPointSourceValue(i);
	  rhs(index) += value;
	}
    }
  

  if (data->BdrDirichletExists())
    {
      for(int j=0; j<data->GetNumberofBoundaryConditions(); j++)
	{
	  int nbdrys = mesh->GetNBdrs();
	  for(int b=0; b<nbdrys; b++)
	    {
	      int bb = b + j*nbdrys;
	      if(data->BdryDirichlet(bb))
		{
		  int numel = mesh->GetNBdryElem(b);
		  for(int i=0; i<=numel; i++)
		    {
		      int index = mesh->GetBdryGindex(b, i) + j*mesh->GetNV();
		      double dval = data->GetDirichletBdryVal(bb, i);
		      mat->ClearRow(index, dval, rhs); 
		    }	   
		}
	    }
	  
	  if(NumLocalDOF>4 && data->GetNumberofBoundaryConditions()==1)
	    {
	      int bdryelindex = 0;
	      for(int b=0; b<nbdrys; b++)
		{
		  int bb = b + j*nbdrys;
		  int numel = mesh->GetNBdryElem(b); 
		  if(data->BdryDirichlet(bb))
		    {
		      //number of elements associated with boundary b
		      for(int i=0; i<numel; i++)
			{
			  GetBdrElementDOF(bdryelindex+i, ind);
			  for(int k=2; k<ind.Size(); k++)
			    { 
			      int index = ind[k];
			      int locbdrindex = numel + (ind.Size()-2)*i + k - 1;
			      double dval = data->GetDirichletBdryVal(bb, locbdrindex);
			      mat->ClearRow(index, dval, rhs); 	
			    }
			}   
		    }
		  bdryelindex += numel;
		}
	    }
	}
    }

  //mat->Finalize();
  //mat->Print();
  //rhs.Print();
}

//============================================================================

void FEM::AssembleMatrix()
{
  // initialize rhs to zero
  rhs = 0.0;

  //mat is already initialized to zero in the constructor so it is not a problem
  // using += without first zeroing the appropriate entries
  Array<int> ind(NumLocalDOF);

  if (data->EllipticPartExists())
    {
      for(int i=0; i<mesh->GetNE(); i++)
	{
	  // get element dofs for fill in
	  GetElementDOF(i, ind);

	  //setup locmat
	  Array<double *> locmat(NumLocalDOF);
	  for(int j=0; j<NumLocalDOF; j++) 
	    { 
	      locmat[j] = new double[NumLocalDOF]; 
	    }

	  //get local system and fill in global matrix
	  ComputeEllipticLocalSystem(i, ind, locmat);
	  
	  for(int j=0; j<NumLocalDOF; j++)
	    {
	      for(int k=0; k<NumLocalDOF; k++)
		{
		  mat->Elem(ind[j], ind[k]) += locmat[j][k];
		}
	    }

	  for(int j=0; j<NumLocalDOF; j++) { delete []locmat[j]; }
	}
    }
  
  if (data->AdvectionPartExists())
    {
      for(int i=0; i<mesh->GetNE(); i++)
	{
	  // get element dofs for fill in
	  GetElementDOF(i, ind);

	  //setup locmat
	  Array<double *> locmat(NumLocalDOF);
	  for(int j=0; j<NumLocalDOF; j++) { locmat[j] = new double[NumLocalDOF]; }

	  //get local system and fill in global matrix
	  ComputeAdvectionLocalSystem(i, ind, locmat);

	  for(int j=0; j<NumLocalDOF; j++)
	    for(int k=0; k<NumLocalDOF; k++)
	      mat->Elem(ind[j], ind[k]) += locmat[j][k];

	  for(int j=0; j<NumLocalDOF; j++) { delete []locmat[j]; }
	}
    }

  if (data->ConservativeAdvectionPartExists())
    {
      for(int i=0; i<mesh->GetNE(); i++)
	{
	  // get element dofs for fill in
	  GetElementDOF(i, ind);

	  //setup locmat
	  Array<double *> locmat(NumLocalDOF);
	  for(int j=0; j<NumLocalDOF; j++) { locmat[j] = new double[NumLocalDOF]; }

	  //get local system and fill in global matrix
	  ComputeConservativeAdvectionLocalSystem(i, ind, locmat);

	  for(int j=0; j<NumLocalDOF; j++)
	    for(int k=0; k<NumLocalDOF; k++)
	      mat->Elem(ind[j], ind[k]) += locmat[j][k];

	  for(int j=0; j<NumLocalDOF; j++) { delete []locmat[j]; }
	}
    }

  if (data->ReactionPartExists())
    {
      for(int i=0; i<mesh->GetNE(); i++)
	{
	  // get element dofs for fill in
	  GetElementDOF(i, ind);

	  //setup locmat
	  Array<double *> locmat(NumLocalDOF);
	  for(int j=0; j<NumLocalDOF; j++) 
	    { 
	      locmat[j] = new double[NumLocalDOF]; 
	    }

	  //get local system and fill in global matrix
	  ComputeReactionLocalSystem(i, ind, locmat);

	  for(int j=0; j<NumLocalDOF; j++)
	    for(int k=0; k<NumLocalDOF; k++)
	      mat->Elem(ind[j], ind[k]) += locmat[j][k];
 
	  for(int j=0; j<NumLocalDOF; j++) { delete []locmat[j]; }
	}
    }

  if (data->ForcePartExists())
    {
      for(int i=0; i<mesh->GetNE(); i++)
	{
	  // get element dofs for fill in
	  GetElementDOF(i, ind);
	  
	  //setup locrhs
	  Array<double> locrhs(NumLocalDOF);
	  
	  // get local system and fill in global matrix
	  ComputeForceLocalSystem(i, ind, locrhs);
	  
	  for(int k=0; k<NumLocalDOF; k++)
	    {
	      rhs(ind[k]) += locrhs[k];
	    }
	}
    }

  if (data->BdrNeumannExists())
    {
      int bdryindex = 0;
      for(int b=0; b<mesh->GetNBdrs(); b++)
	{
	  int numel = mesh->GetNBdryElem(b);
	  bool test;
	  for(int j=0; j<data->GetNumberofBoundaryConditions(); j++)
	    {
	      test = (data->BdryNeumann(b+j*mesh->GetNBdrs())) ? true : false;
	      if(test ==true){ break;} 
	    }
	  if(test)
	    {
	      for(int i=0; i<numel; i++)
		{
		  // get boundary element dofs for fill in
		  GetBdrElementDOF(bdryindex+i, ind);

		  Array<int> tempind(ind.Size());
		  for(int k=0; k<ind.Size(); k++)
		    {
		      tempind[k] = i + k;
		    }
		  
		  //setup locrhs
		  Array<double> locrhs(ind.Size());
		  
		  // get local system and fill in global matrix
		  ComputeNeumannLocalSystem(b, i, tempind, locrhs);
		  for(int k=0; k<ind.Size(); k++)
		    {
		      rhs(ind[k]) += locrhs[k];
		    }
		}
	    }
	  bdryindex += numel;
	}
    }
}

//============================================================================

void FEM::UpdateRHS()
{
  rhs = 0.0;
  if (data->ForcePartExists())
    {
      Array<int> ind(NumLocalDOF);
      for(int i=0; i<mesh->GetNE(); i++)
	{
	  // get element dofs for fill in
	  GetElementDOF(i, ind);
	  
	  //setup locrhs
	  Array<double> locrhs(NumLocalDOF);
	  
	  // get local system and fill in global matrix
	  ComputeForceLocalSystem(i, ind, locrhs);
	  
	  for(int k=0; k<NumLocalDOF; k++)
	    {
	      rhs(ind[k]) += locrhs[k];
	    }
	}

      if (data->BdrNeumannExists())
	{
	  int bdryindex = 0;
	  for(int b=0; b<mesh->GetNBdrs(); b++)
	    {
	      int numel = mesh->GetNBdryElem(b);
	      if(data->BdryNeumann(b))
		{
		  for(int i=0; i<numel; i++)
		    {
		      // get boundary element dofs for fill in
		      GetBdrElementDOF(bdryindex+i, ind);

		      Array<int> tempind(ind.Size());
		      for(int k=0; k<ind.Size(); k++)
			{
			  tempind[k] = i + k;
			}
		  
		      //setup locrhs
		      Array<double> locrhs(ind.Size());
		  
		      // get local system and fill in global matrix
		      ComputeNeumannLocalSystem(b, i, tempind, locrhs);
		      for(int k=0; k<ind.Size(); k++)
			{
			  rhs(ind[k]) += locrhs[k];
			}
		    }
		}
	      bdryindex += numel;
	    }
	}

      if (data->BdrRobinExists())
	{
	  int bdryindex = 0;
	  for(int b=0; b<mesh->GetNBdrs(); b++)
	    {
	      int numel = mesh->GetNBdryElem(b);
	      if(data->BdryRobin(b))
		{
		  for(int i=0; i<numel; i++)
		    {
		      // get boundary element dofs for fill in
		      GetBdrElementDOF(bdryindex+i, ind);
		  		  
		      //setup locrhs
		      Array<double> locrhs(ind.Size());

		      //setup locmat
		      Array<double *> locmat(ind.Size());
		      //setup locind
		      Array<int> tempind(ind.Size());
		      for(int k=0; k<ind.Size(); k++) 
			{ 
			  locmat[k] = new double[ind.Size()]; 
			  tempind[k] = i + k;
			}

		      // get local system and fill in global matrix
		      ComputeRobinLocalSystem(i, b, tempind, locmat, locrhs);

		      for(int j=0; j<ind.Size(); j++)
			{
			  rhs(ind[j]) += locrhs[j];
			  for(int k=0; k<ind.Size(); k++)
			    {
			      mat->Elem(ind[j], ind[k]) += locmat[j][k];
			    }
			}
		      for(int j=0; j<ind.Size(); j++) { delete []locmat[j]; }
		  
		    }
		}
	      bdryindex += numel;
	    }
	}


      if (data->BdrDirichletExists())
	{
	  for(int j=0; j<data->GetNumberofBoundaryConditions(); j++)
	    {
	      int nbdrys = mesh->GetNBdrs();
	      for(int b=0; b<nbdrys; b++)
		{
		  int bb = b + j*nbdrys;
		  if(data->BdryDirichlet(bb))
		    {
		      int numel = mesh->GetNBdryElem(b);
		      for(int i=0; i<=numel; i++)
			{
			  int index = mesh->GetBdryGindex(b, i) + j*mesh->GetNV();
			  double dval = data->GetDirichletBdryVal(bb, i);
			  mat->ClearRow(index, dval, rhs); 
			}	   
		    }
		}
	      
	      if(NumLocalDOF>4 && data->GetNumberofBoundaryConditions()==1)
		{
		  int bdryelindex = 0;
		  for(int b=0; b<nbdrys; b++)
		    {
		      int bb = b + j*nbdrys;
		      int numel = mesh->GetNBdryElem(b); 
		      if(data->BdryDirichlet(bb))
			{
			  //number of elements associated with boundary b
			  for(int i=0; i<numel; i++)
			    {
			      GetBdrElementDOF(bdryelindex+i, ind);
			      for(int k=2; k<ind.Size(); k++)
				{
				  int index = ind[k];
				  int locbdrindex = numel + (ind.Size()-2)*i + k - 1;
				  double dval = data->GetDirichletBdryVal(bb, locbdrindex);
				  mat->ClearRow(index, dval, rhs); 			    
				}
			    }   
			}
		      bdryelindex += numel;
		    }
		}
	    }
	}
    }
}

//============================================================================

void FEM::GetElementVertices(int i, Array<int> &ind)
{
  mesh->GetElementVertices(i, ind);
}

//============================================================================

void FEM::GetBdrElementVertices(int i, Array<int> &ind)
{
  mesh->GetBdrElementVertices(i, ind);
}

//============================================================================

void FEM::Solve(Vector &sol, char *linalgsolv, int maxit, int printit, 
		double rtol, double atol)
{
  
  /*
  ofstream out;
  out.open("matrix.out");
  mat->MPrint(out);
  out.close();
  */
  sol = 0.0;
  if (!strcmp(linalgsolv, "multifrontal")) 
    {
      int mem_fact = 200;
      MultiFrontalMatrixInverse *InvMat = 
	new MultiFrontalMatrixInverse(*mat, mem_fact);
      InvMat->Mult(rhs, sol);
      delete InvMat;
    }
  else 
    {      
      MatrixInverse *prec = new GSSmoother(*mat);
      MatrixInverse *invmat;
      
      if (!strcmp(linalgsolv, "pcg"))
	invmat = new PCGMatrixInverse(*mat, *prec, printit, maxit, rtol, atol);
      if (!strcmp(linalgsolv, "bicgstab"))
	invmat = new BICGSTABMatrixInverse(*mat, *prec, printit, maxit, 
					   rtol, atol);
      if (!strcmp(linalgsolv, "gmres"))
	invmat = new GMRESMatrixInverse(*mat, *prec, printit, maxit, 40, 
					rtol, atol);
      
      invmat->Mult(rhs, sol);
      delete invmat;
      delete prec;
    }
}
