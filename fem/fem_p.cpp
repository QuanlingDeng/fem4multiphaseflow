#include "mpi_fem_header.h"

FEM_p::FEM_p(Mesh_p *_mesh, Data *_data)
{
  mesh = _mesh;
  data = _data;
}

//============================================================================

FEM_p::~FEM_p()
{
  delete mat;
}

//============================================================================

void FEM_p::Assemble()
{

  // initialize rhs to zero
  rhs = 0.0;

  //mat is already initialized to zero in the constructor so it is not a problem
  // using += without first zeroing the appropriate entries

  int startval = mat->Startval();
  int endval = mat->Endval();

  Array<int> ind(NumLocalDOF);

  if (data->EllipticPartExists())
    {
      for(int i=0; i<mesh->GetNE(); i++)
	{
	  // get element dofs for fill in
	  GetElementDOF(i, ind);
	  
	  //setup locmat
	  Array<double *> locmat(NumLocalDOF);
	  for(int j=0; j<NumLocalDOF; j++) { locmat[j] = new double[NumLocalDOF]; }
	  
	  //get local system and fill in global matrix
	  ComputeEllipticLocalSystem(i, ind, locmat);
	  
	  for(int j=0; j<NumLocalDOF; j++)
	    for(int k=0; k<NumLocalDOF; k++)
	      {
		if(startval <= ind[j] && ind[j] < endval)
		  mat->Elem(ind[j]-startval, ind[k]) += locmat[j][k];
	      }
	  for(int j=0; j<NumLocalDOF; j++) { delete []locmat[j]; }
	}
    }
  
  if (data->ReactionPartExists())
    for(int i=0; i<mesh->GetNE(); i++)
      {
	// get element dofs for fill in
	GetElementDOF(i, ind);

	//setup locmat
	Array<double *> locmat(NumLocalDOF);
	for(int j=0; j<NumLocalDOF; j++) { locmat[j] = new double[NumLocalDOF]; }

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


  mat->Finalize();
}

//============================================================================

void FEM_p::GetElementVertices(int i, Array<int> &ind)
{
  mesh->GetElementVertices(i, ind);
}

//============================================================================

void FEM_p::GetBdrElementVertices(int i, Array<int> &ind)
{
  mesh->GetBdrElementVertices(i, ind);
}

//============================================================================

void FEM_p::Solve(Vector &sol, int maxiter, int printopt)
{
  //sol.SetSize(NumGlobalDOF);
  //sol = 0.0;
  /// Call the linear algebra packet to solve the system and store the result in sol.

  PCG_p(*mat, rhs, sol, printopt, maxiter, 1e-7);
}

//============================================================================

void FEM_p::SetDirichletBoundary(Array<bool> &bdry_is_dirichlet, const Vector &Dval)
{

  int matstartval = mat->Startval();
  int matendval = mat->Endval();
  int index;
  Element *el;
  for (int i=0; i<GetNBE(); i++)
    {
      el = GetBdrElement(i);
      int attr = el->GetAttribute();

      if (bdry_is_dirichlet[attr])
	{
	  Array<int> row;
	  GetBdrElementDOF(i, row);
	  for (int ii=0; ii<row.Size(); ii++)
	    {
	      index = row[ii];
	      if(matstartval <= index && index < matendval)
		{
		  mat->EliminateRowCol(index-matstartval);
		}
	      rhs(index) = Dval(index);
	    }
	}
    }
  
  int tn = mpi_Size();
  int mn = mpi_Rank();
  /*
  for(int k=0; k<tn; k++)
    {
      if(mn == k)
	{
	   mat->Print();
	   for(int i=0; i<rhs.Size(); i++)
	     cout << i << " " << rhs(i) << endl;
	}
      mpi_Barrier();
      }*/
}
