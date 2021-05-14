#include "timedependent_header.h"

ParabolicProblem::ParabolicProblem()
{
  // method = _method;
}

//============================================================================

ParabolicProblem::~ParabolicProblem()
{
  
}

//============================================================================

void ParabolicProblem::TimeMarchSolve(Array<double *> &data, 
				      Array<double *> &bdrydata, 
				      Array<bool *> &dirichlet,
				      Array<double *> &bdryreaction,
				      Vector &initsol, Array<double> &time,
				      Vector &sol, Vector &temp)
{
  /*
  method->UpdateData(data);
  method->UpdateBdryData(dirichlet, bdrydata, bdryreaction);
  method->Assemble();

  Vector *force = method->GetForce();
  const SparseMatrix *tempmat = method->GetMatrix();

  Array<double> diag(initsol.Size());
  Vector rhs(initsol.Size());
  sol.SetSize(initsol.Size());
  rhs = initsol;
  sol = initsol;
  int printit = 0;
  int maxit = 10000;
  double rtol = 1.0e-32;
  double atol = 1.0e-32;

  for (int i=1; i<time.Size(); i++)
    {
      double timestep = time[i] - time[i-1];
      rhs = sol;
      SetDiagonalMatrix(dirichlet, timestep, diag, rhs);

      SparseMatrix *mat = new SparseMatrix(tempmat);
      
      for (int k=0; k<rhs.Size(); k++)
	{
	  mat->Elem(k,k) += diag[k];
	  rhs(k) += force->Elem(k);
	}
      
      MatrixInverse *prec = new GSSmoother(*mat);
      //MatrixInverse *InvMat = new PCGMatrixInverse(*mat, *prec, printit, maxit, rtol, atol);
      MatrixInverse *InvMat = new BICGSTABMatrixInverse(*mat, *prec, printit, maxit, rtol, atol);
      InvMat->Mult(rhs, sol);
      delete InvMat;
      delete prec;
      delete mat;

    }
  */
}

//============================================================================

void ParabolicProblem::SetDiagonalMatrix(Array<bool *> &dirichlet, double ts, Array<double> &diag,
					 Vector &rhs)
{
  /*
  Array<int> nums;
  Array<double> length;
  method->GetMeshParameters(nums, length);
  int NumVertices = (nums[0]+1)*(nums[1]+1);
  Vector temp(NumVertices);
  temp = rhs;

  int k = 0;
  for (int j=0; j<=nums[1]; j++)
    {
      double cvy = (j>0 && j<nums[1]) ? length[1] : 0.5 * length[1];
      for (int i=0; i<=nums[0]; i++)
	{
	  double cvx = (i>0 && i<nums[0]) ? length[0] : 0.5 * length[0];
	  double fact = cvx * cvy / ts;
	  diag[k] = fact;
	  rhs(k) = temp(k) * fact;
	  k++;
	}
    }

  // Bottom boundary.
  for (int i=0; i<=nums[0]; i++) { if (dirichlet[0][i]) { diag[i] = 0.0; rhs(i) = 0.0; } }

  // Right boundary
  for (int j=0; j<=nums[1]; j++)
    if (dirichlet[1][j]) { int ind = nums[0]+j*(nums[0]+1); diag[ind] = 0.0; rhs(ind) = 0.0; }

  // Top boundary
  for (int i=0; i<=nums[0]; i++)
    if (dirichlet[2][i]) { int ind = i + nums[1]*(nums[0]+1); diag[ind] = 0.0; rhs(ind) = 0.0; }

  // Left boundary
  for (int j=0; j<=nums[1]; j++)
    if (dirichlet[3][j]) { int ind = j*(nums[0]+1); diag[ind] = 0.0; rhs(ind) = 0.0; }
  */
}

//============================================================================

void ParabolicProblem::TimeMarchSolveConsistentMass(Array<double *> &data, Array<double *> &bdrydata,
						    Array<bool *> &dirichlet, Array<double *> &bdryreaction,
                                                    Vector &initsol, Array<double> &time, Vector &sol, Vector &temp)
{
  /*
  SparseMatrix *tempmat = new SparseMatrix(method->GetSize(), 10);
  SparseMatrix *M = new SparseMatrix(method->GetSize(), 10);
  method->UpdateData(data);
  method->UpdateBdryData(dirichlet, bdrydata, bdryreaction);
  method->AssembleMatrix(tempmat);
  
  Vector rhs(initsol.Size());
  sol.SetSize(initsol.Size());
  rhs = initsol;
  sol = initsol;
  int printit = 0;
  int maxit = 10000;
  double rtol = 1.0e-32;
  double atol = 1.0e-32;

  for (int i=1; i<time.Size(); i++)
    {
      double timestep = time[i] - time[i-1];

      //Assemble ====================================
      Array<int> ind;
      RectangularMatrix locmat(4);

      for (int j=0; j<method->GetNy(); j++)
	{
	  for (int i=0; i<method->GetNx(); i++)
	    {
	      method->GetElementVertices(i, j, ind);
	      for (int k=0; k<4; k++)
		{
		  for (int l=0; l<4; l++) 
		    { 
		      M->Elem(ind[k],ind[l]) = 0.0;
		    }
		}
	    }
	}

      for (int j=0; j<method->GetNy(); j++)
	{
	  for (int i=0; i<method->GetNx(); i++)
	    {
	      method->GetElementVertices(i, j, ind);
	      method->ComputeLocalConsistentMass(i, j, locmat);       
	      for (int k=0; k<4; k++)
		{
		  for (int l=0; l<4; l++) 
		    { 
		      tempmat->Elem(ind[k], ind[l]) += locmat(k,l)/timestep;
		      M->Elem(ind[k],ind[l]) += locmat(k,l)/timestep;
		    }
		}
	    }
	}       
 
      M->Finalize();
      tempmat->Finalize();
      method->SetMatrix(tempmat);
      
      Vector force(initsol.Size());
      method->AssembleForce(force);

      M->Mult(initsol, rhs);

      for (int k=0; k<rhs.Size(); k++)
	{
	  rhs(k) += force(k);
	}

      method->SetForce(rhs);
      method->SetBoundaryData();

      char linalgsolver[256];
      strcpy(linalgsolver, "pcg");
      method->Solve(sol, linalgsolver, maxit, printit, rtol, atol, false);

      for(int i=0; i<initsol.Size(); i++)
	{
	  temp(i) = (initsol(i)-sol(i))/timestep;
	}

      delete M;
      delete tempmat;
    }
  */
}
