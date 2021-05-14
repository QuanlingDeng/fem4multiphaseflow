/*
   Implementation of data types sparse matrix, inverse of sparse matrix

   Author: Bacuta, Ginting, Tomov
   Date:   01/31/2001
*/


#include <iostream>
#include <iomanip>

#include "vector.hpp"
#include "matrix.hpp"
#include "inversemat.hpp"
#include "rectmat.hpp"

using namespace std;

//============================================================================

Matrix::~Matrix() {}

//============================================================================

void Matrix::Finalize() {}

//============================================================================

void Matrix::AddElementMatrix (const Array<int> & dofs, const Matrix & elemmat)
{
  for (int i = 0; i < elemmat.Size(); i++)
    for (int j = 0; j < elemmat.Size(); j++)
       if (dofs[i] >= 0)
          if (dofs[j] >= 0)
             Elem(dofs[i], dofs[j]) += elemmat.Elem(i,j);
          else
             Elem(dofs[i], -1-dofs[j]) -= elemmat.Elem(i,j);
       else
          if (dofs[j] >= 0)
             Elem(-1-dofs[i], dofs[j]) -= elemmat.Elem(i,j);
          else
             Elem(-1-dofs[i], -1-dofs[j]) += elemmat.Elem(i,j);
}

void Matrix::AddElementMatrix (
   const Array<int> & dofs1, const Array<int> & dofs2,
   const RectangularMatrix & elemmat )
{
  for (int i = 0; i < elemmat.Height(); i++)
    for (int j = 0; j < elemmat.Size(); j++)
       if (dofs1[i] >= 0)
          if (dofs2[j] >= 0)
             Elem(dofs1[i], dofs2[j]) += elemmat(i,j);
          else
             Elem(dofs1[i], -1-dofs2[j]) -= elemmat(i,j);
       else
          if (dofs2[j] >= 0)
             Elem(-1-dofs1[i], dofs2[j]) -= elemmat(i,j);
          else
             Elem(-1-dofs1[i], -1-dofs2[j]) += elemmat(i,j);
}

//============================================================================

void Matrix::Print (ostream & out, int width)
{
  // output flags = scientific + show sign 
  out << setiosflags(ios::scientific | ios::showpos); 
  for (int i = 0; i < size; i++) 
    {
      out << "[row " << i << "]\n";
      for (int j = 0; j < size; j++) 
	{
	  out << Elem(i,j) << " ";
	  if ( !((j+1) % width) )
	    out << endl;
	}
      out << endl;
    }
  out << endl;
}

//============================================================================

MatrixInverse::~MatrixInverse() {}

//============================================================================

CGMatrixInverse * Matrix::CG(int printiter, int maxnumiter, 
			     double rtol, double atol) const {
  return new CGMatrixInverse(*this, printiter, maxnumiter, rtol, atol);
}

//============================================================================

PCGMatrixInverse * Matrix::PCG(const MatrixInverse & b, int printiter, 
			    int maxnumiter, double rtol, double atol) const {
  return new PCGMatrixInverse(*this, b, printiter, maxnumiter, rtol, atol);
}

//============================================================================
