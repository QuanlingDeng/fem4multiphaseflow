/*
   Implementation of data types for sparse matrix smoothers

   Author: Constantin, Tony, Victor
   Date:   03/04/2001
*/

#include <iostream>
#include "vector.hpp"
#include "matrix.hpp"
#include "sparsemat.hpp"
#include "sparsesmoothers.hpp"

using namespace std;

//============================================================================

/// Create GSSmoother.
GSSmoother::GSSmoother (const SparseMatrix & a): MatrixInverse(a){
}

//============================================================================

/// Matrix vector multiplication with GS Smoother.
void GSSmoother::Mult ( const Vector & x, Vector & y) const{
  y = 0.;

  ((SparseMatrix *)a)->Gauss_Seidel_forw(x, y);
  ((SparseMatrix *)a)->Gauss_Seidel_back(x, y);
}

//============================================================================

/// Destroys the GS Smoother.
GSSmoother::~GSSmoother(){
}

//============================================================================
//============================================================================

/// Create the diagonal smoother.
DSmoother::DSmoother (const SparseMatrix & a, double s): MatrixInverse(a){
  scale = s;
}

//============================================================================

/// Matrix vector multiplication with Diagonal smoother.
void DSmoother::Mult ( const Vector & x, Vector & y) const{
  for(int i=0; i<x.Size(); i++)
    y(i) = scale * x(i) / a->Elem(i, i);
}

//============================================================================

/// Destroys the Diagonal smoother.
DSmoother::~DSmoother(){
}

//============================================================================
