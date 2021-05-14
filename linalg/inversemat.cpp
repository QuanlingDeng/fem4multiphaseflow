/*
   Implementation of data types inverse of a matrix

   Author: Bacuta, Ginting, Tomov
   Date:   01/31/2001
*/


#include <iostream>
#include <math.h>

#include "../general/array.hpp"

#include "vector.hpp"
#include "matrix.hpp"
#include "rectmat.hpp"
#include "inversemat.hpp"
//#include "cgsolver.cpp"
//#include "pcgsolver.cpp"
#include "krylov.hpp"
#include "bicgstab.cpp"
#include "gmres.cpp"


using namespace std;

//============================================================================
//================ Methods for class CGMatrixInverse =========================
//============================================================================

CGMatrixInverse::CGMatrixInverse (const Matrix & a, int printiter, 
				  int maxnumiter, 
				  double rtol, double atol): MatrixInverse(a){
  print_iter = printiter;
  max_num_iter = maxnumiter;
  rtolerance = rtol;
  atolerance = atol;
}

//============================================================================

void CGMatrixInverse::Mult (const Vector & x, Vector & y) const {
#ifdef DEBUG
  if (( size != x.Size() ) || ( size != y.Size() ))
    cerr << "CGMatrixInverse::Mult\n";
#endif

  CG( *a, x, y, print_iter, max_num_iter, rtolerance, atolerance );
}

//============================================================================

CGMatrixInverse::~CGMatrixInverse (){
}


//============================================================================
//================ Methods for class PCGMatrixInverse ========================
//============================================================================

PCGMatrixInverse::PCGMatrixInverse (const Matrix & a, 
				    const MatrixInverse & a_inverse,
				    int printiter, int maxnumiter, 
				    double rtol,double atol):MatrixInverse(a){
  b = &a_inverse;
  print_iter = printiter;
  max_num_iter = maxnumiter;
  rtolerance = rtol;
  atolerance = atol;
}

//============================================================================

void PCGMatrixInverse::Mult (const Vector & x, Vector & y) const {
#ifdef DEBUG
  if (( size != x.Size() ) || ( size != y.Size() ))
    cerr << "PCGMatrixInverse::Mult\n";
#endif

  PCG( *a, *b, x, y, print_iter, max_num_iter, rtolerance, atolerance );
}

//============================================================================

PCGMatrixInverse::~PCGMatrixInverse (){
}


//============================================================================
//============== Methods for class BICGSTABMatrixInverse =====================
//============================================================================

BICGSTABMatrixInverse::BICGSTABMatrixInverse(
   const Matrix &A, const MatrixInverse &M, int printiter,
   int maxnumiter, double rtol, double atol ) : MatrixInverse(A)
{
   Prec = &M;
   print_iter = printiter;
   max_num_iter = maxnumiter;
   rtolerance = rtol;
   atolerance = atol;
}

//============================================================================

void BICGSTABMatrixInverse::Mult(const Vector & x, Vector & y) const
{
   int ok, maxit;
   double tol;

   maxit = max_num_iter;
   tol = rtolerance;

   ok = BiCGSTAB( *a, y, x, *Prec, maxit, tol, atolerance, print_iter);

   if (ok != 0)
      cout << "BICGSTAB: NO Convergence! ( error code = "
           << ok << " )" << endl;
}

//============================================================================
//================ Methods for class GMRESMatrixInverse ======================
//============================================================================

GMRESMatrixInverse::GMRESMatrixInverse(
   const Matrix &A, const MatrixInverse &M, int printiter,
   int maxnumiter, int restart_after, double rtol, double atol ) : MatrixInverse(A)
{
   Prec = &M;
   print_iter = printiter;
   max_num_iter = maxnumiter;
   m = restart_after;
   rtolerance = rtol;
   atolerance = atol;
}

//============================================================================

void GMRESMatrixInverse::Mult(const Vector & x, Vector & y) const
{
   int ok, maxit;
   double tol, atol;

   maxit = max_num_iter;
   tol   = rtolerance;
   atol  = atolerance;

   ok = GMRES( *a, y, x, *Prec, maxit, m, tol, atol, print_iter);

   if (ok != 0)
      cerr << "GMRES: NO Convergence!" << endl;
}

//============================================================================
//================ Methods for class MultiFrontalMatrixInverse ===============
//============================================================================

MultiFrontalMatrixInverse::MultiFrontalMatrixInverse (const SparseMatrix &A,
						      int mem_mult )
{

  int dim, nel, *I, *J;
  double *val;
  dim = A.Size();
  I = A.GetI();
  J = A.GetJ();
  val = A.GetValue();
  nel = I[dim];
  ind = new int[2*nel];
  for( int i=0; i<dim; i++ )
    {
      for( int k=I[i]; k<I[i+1]; k++ )
	{
	  ind[k] = i;
	  ind[k+nel] = J[k];
	}
    }
  linsolver = new multifrontal( dim, nel, ind, val, mem_mult );
}

//============================================================================

void MultiFrontalMatrixInverse::Mult( const Vector &x, Vector &y ) const
{

  linsolver->LU_fact();
  double *xx = x.GetData();
  double *yy = new double[x.Size()];
  linsolver->Solve( xx, yy );
  for( int i=0; i<x.Size(); i++ )
    y(i) = yy[i];
  delete []yy;
}
