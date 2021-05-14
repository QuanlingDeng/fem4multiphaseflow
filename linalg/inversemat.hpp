#ifndef FILE_INVERSEMAT
#define FILE_INVERSEMAT

/*
   Data types for inverse of a matrix

   Author: Bacuta, Ginting, Tomov
   Date:   01/31/2001
*/

using namespace std;

#include "./multifrontal/multifrontal.h"
#include "sparsemat.hpp"

class MatrixInverse;

/// Data type for matrix inverse (done by CG solver).
   
class CGMatrixInverse : public MatrixInverse
{
private:
  
  /// Pointer to the matrix.
  int print_iter;
  int max_num_iter;
  double rtolerance, atolerance;

public:

  /// Creates the matrix inverse. 
  CGMatrixInverse (const Matrix & a, int printiter = 0, 
                   int maxnumiter = 1000, 
	           double rtol = 10e-12, double atol = 10e-24);
  
  /// Matrix vector multiplication with CG Matrix Inverse. 
  virtual void Mult (const Vector & x, Vector & y) const;

  /// Destroys the matrix inverse.  
  virtual ~CGMatrixInverse ();
};


//============================================================================
//============================================================================

/// Data type for matrix inverse (done by PCG solver).
   
class PCGMatrixInverse : public MatrixInverse
{
private:
  
  int print_iter;
  int max_num_iter;  
  double rtolerance, atolerance;

  /// Pointer to the Preconditioner.  
  const MatrixInverse *b;

public:

  /// Creates matrix inverse. 
  PCGMatrixInverse (const Matrix & a, const MatrixInverse & a_inverse,
                    int printiter = 0, int maxnumiter = 1000, 
		    double rtol = 10e-12, double atol = 10e-24);
  
  /// Matrix vector multiplication with PCG Matrix Inverse. 
  virtual void Mult (const Vector & x, Vector & y) const;

  /// Destroys the matrix inverse.  
  virtual ~PCGMatrixInverse ();
};

//============================================================================
//============================================================================

/// Data type for matrix inverse (done by GMRES solver).

class GMRESMatrixInverse : public MatrixInverse
{
   private:

      int print_iter;
      int max_num_iter;
      int m;
      double rtolerance, atolerance;

      const MatrixInverse * Prec;

   public:

      GMRESMatrixInverse (const Matrix &, const MatrixInverse &,
                          int printiter = 0, int maxnumiter = 10000,
			  int restart_after = 40, 
			  double rtol = 10e-12, double atol = 10e-24);

      virtual void Mult (const Vector & x, Vector & y) const;

      virtual ~GMRESMatrixInverse () { };
};

//============================================================================
//============================================================================

/// Data type for matrix inverse (done by BICGSTAB solver).

class BICGSTABMatrixInverse : public MatrixInverse
{
   private:

      int print_iter;
      int max_num_iter;
      double rtolerance, atolerance;

      const MatrixInverse * Prec;

   public:

      BICGSTABMatrixInverse (const Matrix &, const MatrixInverse &,
                             int printiter = 0, int maxnumiter = 1000,
                             double rtol = 10e-12, double atol = 10e-24);

      virtual void Mult (const Vector & x, Vector & y) const;

      virtual ~BICGSTABMatrixInverse () { };
};

//============================================================================
//============================================================================

/// Data type for matrix inverse (done by MultiFrontal (direct) solver).

class MultiFrontalMatrixInverse
{
  private:
    multifrontal *linsolver;
    int *ind;

  public:
  MultiFrontalMatrixInverse (const SparseMatrix &, int mem_mult );

  void Mult (const Vector & x, Vector & y) const;
  
  ~MultiFrontalMatrixInverse ()
  {
    delete []ind;
    delete linsolver;
  }

};

#endif
