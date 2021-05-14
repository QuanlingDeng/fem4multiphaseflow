#ifndef KRYLOV_METHODS
#define KRYLOV_METHODS

#include "vector.hpp"
#include "matrix.hpp"

/** Conjugate gradient method. Given Matrix A, vector b and initial guess
    x, iteratively solve A x = b. When the default arguments are used
    CG doesn't print current residual and number of iterations, maximum
    number of iterations is 1000, the relative tolerance is 10e-12 and
    the absolute tolerance is 10e-24                                    **/ 

void CG ( const Matrix &A, const Vector &b, Vector &x,
          int print_iter=0, int max_num_iter=1000, 
	  double RTOLERANCE=10e-12, double ATOLERANCE=10e-24);


/** Preconditioned conjugate gradient method. Given Matrix A, preconditioner
    Matrix B, vector b and initial guess x, iteratively solve A x = b. 
    When the default arguments are used PCG doesn't print current residuals
    and number of iterations, maximum number of iterations is 1000, the 
    relative tolerance is 10e-12 and the absolute tolerance is 10e-24. 
    Remark : if no better initial guess is available, the user may set
             it as B b (since not done in PCG routine).                 **/ 

void PCG ( const Matrix &A, const MatrixInverse &B, const Vector &b,Vector &x,
           int print_iter=0, int max_num_iter=1000, 
	   double RTOLERANCE=10e-12, double ATOLERANCE=10e-24);

#endif
