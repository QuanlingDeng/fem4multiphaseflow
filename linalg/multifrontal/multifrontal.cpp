#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "multifrontal.h"

using namespace std;

// Three of the fortran subroutines that are called

//=========================================================================

extern "C" int umd21i_(int *keep, double *cntl, int *icntl) ;

//=========================================================================

extern "C" int umd2fa_(int *N, int *nelem, int *job, 
		       int *transa, int *lvalue, 
		       int *lindex, double *value,
		       int *index, int *keep, double *cntl, 
		       int *icntl, int *info, double *rinfo);

//=========================================================================

extern "C" int umd2so_(int *N, int *job, int *transa, 
		       int *lvalue, int *lindex, double *value,
		       int *index, int *keep, double *B, double *X, 
		       double *work, double *cntl, int *icntl, 
		       int *info, double *rinfo);

//=========================================================================

multifrontal::multifrontal( int dim, int nel, int *ind, double *val, 
			    int mem_mult )
{
  n = N = dim;
  nelem = Nelem = nel;
  Index = ind;
  Value = val; 
  size_factor = mem_mult;
  work = 0;
  index = 0;
  value = 0;
}

//========================================================================

int multifrontal::LU_fact()
{
  int wlen, alen, i; 
  int acopy = Nelem + N + 1;
  job = 0;     // umfpack may overwrite input matrix
  transa = 0;  // factorize A, not its transpose

  // set "index" with extra space, and fortran indexing
  wlen = 14 * N + 8 ; 
  alen = 3 * Nelem + 22 * N ;
  lindex =  wlen + alen + acopy; 
  index = new int[lindex];

  for (i=0; i<Nelem; i++) {
    index[i] = Index[i] + 1;             // need fortran index
    index[i+Nelem] = Index[i+Nelem] + 1;
  }
 
  // set "value" with extra space
  lvalue = size_factor * Nelem;
  value = new double[lvalue];

  for (i=0; i<Nelem; i++)
    value[i] = Value[i]; 

  // umfpack calls
  umd21i_( keep, cntl, icntl ) ;

  icntl[5] = 1; // symmetric sparsity pattern, so prefer diagonal pivot
  icntl[2] = 2; // printing verbosity level

  umd2fa_( &n, &nelem, &job, &transa, &lvalue, &lindex, value, 
	   index, keep, cntl, icntl, info, rinfo );

  
  if (info[0]==0 || info[0]==2) {}
  else  {
    printf( "Factorization unsuccessfull. Info=%d\n",  info[0] );
    if (info[0]==-3) 
      printf("Change parameter size_factor and retry.\n");
    fflush( stdout );
  }

  // set parameters for future solve calls
  work = new double[2*n];

  return (info[0]);
}

//========================================================================

void multifrontal::Solve( const double *b, double *x )
{
  umd2so_( &n, &job, &transa, &lvalue, &lindex, value, index, 
	   keep, (double *)b, x, work, cntl, icntl, info, rinfo ); 

  if (info[0] !=0) {
    if (info[0] == -7) 
      printf( "In splu_solv: Corrupted LU factors given!\n" );
    printf( "Unsuccessfull solve. Quitting.\n" ); exit(-1);
  }
}

//=======================================================================

multifrontal::~multifrontal()
{ 
  delete []index; index = 0;
  delete []value; value = 0;
  delete []work;  work  = 0;
}
