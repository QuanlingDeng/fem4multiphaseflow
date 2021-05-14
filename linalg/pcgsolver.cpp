#include <iostream>
#include <iomanip>
#include "vector.hpp"
#include "matrix.hpp"
#include "sparsemat.hpp"

/*
  Author: Q. Deng
  Data: 03/16/2014
*/

using namespace std;

//============================================================================

void PCG ( const Matrix &A, const MatrixInverse &B, 
	   const Vector &b, Vector &x,
	   int print_iter=0, int max_num_iter=1000, 
	   double RTOLERANCE=10e-12, double ATOLERANCE=10e-24){
  
  int i, dim = x.Size();
  double r0, den, nom, nom0, betanom, alpha, beta;
  Vector r(dim), d(dim), z(dim);

  A.Mult(x, r);                               //    r = A x
  subtract(b, r, r);                          //    r = b  - r
  B.Mult(r, z);                               //    z = B r
  d = z;
  nom0 = nom = z * r;
  A.Mult(d, z);
  den = z * d;

  if ( (r0 = nom * RTOLERANCE) < ATOLERANCE) r0 = ATOLERANCE;
  if (nom < r0) 
    return;

  if (den <= 0.0) {
    cout << "Negative or zero denominator in step 0 of PCG. Exiting!\n";
    return;
  }
  
  if (print_iter == 1)
    cout << "   Iteration : " << setw(3) << 0 << "  (B r, r) = "
	 << nom << endl;

  // start iteration
  for(i= 1; i < max_num_iter ; i++) {
    alpha = nom/den;
    add(x, alpha, d, x);                  //  x = x + alpha d
    add(r,-alpha, z, r);                  //  r = r - alpha z
    B.Mult( r, z);                        //  z = B r
    betanom = r * z;
    
    if (print_iter == 1)
      cout << "   Iteration : " << setw(3) << i << "  (B r, r) = "
	   << betanom << endl;

    if ( betanom < r0) {
      if (print_iter == 2)
	cout << "Number of PCG iterations: " << i << endl;
      else
	if (print_iter == 3)
	  cout << "(B r_0, r_0) = " << nom0 << endl
	       << "(B r_N, r_N) = " << betanom << endl
	       << "Number of PCG iterations: " << i << endl;
      break;
    }
    
    beta = betanom/nom;
    add(z, beta, d, d);                   //  d = z + beta d
    A.Mult(d, z);
    den = d * z;
    nom = betanom;
  }
}

//============================================================================
