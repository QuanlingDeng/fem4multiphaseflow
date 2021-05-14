#include <iostream>
#include <iomanip>
#include "vector.hpp"
#include "matrix.hpp"
#include "sparsemat.hpp"

/*
  Author: Q. Deng
  Data: 03/01/2014
*/

using namespace std;

//============================================================================

void CG( const Matrix &A, const Vector &b, Vector &x,
	 int print_iter=0, int max_num_iter=1000, 
	 double RTOLERANCE=10e-12, double ATOLERANCE=10e-24){

  int i, dim = x.Size();
  double den, nom, nom0, betanom, alpha, beta, r0;
  Vector r(dim), d(dim), z(dim);

  A.Mult( x, r);
  subtract( b, r, r);                         // r = b - A x
  d = r;
  nom0 = nom = r * r;
  A.Mult( r, z);
  den = r * z;

  if ( (r0 = nom * RTOLERANCE) < ATOLERANCE) r0 = ATOLERANCE;
  if (nom < r0)
    return;
   
  if (den <= 0.0) {
    cout <<"Operator A is not postive definite. (Ar,r) = " 
	 << den << endl;
    return;
  }
  
  if (print_iter == 1)
    cout << "   Iteration : " << setw(3) << 0 << "  (r, r) = "
	 << nom << endl;

  // start iteration                          //  d = r, z = A r
  for(i = 1; i<max_num_iter ;i++) {   
    alpha= nom/den;                           // alpha = (r_o,r_o)/(Ar_o,r_o)
    
    add(x, alpha, d, x);                      //   x =   x + alpha * d 
    add(r, -alpha, z, r);                     // r_n = r_o - alpha * z 
    betanom = r*r;                            // betanom = (r_o, r_o)

    if (print_iter == 1)
      cout << "   Iteration : " << setw(3) << i << "  (r, r) = "
	   << betanom << endl;

    if ( betanom < r0) {
      if (print_iter == 2)
	cout << "Number of CG iterations: " << i << endl;
      else
	if (print_iter == 3)
	  cout << "(r_0, r_0) = " << nom0 << endl
	       << "(r_N, r_N) = " << betanom << endl
	       << "Number of CG iterations: " << i << endl;
      break;
    }
    
    beta = betanom/nom;                       // beta = (r_n,r_n)/(r_o,r_o)
    add(r, beta, d, d);                       //    d = r_n + beta * d
    A.Mult(d, z);                             //    z = A d
    den = d * z;                              //  den = (d , A d)
    nom = betanom;                            //  nom = (r_n, r_n)
  }
}

//============================================================================
