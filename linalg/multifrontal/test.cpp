#include <stdio.h>
#include <iostream.h>
#include "multifrontal.h"

extern "C" void exit( int status );

main () 
{

  int N;

  cout << "" << endl;
  cout << "   Enter matrix dimension : ... ";
  cin >> N;
  cout << "" << endl;

  int NELEM = 3*N-2;
  int i,l;
  int *index;
  double *value, *x, *ax;
  index = new int[2*NELEM]; 
  value = new double[NELEM];
  x = new double[N];
  ax = new double[N];
  
  l=0;
  for (i=0; i<N; i++) {
    index[l] = i;
    index[NELEM+l] = i;
    value[l] = 2.0;
    l++;
  }
  for (i=0; i<N-1; i++) {
    index[l] = i;
    index[l+NELEM] = i+1;
    value[l] = -1.0;
    l++;
  }
  for (i=0; i<N-1; i++) {
    index[l] = i+1;
    index[l+NELEM] = i;
    value[l] = -1.0;
    l++;
  }

  cout << "   Setting up the matrix" << endl;
  multifrontal s( N, NELEM, index, value, 3 );

  cout << "   Setting up the vector x = (1, 1, ..., 1)" << endl;
  for (i=0; i<N; i++) ax[i] = 1.0;

  cout << "   Factorize the matrix" << endl; 
  s.LU_fact();

  cout << "   Solve for x :..." << endl;
  s.Solve( ax, x );

  for (i=0 ; i<N; i++) printf("   x(%d) = %e\n", i, x[i] );

  delete []index;
  delete []value;
  delete []x;
  delete []ax;
    
}

