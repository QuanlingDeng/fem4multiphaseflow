/*
   Implementation of data types tridiagonal matrix, inverse tridiagonal matrix

   Author: Tzanio, Joachim
   Date:   01/18/2001
*/


#include <iostream>

#include "vector.hpp"
#include "matrix.hpp"
#include "tridiagmat.hpp"

using namespace std;

TriDiagonalMatrix::TriDiagonalMatrix (int s) : Matrix(s)
{
  l = new double[size-1];
  d = new double[size];
  u = new double[size-1];

  for (int i = 0; i <= size-2; i++)
    l[i] = d[i] = u[i]=0;
  d[size-1]=0;
}

double & TriDiagonalMatrix::Elem (int i, int j)
{
  return operator()(i,j);
}

const double & TriDiagonalMatrix::Elem (int i, int j) const
{
  return operator()(i,j);
}

void TriDiagonalMatrix::Mult (const Vector & x, Vector & y) const
{
#ifdef DEBUG
  if (( size != x.Size() ) || ( size != y.Size() ))
    cerr << "TriDiagonalMatrix::Mult";
#endif
  int i;

  y(0) = 
    d[0] * x(0) + 
    u[0] * x(1);
  for (i = 1; i <= (size - 2); i++) 
    y(i) =  
      l[i-1] * x(i-1) + 
      d[i]   * x(i)   + 
      u[i]   * x(i+1);
  y(size-1) = 
    l[size-2] * x(size-2) + 
    d[size-1] * x(size-1);
}

MatrixInverse * TriDiagonalMatrix::Inverse() const
{
  return new TriDiagonalMatrixInverse(*this);
}

TriDiagonalMatrix::~TriDiagonalMatrix ()
{
  delete [] l;
  delete [] d;
  delete [] u;
}



TriDiagonalMatrixInverse::TriDiagonalMatrixInverse (const TriDiagonalMatrix & mat) : MatrixInverse(mat)
{
  l = new double[size-1];
  d = new double[size];
  u = new double[size-1];

  Rebuild(mat);
}

void TriDiagonalMatrixInverse::Rebuild(const TriDiagonalMatrix & mat)
{
  // perform LU factorization

  d[0] = mat.d[0];
  u[0] = mat.u[0];
  
  for (int i=1; i <= size-2; i++)
    {
#ifdef DEBUG
      if (!d[i-1]) 
	cerr << "TriDiagonalMatrixInverse::TriDiagonalMatrixInverse - LU breaks down";
#endif      
 
      l[i-1] = mat.l[i-1] / d[i-1];
      d[i]   = mat.d[i] - u[i-1] * l[i-1];
      u[i]   = mat.u[i];
    }
 
#ifdef DEBUG
  if (!d[size-2]) 
    cerr << "TriDiagonalMatrixInverse::TriDiagonalMatrixInverse - LU breaks down";
#endif

  l[size-2] = mat.l[size-2] / d[size-2];
  d[size-1] = mat.d[size-1] - u[size-2] * l[size-2];
}

void TriDiagonalMatrixInverse::Mult (const Vector & x, Vector & y) const
{
#ifdef DEBUG
  if (( size != x.Size() ) || ( size != y.Size() ))
    cerr << "TriDiagonalMatrixInverse::Mult";
#endif

  int i;

  // y <- L^{-1} x
  y(0) = x(0);
  for (i = 1; i < size; i++)
    y(i) = x(i) - l[i-1] * y(i-1);

  // y <- U^{-1} y 
  y(size-1) = y(size-1) / d[size-1];
  for (i = size-2; i >= 0; i--)
    y(i) = (y(i) - u[i] * y(i+1)) / d[i];
}

TriDiagonalMatrixInverse::~TriDiagonalMatrixInverse ()
{
  delete [] l;
  delete [] d;
  delete [] u;
}
