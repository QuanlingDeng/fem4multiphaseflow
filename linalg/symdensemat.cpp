/*
   Implementation of data types symmetric dense matrix &
                                inverse symmetric dense matrix

   Authors: Joanna, Taejong, Veselin
   Date:   01/29/2001
*/

// #define DEBUG

#include <iostream>

#include "vector.hpp"
#include "matrix.hpp"
#include "symdensemat.hpp"

using namespace std;


SymDenseMatrix::SymDenseMatrix (int s) : Matrix(s)
{
   int i;

   data = new double[i=(s+1)*s/2];

   for ( ; i > 0; )
      data[--i] = 0.0;
}

double & SymDenseMatrix::Elem(int i, int j)
{
   return (*this)(i, j); // Use the inline function
}

const double & SymDenseMatrix::Elem(int i, int j) const
{
   return (*this)(i, j); // Use the inline function
}

void SymDenseMatrix::Mult (const Vector & x, Vector & y) const
{
#ifdef DEBUG
   if (( size != x.Size() ) || ( size != y.Size() ))
      cerr << "SymDenseMatrix::Mult";
#endif

  int i,b;
  double *ptr;

  for (i = 0; i < size; i++)
     y(i) = data[i] * x(i);
  
  ptr = data + size;
  for (b = 1; b <size; ptr += size-b,b++)
     for (i = b; i < size; i++)
     {
        y(i)   += ptr[i-b] * x(i-b);
        y(i-b) += ptr[i-b] * x(i);
     }
}

/// Returns a pointer to (approximation) of the matrix inverse.
MatrixInverse * SymDenseMatrix::Inverse() const
{
  return new SymDenseMatrixInverse(*this);
}

SymDenseMatrix::~SymDenseMatrix ()
{
   delete [] data;
}


inline double sqr(double x)
{
   return x*x;
}

SymDenseMatrixInverse::SymDenseMatrixInverse(const SymDenseMatrix & mat)
   : MatrixInverse(mat)
{
   int i, j, k, ofs, ofsi, ofsj;

   data = new double[(size+1)*size/2];

   // perform L D L^T factorization

   for (j = 0; j < size; j++)
   {
      data[j] = mat.data[j];
      for (k = j, ofs = j; k > 0; )
      {
         k--; ofs += size-j+k;
         data[j] -= data[k] * sqr(data[ofs]);
      }
#ifdef DEBUG
      // check if L D L^T factorization fails
      if (data[j] == 0.0)
      {
         cerr << "SymDenseMatrixInverse::SymDenseMatrixInverse" << endl;
         data[j] = 1.0;
      }
#endif
      for (i = j+1, ofs = j+size; i < size; ofs += size-i+j, i++)
      {
         data[ofs] = mat.data[ofs];
         for (k = j, ofsj = j, ofsi = ofs; k > 0; )
         {
            k--;
            ofsj += size-j+k;
            ofsi += size-i+k;
            data[ofs] -= data[ofsj] * data[k] * data[ofsi];
         }
         data[ofs] /= data[j];
      }
   }
}

void SymDenseMatrixInverse::Mult (const Vector & x, Vector & y) const
{
   int i, j;
   double sum;

#ifdef DEBUG
   if (( size != x.Size() ) || ( size != y.Size() ))
      cerr << "SymDenseMatrixInverse::Mult";
#endif

   // y <- L^{-1} x
   for (i = 0; i < size; i++)
   {
      sum = 0.0;
      for(j = 0; j < i; j++)
         sum += data[(i-j)*(2*size-i+j+1)/2 +j] * y(j);
      y(i) = x(i) - sum;
   }

   // y <- D^{-1} y
   for (i = 0; i < size; i++)
   {
      y(i) /= data[i];
   }

   // y<- LT^{-1} y
   for (j = size-1; j >= 0; j--)
   {
      sum = 0.0;
      for(i = j+1; i < size; i++)
         sum += data[(i-j)*(2*size-i+j+1)/2 +j] * y(i);
      y(j) -= sum;
   }
}

SymDenseMatrixInverse::~SymDenseMatrixInverse ()
{
   delete [] data;
}
