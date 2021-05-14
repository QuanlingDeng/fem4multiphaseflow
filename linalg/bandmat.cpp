/*
   Implementation of data types band matrix, invesre band matrix

   Authors: Joanna, Taejong, Veselin
   Date:   01/29/2001
*/

#include <iostream>
#include <stdlib.h>

#include "vector.hpp"
#include "matrix.hpp"
#include "bandmat.hpp"

using namespace std;

BandMatrix::BandMatrix (int s, int h) : Matrix(s)
{
   int i;

   hbw  = h;
   data = new double[i=(2*h+1)*s];

   for ( ; i > 0; )
      data[--i] = 0.0;
}

double & BandMatrix::Elem(int i, int j)
{
   return (*this)(i, j);
}

const double & BandMatrix::Elem(int i, int j) const
{
   return (*this)(i, j);
}

void BandMatrix::Mult (const Vector & x, Vector & y) const
{
#ifdef DEBUG
   if (( size != x.Size() ) || ( size != y.Size() ))
      cerr << "BandMatrix::Mult";
#endif

  int b, i;
  double *ptr;

  for (i = 0; i < size; i++)
     y(i) = data[i] * x(i);
  ptr = data + size;
  for (b = 1; b <= hbw; b++, ptr += size)
     for (i = b; i < size; i++)
        y(i) += ptr[i-b] * x(i-b);
  for (b = 1; b <= hbw; b++, ptr += size)
     for (i = b; i < size; i++)
        y(i-b) += ptr[i-b] * x(i);
}

/// Returns a pointer to (approximation) of the matrix inverse.
MatrixInverse * BandMatrix::Inverse() const
{
   return new BandMatrixInverse(*this);
}

BandMatrix::~BandMatrix ()
{
   delete [] data;
}


BandMatrixInverse::BandMatrixInverse (const BandMatrix & mat)
   : MatrixInverse(mat)
{
   int i, j, k, u_ofs, ofs;
   double l_ji, u_ij, *ptr;

   hbw  = mat.HalfBandWidth();
   data = new double[(2*hbw+1)*size];

   // perform LU factorization

   u_ofs = hbw * size;
   for (i = 0; i < size; i++)
   {
      data[i] = mat.data[i];
      k = i-hbw; if (k < 0) k = 0;
      ptr = data + (i-k)*size+k;
      for ( ; k < i; k++, ptr -= (size-1))
         data[i] -= ptr[0] * ptr[u_ofs];
#ifdef DEBUG
      // check if LU factorization fails
      if (data[i] == 0.0)
      {
         cerr << "BandMatrixInverse::BandMatrixInverse" << endl;
         data[i] = 1.0;
      }
#endif
      for (j = i+1; (j <= i+hbw) && (j < size); j++)
      {
         ofs = (j-i)*size;
         l_ji = mat.data[ofs+i];
         u_ij = mat.data[u_ofs+ofs+i];
         k = j-hbw; if (k < 0) k = 0;
         ptr = data + (i-k)*size+k;
         for ( ; k < i; k++, ptr -=(size-1))
         {
            l_ji -= ptr[ofs] * ptr[u_ofs];
            u_ij -= ptr[0]   * ptr[u_ofs+ofs];
         }
         u_ij /= data[i];

         data[ofs+i]       = l_ji;
         data[u_ofs+ofs+i] = u_ij;
      }
   }
}

void BandMatrixInverse::Mult (const Vector & x, Vector & y) const
{
   int i, j;
   double sum;
#ifdef DEBUG
   if (( size != x.Size() ) || ( size != y.Size() ))
      cerr << "BandMatrixInverse::Mult";
#endif

   // y <- L^{-1} x
   for(i = 0; i < size; i++)
   {
      sum = 0.0;
      if (i <= hbw)
         j = 0;
      else
         j = i - hbw;
      for ( ; j < i; j++)
         sum += data[(i-j) * size + j] * y(j);
      y(i) = (x(i) - sum) / data[i];
   }

   // y <- U^{-1} y
   for(i = size-1; i >= 0; i--)
   {
      sum = 0.0;
      if (i+hbw < size)
         j = i + hbw;
      else
         j = size-1;
      for ( ; j > i; j--)
         sum += data[(hbw + j - i) * size + i] * y(j);
      y(i) -= sum;
   }
}

BandMatrixInverse::~BandMatrixInverse ()
{
   delete [] data;
}
