/* 
   Implementation of data types rectangular matrix, inverse dense matrix

   Author: Kristi, Yanqiu
   Date:   01/24/2001
*/


#include <iostream>
#include <iomanip>

#include "vector.hpp"
#include "matrix.hpp"
#include "rectmat.hpp"

using namespace std;

RectangularMatrix::RectangularMatrix () : Matrix(0)
{
   height = 0;
   data = 0;
}

RectangularMatrix::RectangularMatrix (int s) : Matrix(s)
{
  height = s;

  data = new double[s*s];

  for (int i = 0; i < s; i++)
    for (int j = 0; j < s; j++)
      (*this)(i,j) = 0;
}

RectangularMatrix::RectangularMatrix (int m, int n) : Matrix(n)
{
  height = m;

  data = new double[m*n];

  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      (*this)(i,j) = 0;
}

RectangularMatrix::RectangularMatrix (const RectangularMatrix & mat, char ch) : Matrix(mat.height)
{
  if (ch == 't')
    {
      height = mat.size;
      
      data = new double[size*height];
      
      for (int i = 0; i < height; i++)
	for (int j = 0; j < size; j++)
	  (*this)(i,j) = mat(j,i);  
    }
}

void RectangularMatrix::SetSize(int s)
{
   if (Size() == s && Height() == s)
      return;
   if (data != 0)
      delete [] data;
   size = height = s;
   if (s > 0)
   {
      int ss = s*s;
      data = new double[ss];
      for (int i = 0; i < ss; i++)
         data[i] = 0.0;
   }
   else
      data = 0;
}

void RectangularMatrix::SetSize(int h, int w)
{
   if (Size() == w && Height() == h)
      return;
   if (data != 0)
      delete [] data;
   size = w;
   height = h;
   if (h > 0 && w > 0)
   {
      int hw = h*w;
      data = new double[hw];
      for (int i = 0; i < hw; i++)
         data[i] = 0.0;
   }
   else
      data = 0;
}

double & RectangularMatrix::Elem (int i, int j)
{
  return (*this)(i,j);
}

const double & RectangularMatrix::Elem (int i, int j) const
{
  return (*this)(i,j);
}

void RectangularMatrix::Mult (const Vector & x, Vector & y) const
{
#ifdef DEBUG
  if ( size != x.Size() || height != y.Size() )
    {
      cerr << "RectangularMatrix::Mult\n";
      return;
    }
#endif

  for (int i = 0; i < height; i++)
    {
      y(i) = 0;
      for (int j = 0; j < size; j++)
	y(i) += (*this)(i,j) * x(j);
    }
}

void RectangularMatrix::MultTranspose (const Vector & x, Vector & y) const
{
#ifdef DEBUG
  if ( height != x.Size() || size != y.Size() )
    {
      cerr << "RectangularMatrix::MultTranspose\n";
      return;
    }
#endif

  for (int i = 0; i < size; i++)
    {
      y(i) = 0;
      for (int j = 0; j < height; j++)
	y(i) += (*this)(j,i) * x(j);
    }
}

MatrixInverse * RectangularMatrix::Inverse() const
{
  return new DenseMatrixInverse(*this);
}

double RectangularMatrix::Det()
{
  double det;

#ifdef DEBUG
  if ( Height() != 1 || Size()!= 1 )
    if ( Height() != 2 || Size()!= 2 )
      if ( Height() != 3 || Size()!= 3 )
        cerr << "RectangularMatrix::Det" << endl;
#endif

  switch (Size()) {
  case 1:
    det = (*this)(0,0);
    break;
  case 2:
    det = (*this)(0,0) * (*this)(1,1) - (*this)(0,1) * (*this)(1,0);
    break;
  case 3:
    det = (*this)(0,0)*
      ( (*this)(1,1)*(*this)(2,2)-(*this)(1,2)*(*this)(2,1) ) +
      (*this)(0,1)*
      ( (*this)(2,0)*(*this)(1,2)-(*this)(1,0)*(*this)(2,2) ) +
      (*this)(0,2)*
      ( (*this)(1,0)*(*this)(2,1)-(*this)(2,0)*(*this)(1,1) );
    break;
  }
  return det;
}

void RectangularMatrix::Add(const double c, RectangularMatrix & A)
{
   for(int i=0;i<Height();i++)
      for(int j=0;j<Size();j++)
         (*this)(i,j) += c * A(i,j);
}

RectangularMatrix &RectangularMatrix::operator=(double c)
{
   int s=Size()*Height();
   if (data != 0)
      for(int i = 0; i < s; i++)
         data[i] = c;
   return *this;
}

RectangularMatrix &RectangularMatrix::operator+=(RectangularMatrix &m)
{
   int i, j;

   for (i = 0; i < height; i++)
      for (j = 0; j < size; j++)
         (*this)(i, j) += m(i, j);

   return *this;
}

RectangularMatrix &RectangularMatrix::operator-=(RectangularMatrix &m)
{
   int i, j;

   for (i = 0; i < height; i++)
      for (j = 0; j < size; j++)
         (*this)(i, j) -= m(i, j);

   return *this;
}

void RectangularMatrix::Neg()
{
   int i, hw = Height() * Size();

   for (i = 0; i < hw; i++)
      data[i] = -data[i];
}

void RectangularMatrix::Invert()
{
#ifdef DEBUG
   if (Size() != Height())
   {
      cerr << "RectangularMatrix::Invert()" << endl;
      return;
   }
#endif

   int c, i, j, n = Size();
   double a, b;

   for (c = 0; c < n; c++)
   {
#ifdef DEBUG
      if ((*this)(c, c) == 0.0)
         cerr << "RectangularMatrix::Invert() : division by zero" << endl;
#endif
      a = (*this)(c, c) = 1.0 / (*this)(c, c);
      for (j = 0; j < c; j++)
         (*this)(c, j) *= a;
      for (j = c+1; j < n; j++)
         (*this)(c, j) *= a;
      for (i = 0; i < c; i++)
      {
         (*this)(i, c) = a * (b = -(*this)(i, c));
         for (j = 0; j < c; j++)
            (*this)(i, j) += b * (*this)(c, j);
         for (j = c+1; j < n; j++)
            (*this)(i, j) += b * (*this)(c, j);
      }
      for (i = c+1; i < n; i++)
      {
         (*this)(i, c) = a * (b = -(*this)(i, c));
         for (j = 0; j < c; j++)
            (*this)(i, j) += b * (*this)(c, j);
         for (j = c+1; j < n; j++)
            (*this)(i, j) += b * (*this)(c, j);
      }
   }
}

void RectangularMatrix::Print (ostream & out, int width)
{
  // output flags = scientific + show sign 
  out << setiosflags(ios::scientific | ios::showpos); 
  for (int i = 0; i < height; i++) 
    {
      out << "[row " << i << "]\n";
      for (int j = 0; j < size; j++) 
	{
	  out << Elem(i,j) << " ";
	  if ( !((j+1) % width) )
	    out << endl;
	}
      out << endl;
    }
  out << endl;
}

RectangularMatrix::~RectangularMatrix()
{
   if (data != 0)
      delete [] data;
}



void Mult (const RectangularMatrix & b, 
	   const RectangularMatrix & c, 
	   RectangularMatrix & a)
{
#ifdef DEBUG
  if ( a.height != b.height || a.size != c.size || b.size != c.height )
    cerr << "Mult\n";
#endif

  for (int i = 0; i < a.height; i++)
    for (int j = 0; j < a.size; j++)
      {
	a(i,j) = 0;
	for (int k = 0; k < b.size; k++)
	  a(i,j) += b(i,k) * c(k,j);
      }
}

void CalcInverse(RectangularMatrix & a, RectangularMatrix & inva)
{
#ifdef DEBUG
   if ( (a.Size() != a.Height()) || ( (a.Height()!= 1) && (a.Height()!= 2)
                                      && (a.Height()!= 3) ) )
      cerr << "RectangularMatrix::CalcInverse(...)" << endl;
#endif

   double t = 1. / a.Det() ;

   switch (a.Height()) {
   case 1:
     inva(0,0) = 1.0 / a(0,0);
     break;
   case 2:
     inva(0,0) = a(1,1) * t ;
     inva(0,1) = -a(0,1) * t ;
     inva(1,0) = -a(1,0) * t ;
     inva(1,1) = a(0,0) * t ;
     break;
   case 3:
     inva(0,0) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))*t;
     inva(0,1) = (a(0,2)*a(2,1)-a(0,1)*a(2,2))*t;
     inva(0,2) = (a(0,1)*a(1,2)-a(0,2)*a(1,1))*t;
     
     inva(1,0) = (a(1,2)*a(2,0)-a(1,0)*a(2,2))*t;
     inva(1,1) = (a(0,0)*a(2,2)-a(0,2)*a(2,0))*t;
     inva(1,2) = (a(0,2)*a(1,0)-a(0,0)*a(1,2))*t;
     
     inva(2,0) = (a(1,0)*a(2,1)-a(1,1)*a(2,0))*t;
     inva(2,1) = (a(0,1)*a(2,0)-a(0,0)*a(2,1))*t;
     inva(2,2) = (a(0,0)*a(1,1)-a(0,1)*a(1,0))*t;
     break;
   }
}

void MultAAt( RectangularMatrix & a, RectangularMatrix & aat)
{
  double temp ;
  for(int i=0;i<a.Height();i++)
    for(int j=0;j<a.Height();j++)
    {
      temp = 0. ;
      for(int k=0;k<a.Size();k++)
	temp += a(i,k) * a(j,k) ;
      aat(i,j) = temp;
    }
}

void MultABt ( RectangularMatrix & A, RectangularMatrix & B,
               RectangularMatrix & ABt )
{
   int i, j, k;
   double d;

   for (i = 0; i < A.Height(); i++)
      for (j = 0; j < B.Height(); j++)
      {
         d = 0.0;
         for (k = 0; k < A.Size(); k++)
            d += A(i, k) * B(j, k);
         ABt(i, j) = d;
      }
}

void MultAtB ( RectangularMatrix & A, RectangularMatrix & B,
               RectangularMatrix & AtB )
{
   int i, j, k;
   double d;

   for (i = 0; i < A.Size(); i++)
      for (j = 0; j < B.Size(); j++)
      {
         d = 0.0;
         for (k = 0; k < A.Height(); k++)
            d += A(k, i) * B(k, j);
         AtB(i, j) = d;
      }
}

void MultVVt( Vector & v, RectangularMatrix & vvt)
{
   for(int i = 0; i < v.Size(); i++)
      for(int j = 0; j <= i; j++)
      {
         vvt(i,j) = vvt(j,i) = v(i) * v(j);
      }
}

void MultVWt(Vector &v, Vector &w, RectangularMatrix &VWt)
{
   int i, j;

#ifdef DEBUG
   if (v.Size() != VWt.Height() || w.Size() != VWt.Size())
      cerr << "MultVWt (...)" << endl;
#endif

   for (i = 0; i < v.Size(); i++)
      for (j = 0; j < w.Size(); j++)
         VWt(i, j) = v(i) * w(j);
}





DenseMatrixInverse::DenseMatrixInverse (const RectangularMatrix & mat) : MatrixInverse(mat)
{
#ifdef DEBUG
  if ( mat.height != mat.size )
    cerr << "DenseMatrixInverse::DenseMatrixInverse";
#endif
  data = new double[size*size];

  int i,j,k;

  // perform LU factorization.
  for (i = 0; i < size; i++)
    {      
#ifdef DEBUG
      if ( i>0 && !data[i-1+size*(i-1)] )
	cerr << "DenseMatrixInverse::DenseMatrixInverse";
#endif
      for (j = 0; j < i; j++)
	{
	  data[i+size*j] = mat(i,j);
	  for (k = 0; k < j; k++)
	    data[i+size*j] -= data[i+size*k] * data[k+size*j];
	  data[i+size*j] /= data[j+size*j];
	}
      
      for (j = i; j < size; j++)
	{
	  data[i+size*j] = mat(i,j);
	  for (k = 0; k < i; k++)
	    data[i+size*j] -= data[i+size*k] * data[k+size*j];
	}
    }
}

void DenseMatrixInverse::Mult (const Vector & x, Vector & y) const
{
  int i,j;
  
  // y <- L^{-1} x
  for (i = 0; i < size; i++)
    {
      y(i) = x(i);
      for (j = 0; j < i; j++)
	y(i) -= data[i+size*j] * y(j);
    }

  // y <- U^{-1} y
  for (i = size-1; i >= 0; i--)
    {
      for (j = i+1; j < size; j++)
	y(i) -= data[i+size*j] * y(j);
#ifdef DEBUG
      if ( !(data[i+size*i]) )
	cerr << "DenseMatrixInverse::Mult";
#endif
      y(i) /= data[i+size*i] ;
    }
}

DenseMatrixInverse::~DenseMatrixInverse()
{
  delete [] data;
}
