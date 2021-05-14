/*
   Implementation of data types sparse matrix, inverse of sparse matrix

   Author: Bacuta, Ginting, Tomov
   Date:   01/31/2001
*/


#include <iostream>
#include <iomanip>

#include "linalg_header.h"

using namespace std;

//============================================================================

/// Creates sparse matrix with fixed number of elements per row.
SparseMatrix::SparseMatrix (int dim, int elements_per_row, int is_rect)
   : Matrix(dim)
{
  int i, j, sum = dim * elements_per_row;

  I = new int[size+1];
  J = new int[sum];
  A = new double[sum];

  for (i = 0; i < sum; i++)
    A[i] = 0.0;

  I[0] = 0;
  for( i=1; i<=size; i++ ){ 
    I[i] = I[i-1] + elements_per_row;
    J[I[i-1]] = (is_rect) ? (-1) : (i-1);
    for( j=I[i-1]+1; j<I[i]; j++ ) { J[j] = -1; }
  }
  width = (is_rect) ? (is_rect) : (size);

#ifdef DEBUG 
  if( I[size] != sum )
     cerr << "SparseMatrix::SparseMatrix()" << endl;
#endif
}

//============================================================================

/// Creates sparse matrix with variable number of elements per row.
SparseMatrix::SparseMatrix (int dim, const int *elements_per_row, int is_rect
			    ) : Matrix( dim)
{
  int i, j, sum = 0;
  I = new int[size+1];

  for( i=0; i<size; i++ ){ sum += elements_per_row[i]; }
  J = new int[sum];

  A = new double[sum];
  for (i = 0; i < sum; i++)
    A[i] = 0.0;

  I[0] = 0;
  for( i=1; i<=size; i++ ){ 
    I[i] = I[i-1] + elements_per_row[i-1];
    J[I[i-1]] = (is_rect) ? (-1) : (i-1);
    for( j=I[i-1]+1; j<I[i]; j++ ) { J[j] = -1; }
  }
  width = (is_rect) ? (is_rect) : (size);

#ifdef DEBUG 
  if( I[size] != sum )
     cerr << "SparseMatrix::SparseMatrix()" << endl;
#endif
}

//============================================================================

SparseMatrix::SparseMatrix(const SparseMatrix *mat, double c) : Matrix(mat->Size())
{
  int *Itemp = mat->GetI();
  int *Jtemp = mat->GetJ();
  double *Atemp = mat->GetValue();

  I = new int[size+1];
  for (int i=0; i<=size; i++) { I[i] = Itemp[i]; }
  
  int sum = I[size];
  J = new int[sum];
  A = new double[sum];
  for (int i=0; i<sum; i++)
    {
      J[i] = Jtemp[i];
      A[i] = c * Atemp[i];
    }
}

//============================================================================

double& SparseMatrix::Elem (int i, int j){
  return operator()(i,j);
}

//============================================================================

const double& SparseMatrix::Elem (int i, int j) const{
  return operator()(i,j);
}

//============================================================================

double& SparseMatrix::operator() (int i, int j)
{
  int k, end = I[i+1];

#ifdef DEBUG
 if ( i>=size || i<0 || j>=width || j<0 )
    cerr << "SparseMatrix::operator() #1" << endl;
#endif

  for(k=I[i]; k<end; k++){
    if (J[k]==j)
      return A[k];
    else if (J[k]==-1){
      J[k] = j; 
      return A[k];
    }
  }

// #ifdef DEBUG
  cerr << "SparseMatrix::operator() #2" << endl;
// #endif
  return A[0];
}

//============================================================================

const double& SparseMatrix::operator() (int i, int j) const
{
  int k, end = I[i+1];

#ifdef DEBUG
 if ( i>=size || i<0 || j>=width || j<0 )
    cerr << "SparseMatrix::operator() #1" << endl;
#endif

  for(k=I[i]; k<end; k++){
    if (J[k]==j)
      return A[k];
    else if (J[k]==-1){
      J[k] = j; 
      return A[k];
    }
  }

// #ifdef DEBUG
  cerr << "SparseMatrix::operator() #2" << endl;
// #endif
  return A[0];
}

//============================================================================

void SparseMatrix::Mult (const Vector & x, Vector & y) const
{
   y = 0.0;
   AddMult (x, y);
}

void SparseMatrix::AddMult (const Vector & x, Vector & y) const
{
#ifdef DEBUG
  if (( width != x.Size() ) || ( size != y.Size() ))
     cerr << "SparseMatrix::AddMult" << endl;
#endif

  int i, j, end;

  for(i=0; i<size; i++){
    end = I[i+1];
    for(j=I[i]; j<end; j++){

#ifdef DEBUG
   if (J[j] == -1)
      cerr << "SparseMatrix::AddMult (SparseMatrix::Finalize not called)"
           << endl;
#endif

       y(i) += A[j]*x(J[j]);
    }
  }
}

double SparseMatrix::InfinityNorm()
{
  Vector x(size);
  Vector y(size);
  x = 1.0;
  y = 0.0;

  int i, j, end;

  for(i=0; i<size; i++){
    end = I[i+1];
    for(j=I[i]; j<end; j++){
      y(i) += fabs(A[j]*x(J[j]));
    }
  }

  return y.InfinityNorm();
}

void SparseMatrix::MultTranspose (const Vector & x, Vector & y) const
{
   y = 0.0;
   AddMultTranspose (x, y);
}

void SparseMatrix::AddMultTranspose (const Vector & x, Vector & y) const
{
#ifdef DEBUG
  if (( size != x.Size() ) || ( width != y.Size() ))
     cerr << "SparseMatrix::AddMultTranspose" << endl;
#endif

  int i, j, end;

  for(i=0; i<size; i++){
    end = I[i+1];
    for(j=I[i]; j<end; j++){

#ifdef DEBUG
   if (J[j] == -1)
      cerr << "SparseMatrix::AddMultTranspose (SparseMatrix::Finalize not called)"
           << endl;
#endif

       y(J[j]) += A[j]*x(i);
    }
  }
}

//============================================================================

void SparseMatrix::Finalize()
{
  int i, j, end, sum = 0, n = 0, newI = 0;

  for(i=0; i<I[size]; i++)
    if (J[i] != -1)
      sum++;
  
  if (sum != I[size])
    {
      int *NewJ = new int[sum];
      double *NewA = new double[sum];
    
      for(i=0; i<size; i++)
	{
	  end = I[i+1];
	  for(j=I[i]; j<end; j++)
	    {
	      if (J[j] == -1) break;
	      NewJ[ n   ] = J[j];
	      NewA[ n++ ] = A[j];
	    }
	  I[i] = newI;
	  newI = n;
	}
      I[size] = sum;

      delete [] J;
      delete [] A;
  
      J = NewJ;
      A = NewA;
#ifdef DEBUG
      if (sum != n)
	cerr << "SparseMatrix::Finalize" << endl;
#endif
    }
}

//============================================================================

MatrixInverse * SparseMatrix::Inverse() const
{
  return (MatrixInverse*) CG();
}

void SparseMatrix::EliminateRow (int row, const double sol, Vector &rhs)
{
   int k, end;

#ifdef DEBUG
   if ( row >= size || row < 0 )
      cerr << "SparseMatrix::EliminateRow ()" << endl;
#endif

   end = I[row+1];

   for(k = I[row]; (k < end) && (J[k] != -1); k++)
   {
      rhs(J[k]) -= sol * A[k];
      J[k] = -1;
   }
}



void SparseMatrix::EliminateRowCol (int rc, const double sol, Vector &rhs)
{
   int k, end;

#ifdef DEBUG
   if ( rc >= size || rc < 0 )
      cerr << "SparseMatrix::EliminateRowCol ()" << endl;
#endif

   end = I[rc+1];
   k = I[rc];
   A[k] = 1.0;  // J[I[rc]] should be equal to rc !!!
   rhs(rc) = sol;
   for (k++; (k < end) && (J[k] != -1); k++)
   {
      int l = I[J[k]], oend = I[J[k]+1];
#ifdef DEBUG
      while (l < oend  && J[l] != rc) l++;
      if (l == oend)
         cerr << "SparseMatrix::EliminateRowCol ()" << endl;
#else
      while (J[l] != rc) l++;
#endif
      rhs(J[k]) -= sol * A[l];
      while ( l < oend && (J[l] = J[l+1]) != -1 ) A[l] = A[l+1], l++;
      J[k] = -1;
   }
}

void SparseMatrix::ClearRow (int rc, const double sol, Vector &rhs)
{
   int k, end;

   end = I[rc+1];
   k = I[rc];
   A[k] = 1.0;  // J[I[rc]] should be equal to rc !!!
   rhs(rc) = sol;
   
   for (k++; (k < end); k++)
   {
     A[k] = 0.0;
   }
   
   
}

//============================================================================

void SparseMatrix::Gauss_Seidel_forw(const Vector &x, Vector &y) const{
  int i, j, end;
  double sum;

  for(i=0; i<size; i++){
    end = I[i+1];
    sum = 0.;
    for(j=(I[i]+1); j<end ; j++)
      sum += A[j] * y( J[j] );
    
    y(i) = (x(i) - sum)/A[I[i]];  
  }
}

//============================================================================

void SparseMatrix::Gauss_Seidel_back(const Vector &x, Vector &y) const{
  int i, j, end;
  double sum = 0.;
  
  for(i=size-1; i >= 0; i--){
    end = I[i+1];
    sum = 0.;
    for(j=(I[i]+1); j<end ; j++)
      sum += A[j] * y( J[j] );
    
    y(i) = (x(i) - sum)/A[I[i]];  
  }
}

//============================================================================

void SparseMatrix::GetSubMatrix(Array<int> &rows, Array<int> &cols,
                                RectangularMatrix &subm)
{
   int i, j, gi, gj, k, end, s, t;
   double a;

   subm.SetSize (rows.Size(), cols.Size());

   for (i = 0; i < rows.Size(); i++)
   {
      if ((gi=rows[i]) < 0) gi = -1-gi, s = -1; else s = 1;
#ifdef DEBUG
      if (gi > size)
         cerr << "SparseMatrix::GetSubMatrix(...) #1" << endl;
#endif
      end = I[gi+1];
      for (j = 0; j < cols.Size(); j++)
      {
         if ((gj=cols[j]) < 0) gj = -1-gj, t = -s; else t = s;
#ifdef DEBUG
         if (gj > width)
            cerr << "SparseMatrix::GetSubMatrix(...) #2" << endl;
#endif
         k = I[gi];
         while (k < end && J[k] != gj && J[k] != -1)  k++;
         a = (k < end && J[k] == gj) ? (A[k]) : (0.0);
         if (t < 0)  a = -a;
         subm(i, j) = a;
      }
   }
}

//============================================================================

void SparseMatrix::SetSubMatrix(Array<int> &rows, Array<int> &cols,
                                RectangularMatrix &subm)
{
   int i, j, gi, gj, k, end, s, t;

   for (i = 0; i < rows.Size(); i++)
   {
      if ((gi=rows[i]) < 0) gi = -1-gi, s = -1; else s = 1;
#ifdef DEBUG
      if (gi > size)
         cerr << "SparseMatrix::SetSubMatrix(...) #1" << endl;
#endif
      end = I[gi+1];
      for (j = 0; j < cols.Size(); j++)
      {
         if ((gj=cols[j]) < 0) gj = -1-gj, t = -s; else t = s;
#ifdef DEBUG
         if (gj > width)
            cerr << "SparseMatrix::SetSubMatrix(...) #2" << endl;
#endif
         k = I[gi];
         while (k < end && J[k] != gj && J[k] != -1)  k++;
         if (k < end)
            J[k] = gj, A[k] = (t >= 0) ? (subm(i, j)) : (-subm(i, j));
         else
            cerr << "SparseMatrix::SetSubMatrix(...) #3" << endl;
      }
   }
}

//============================================================================

void SparseMatrix::AddSubMatrix(const Array<int> &rows, const Array<int> &cols,
                                const RectangularMatrix &subm)
{
   int i, j, gi, gj, k, end, s, t;
   double a;

   for (i = 0; i < rows.Size(); i++)
   {
      if ((gi=rows[i]) < 0) gi = -1-gi, s = -1; else s = 1;
#ifdef DEBUG
      if (gi > size)
         cerr << "SparseMatrix::AddSubMatrix(...) #1" << endl;
#endif
      end = I[gi+1];
      for (j = 0; j < cols.Size(); j++)
      {
         if ((gj=cols[j]) < 0) gj = -1-gj, t = -s; else t = s;
#ifdef DEBUG
         if (gj > width)
            cerr << "SparseMatrix::AddSubMatrix(...) #2" << endl;
#endif
         a = subm(i, j);
         if (t < 0)  a = -a;
         k = I[gi];
         while (k < end && J[k] != gj && J[k] != -1)  k++;
         if (k < end)
            if (J[k] == gj)
               A[k] += a;
            else
               J[k] = gj, A[k] = a;
         else
            cerr << "SparseMatrix::AddSubMatrix(...) #3" << endl;
      }
   }
}

//============================================================================

void SparseMatrix::AddElementMatrix (const Array<int> & dofs1,
                                     const Array<int> & dofs2,
                                     const RectangularMatrix & elemmat)
{
   AddSubMatrix (dofs1, dofs2, elemmat);
}

//============================================================================

void SparseMatrix::Print(ostream & out, int _width){
  int i, j;

  out << setiosflags(ios::scientific | ios::showpos);
  for(i = 0; i < size; i++) {
      out << "[row " << i << "]\n";
      for (j = I[i]; j < I[i+1]; j++) 
	{
	  out << "(" << setw(3) << J[j] << ","<< A[j] << ") ";
	  if ( !((j+1-I[i]) % _width) )
	    out << endl;
	}
      out << endl;
    }
  out << endl;
}

void SparseMatrix::MPrint(ostream & out){
  int i, j;

  out << setiosflags(ios::scientific | ios::showpos);
  for(i = 0; i < size; i++) {
      for (j = I[i]; j < I[i+1]; j++) 
	{
	  out << i+1 << " " << J[j]+1 << " "<< A[j] << endl;
	}
    }
}

void SparseMatrix::PrintCSR(ostream & out)
{
   int i;
   ios::fmtflags old_fmt = out.setf(ios::scientific);
   int old_prec = out.precision(14);

   out << size << '\n';  // number of rows

   for (i = 0; i <= size; i++)
      out << I[i]+1 << '\n';

   for (i = 0; i < I[size]; i++)
      out << J[i]+1 << '\n';

   for (i = 0; i < I[size]; i++)
      out << A[i] << '\n';

   out.precision(old_prec);
   out.setf(old_fmt);
}

//============================================================================

SparseMatrix::~SparseMatrix (){
  delete [] I;
  delete [] J;
  delete [] A;
}
