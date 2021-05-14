/*
   Implementation of data types Table.

   Author: Bacuta, Ginting, Tomov
   Date:   03/10/2001
*/


using namespace std;

#include <iostream>
#include <iomanip>

#include "array.hpp"
#include "table.hpp"
#include "mergesort.hpp"


//============================================================================

/// Creates Table with fixed number of connections.
Table::Table (int dim, int connections_per_row) {
  int i, j, sum = dim * connections_per_row;

  size = dim;
  I = new int[size+1];
  J = new int[sum];

  I[0] = 0;
  for( i=1; i<=size; i++ ){ 
    I[i] = I[i-1] + connections_per_row; 
    for( j=I[i-1]; j<I[i]; j++ ) { J[j] = -1; }
  } 

#ifdef DEBUG 
  if( I[size] != sum )
    cerr << "Table::Table()\n";
#endif
}

//============================================================================

/// Creates Table with variable number of connections.
Table::Table (int dim, const int *connections_per_row) {
  int i, j, sum = 0;
  size = dim;
  I = new int[size+1];

  for( i=0; i<size; i++ ){ sum += connections_per_row[i]; }
  J = new int[sum];

  I[0] = 0;
  for( i=1; i<=size; i++ ){ 
    I[i] = I[i-1] + connections_per_row[i-1];
    for( j=I[i-1]; j<I[i]; j++ ) { J[j] = -1; }
  } 

#ifdef DEBUG 
  if( I[size] != sum )
    cerr << "Table::Table()\n";
#endif
}

//============================================================================

Table::Table (int nrows, int *partitioning)
{
   size = nrows;

   I = new int[size+1];
   J = new int[size];

   for (int i = 0; i < size; i++)
   {
      I[i] = i;
      J[i] = partitioning[i];
   }
   I[size] = size;
}

//============================================================================

void Table::MakeI (int nrows)
{
   size = nrows;
   I = new int[nrows+1];
   J = NULL;

   for (int i = 0; i <= nrows; i++)
      I[i] = 0;
}

//============================================================================

void Table::MakeJ()
{
   int i, j, k;

   for (k = i = 0; i < size; i++)
      j = I[i], I[i] = k, k += j;

   J = new int[I[size]=k];
}

//============================================================================

void Table::ShiftUpI()
{
   for (int i = size; i > 0; i--)
      I[i] = I[i-1];
   I[0] = 0;
}

//============================================================================

void Table::SetSize(int dim, int connections_per_row)
{
   int i, j, sum = dim * connections_per_row;

   SetDims (dim, sum);

   if (size > 0)
   {
      I[0] = 0;
      for( i=1; i<=size; i++ )
      {
         I[i] = I[i-1] + connections_per_row; 
         for( j=I[i-1]; j<I[i]; j++ ) { J[j] = -1; }
      }
   }

#ifdef DEBUG 
    if( I[size] != sum )
      cerr << "Table::SetSize()\n";
#endif
}

//============================================================================

void Table::SetDims(int rows, int nnz)
{
  int i, j;

  j = (I) ? (I[size]) : (0);
  if (size != rows)
  {
     size = rows;
     if (I) delete [] I;
     I = (rows > 0) ? (new int[rows+1]) : (NULL);
  }

  if (j != nnz)
  {
     if (J) delete [] J;
     J = (nnz > 0) ? (new int[nnz]) : (NULL);
  }

  if (size > 0)
  {
     I[0] = 0;
     I[size] = nnz;
  }
}

//============================================================================

int Table::operator() (int i, int j) const
{
  if ( i>=size || i<0 )
    return -1;

  int k, end = I[i+1];

  for(k=I[i]; k<end; k++)
    if (J[k]==j)
      return k;
    else if (J[k]==-1){
#ifdef DEBUG
      //      cerr << "Connection (i, j) = ("  << i << ", " << j 
      //   << ") wasn't found in the table! Entry -1 found! \n";
#endif 
      return -1;
    }

#ifdef DEBUG
  // cerr << "Connection (i, j) = ("  << i << ", " << j 
  //   << ") wasn't found in the table!\n";
#endif 
  return -1;
}

//============================================================================

void Table::GetRow(int i, Array<int> &row) const{
  int j, n = 0;

  for(j=I[i]; j<I[i+1]; j++)
    if (J[j]!=-1) n++;
    else break;

  row.SetSize( n);
  for(j=0; j<n; j++)
    row[j] = J[I[i]+j];
}

//============================================================================

void Table::SetRow(int i, Array<int> &row)
{
   int j, k, n;

#ifdef DEBUG
   if (i < 0 || i >= size)
   {
      cerr << "Table::SetRow (...) 1" << endl;
      return;
   }
#endif

   n = row.Size();
#ifdef DEBUG
   if (n > I[i+1]-I[i])
   {
      cerr << "Table::SetRow (...) 2" << endl;
      n = I[i+1]-I[i];
   }
#endif
   for (j = I[i], k = 0; k < n; j++, k++)
      J[j] = row[k];
   for (n = I[i+1]; j < n; j++)
      J[j] = -1;
}

//============================================================================

void Table::MakeTransposeOf (int *ii, int *jj, int nrows, int _ncols)
{
   if (nrows == 0)
   {
      SetDims (0, 0);
      return;
   }

   int i, j, ncol, nnz = ii[nrows];

   if (_ncols < 0)
   {
      ncol = 0;
      for (i = 0; i < nnz; i++)
         if (ncol < (j=jj[i]))
            ncol = j;
      ncol++;
   }
   else
      ncol = _ncols;

   SetDims (ncol, nnz);

   int *i_tr, *j_tr;

   i_tr = I;
   j_tr = J;

   for (i = 0; i <= ncol; i++)
      i_tr[i] = 0;
   for (i = 0; i < nnz; i++)
      i_tr[jj[i]+1]++;
   for (i = 1; i < ncol; i++)
      i_tr[i+1] += i_tr[i];

#ifdef DEBUG
   if (i_tr[ncol] != nnz)
      cerr << "Table::MakeTransposeOf (...) 1" << endl;
   if (i_tr[0] != 0)
      cerr << "Table::MakeTransposeOf (...) 2" << endl;
#endif

   for (i = 0; i < nrows; i++)
      for (j = ii[i]; j < ii[i+1]; j++)
         j_tr[i_tr[jj[j]]++] = i;
   for (i = ncol; i > 0; i--)
      i_tr[i] = i_tr[i-1];
   i_tr[0] = 0;
}

//============================================================================

void Table::GetTranspose (Table &tr, int _ncols)
{
   tr.MakeTransposeOf (I, J, size, _ncols);
}

//============================================================================

void Table::SortColumnIndexes()
{
   int i, j, max, *aux;

   max = 0;
   for (i = 0; i < size; i++)
      if ( max < (j=I[i+1]-I[i]) )
         max = j;

   aux = new int[max];

   for (i = 0; i < size; i++)
      mergesort (I[i+1]-I[i], J+I[i], aux);

   delete [] aux;
}

//============================================================================

int Table::Push( int i, int j ){
  int k, end = I[i+1];

#ifdef DEBUG
  if ( i>=size || i<0 )
    cerr << "Table::Push() i = " << i << "\n";
#endif

  for(k=I[i]; k<end; k++)
    if (J[k]==j)
      return k;
    else if (J[k]==-1){
      J[k] = j; 
      return k;
    }

// #ifdef DEBUG
  cerr << "Unsuccessful try to push connection (i,j) = (" 
       << i << ", " << j << ")\n";
// #endif

  return -1;
}

//============================================================================

void Table::Finalize(){
  int i, j, end, sum = 0, n = 0, newI = 0;

  for(i=0; i<I[size]; i++)
    if (J[i] != -1)
      sum++;
  
  if (sum != I[size]){
    int *NewJ = new int[sum];
    
    for(i=0; i<size; i++){
      end = I[i+1];
      for(j=I[i]; j<end; j++){
	if (J[j] == -1) break;
	NewJ[ n++ ] = J[j];
      }
      I[i] = newI;
      newI = n;
    }
    I[size] = sum;

    delete [] J;
  
    J = NewJ;

#ifdef DEBUG
  if (sum != n)
    cerr << "Table::Finalize\n";
#endif
  }
}

//============================================================================

int Table::Width() const {
  int width = 0;
  for (int k = 0; k < I[size]; k++)
    if (J[k] > width) width = J[k];
  return width + 1;  
}

//============================================================================

void Table::Print(ostream & out, int width){
  int i, j;

  out << setiosflags(ios::scientific | ios::showpos);
  for(i = 0; i < size; i++) {
      out << "[row " << i << "]\n";
      for (j = I[i]; j < I[i+1]; j++) 
	{
	  out << setw(3) << J[j] << "  ";
	  if ( !((j+1-I[i]) % width) )
	    out << endl;
	}
      out << endl;
    }
  out << endl;
}

//============================================================================

void Table::Save(ostream & out)
{
   int i;

   out << size << '\n';

   for (i = 0; i <= size; i++)
      out << I[i] << '\n';
   for (i = 0; i < I[size]; i++)
      out << J[i] << '\n';
}

//============================================================================

Table::~Table (){
  if (I) delete [] I;
  if (J) delete [] J;
}

//============================================================================

STable::STable (int dim, int connections_per_row) : 
  Table(dim, connections_per_row) {
}

//=============================================================================

STable::STable (int dim, const int *connections_per_row) :
  Table(dim, connections_per_row) {
}

//=============================================================================

int STable::operator() (int i, int j) const
{
  if (i < j)
    return Table::operator()(i,j);
  else
    return Table::operator()(j,i);
}

//=============================================================================

int STable::Push( int i, int j ){
  if (i < j)
    return Table::Push(i, j);
  else
    return Table::Push(j, i);
}

//=============================================================================
