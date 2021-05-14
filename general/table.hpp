#ifndef FILE_TABLE
#define FILE_TABLE

/*
   Data types for Table.

   Author: Bacuta, Ginting, Tomov
   Date:   03/10/2001
*/


/** Data type Table. Table stores the connectivity of elements of TYPE I
    to elements of TYPE II, for example, it may be Element-To-Face
    connectivity table, etc.                                               **/
class Table
{
protected:

  /// size is the number of TYPE I elements.  
  int size;

  /** Arrays for the connectivity information in the CSR storage.
      I is of size "size+1", J is of size the number of connections
      between TYPE I to TYPE II elements (actually stored I[size]).        **/
  int *I, *J;

public:
  /// Creates an empty table
  Table() { size = 0; I = J = NULL; }

  /// Creates table with fixed number of connections.
  Table (int dim, int connections_per_row = 3);

  /// Creates table with variable number of connections.
  Table (int dim, const int *connections_per_row);

  Table (int nrows, int *partitioning);

  /// Next 6 are used together with the default constructor
  void MakeI (int nrows);
  void AddAColumnInRow (int r) { I[r]++; };
  void AddColumnsInRow (int r, int ncol) { I[r] += ncol; };
  void MakeJ();
  void AddConnection (int r, int c) { J[I[r]++] = c; };
  void ShiftUpI();

  /// Set the size and the number of connections for the table.
  void SetSize(int dim, int connections_per_row);

  /** Set the rows and the number of all connections for the table.
      Does NOT initialize the whole array I ! (I[0]=0 and I[rows]=nnz only) */
  void SetDims(int rows, int nnz);

  /// Returns the number of TYPE I elements.
  inline int Size() const { return size; }

  /** Returns the number of connections in the table. If Finalize() is
      not called, it returns the number of possible connections established
      by the used constructor. Otherwise, it is exactly the number of
      established connections before calling Finalize().                   **/
  inline int Size_of_connections() const { return I[size]; } 

  /** Returns index of the connection between element i of TYPE I and
      element j of TYPE II. If there is no connection between element i 
      and element j established in the table, then the return value is -1. **/
  int operator() (int i, int j) const;

  /// Return row i in array row.
  void GetRow(int i, Array<int> &row) const;

  /// Set row 'i' as in 'row'
  void SetRow(int i, Array<int> &row);

  int RowSize (int i) { return I[i+1]-I[i]; };

  int *GetRow (int i) { return J+I[i]; };

  int *GetI () { return I; };
  int *GetJ () { return J; };

  void MakeTransposeOf (int *ii, int *jj, int nrows, int _ncols = -1);

  void MakeTransposeOf (Table &A, int _ncols = -1)
         { MakeTransposeOf (A.GetI(), A.GetJ(), A.Size(), _ncols); };

  /// Return the transpose table in 'tr'. Must be Finalized().
  void GetTranspose (Table &tr, int _ncols = -1);

  void SortColumnIndexes();

  /** Establish connection between element i and element j in the table.
      The return value is the index of the connection. It returns -1 if it
      fails to establish the connection. Possibilities are there is not
      enough memory on row i to establish connection to j, an attempt to
      establish new connection after calling Finalize().                   **/
  int Push( int i, int j );

  /** Finalize the table initialization. The function may be called
      only once, after the table has been initialized, in order to densen
      array J (by getting rid of -1's in array J). Calling this function
      will "freeze" the table and function Push will work no more. 
      Note: The table is functional even without calling Finalize().       **/
  void Finalize();

  /// Returns the number of TYPE II elements (after Finalize() is called).
  int Width() const;

  /// Prints the table to stream out.
  void Print(ostream & out = cout, int width = 4);

  void Save(ostream & out);

  /// Destroys Table.  
  ~Table();  
};



/** Data type STable. STable is similar to Table, but it's for symmetric
    connectivity, i.e. TYPE I is equivalent to TYPE II. In the first 
    dimension we put the elements with smaller index.                      **/

class STable : public Table
{
public:
  /// Creates table with fixed number of connections.
  STable (int dim, int connections_per_row = 3);

  /// Creates table with variable number of connections.
  STable (int dim, const int *connections_per_row);

  /** Returns index of the connection between element i of TYPE I and
      element j of TYPE II. If there is no connection between element i 
      and element j established in the table, then the return value is -1. **/
  int operator() (int i, int j) const;

  /** Establish connection between element i and element j in the table.
      The return value is the index of the connection. It returns -1 if it
      fails to establish the connection. Possibilities are there is not
      enough memory on row i to establish connection to j, an attempt to
      establish new connection after calling Finalize().                   **/
  int Push( int i, int j );

  /// Destroys STable.  
  ~STable() {}  
};

#endif
