#ifndef FILE_SPARSEMAT
#define FILE_SPARSEMAT

/*
   Data types for sparse matrix, sparse matrix inverse

   Author: Bacuta, Ginting, Tomov
   Date:   01/31/2001
*/

using namespace std;

/// Data type sparse matrix 
class SparseMatrix : public Matrix 
{
private:
  
  /** Arrays for the connectivity information in the CSR storage.
      I is of size "size+1", J is of size the number of nonzero entries
      in the Sparse matrix (actually stored I[size]) **/
  int *I, *J, width;

  /// The nonzero entries in the Sparse matrix with size I[size].
  double *A;

public:
  /// Creates sparse matrix with fixed number of elements per row.
  SparseMatrix (int dim, int elements_per_row = 3, int is_rect = 0);

  /// Creates sparse matrix with variable number of elements per row.
  SparseMatrix (int dim, const int *elements_per_row, int is_rect = 0);

  /// Copy Constructor. Assume finalize has been called.
  SparseMatrix (const SparseMatrix *mat, double c=1.0);

  int Width() { return width; };

  int *GetI() const { return I; }
  int *GetJ() const { return J; }
  double *GetValue() const { return A; }

  /// Returns reference to a_{ij}.  Index i, j = 0 .. size-1 
  virtual double& Elem (int i, int j);

  /// Returns constant reference to a_{ij}.  Index i, j = 0 .. size-1 
  virtual const double& Elem (int i, int j) const;

  /// Returns reference to A[i][j].  Index i, j = 0 .. size-1 
  double& operator() (int i, int j);

  /// Returns reference to A[i][j].  Index i, j = 0 .. size-1 
  const double& operator() (int i, int j) const;

  /// Matrix vector multiplication.   
  virtual void Mult (const Vector & x, Vector & y) const;

  /// y += A * x
  void AddMult (const Vector & x, Vector & y) const;

  /// Compute Infinity Norm.
  double InfinityNorm();

  /// Multiply a vector with the transposed matrix. y = At * x
  void MultTranspose (const Vector & x, Vector & y) const;

  /// y += At * x
  void AddMultTranspose (const Vector & x, Vector & y) const;

  /// Returns a pointer to approximation of the matrix inverse.
  virtual MatrixInverse * Inverse() const;

  /// Eliminates a column from the transpose matrix.
  void EliminateRow (int row, const double sol, Vector &rhs);

  /** Eliminates the column 'rc' to the 'rhs', deletes the row 'rc' and
      replaces the element (rc,rc) with 1.0; assumes that element (i,rc)
      is assembled if and only if the element (rc,i) is assembled.         */
  void EliminateRowCol (int rc, const double sol, Vector &rhs);

  /** Sets matrix entry (rc,rc) to one and 'rhs(rc)' to sol  */
  void ClearRow (int rc, const double sol, Vector &rhs);

  /// Gauss-Seidel forward and backward iterations over a vector x.
  void Gauss_Seidel_forw(const Vector &x, Vector &y) const;
  void Gauss_Seidel_back(const Vector &x, Vector &y) const;

  /** Finalize the matrix initialization. The function should be called
      only once, after the matrix has been initialized. It's densenning
      the J and A arrays (by getting rid of -1s in array J).              **/
  virtual void Finalize();

  void GetSubMatrix(Array<int> &rows, Array<int> &cols,
                    RectangularMatrix &subm);

  void SetSubMatrix(Array<int> &rows, Array<int> &cols,
                    RectangularMatrix &subm);

  void AddSubMatrix(const Array<int> &rows, const Array<int> &cols,
                    const RectangularMatrix &subm);

  /// Add (element) submatrix to the matrix.
  virtual void AddElementMatrix (const Array<int> & dofs,
                                 const Matrix & elemmat)
         { Matrix::AddElementMatrix(dofs, elemmat); };

  /// Add (element) submatrix to the matrix.
  virtual void AddElementMatrix (const Array<int> & dofs1,
                                 const Array<int> & dofs2,
                                 const RectangularMatrix & elemmat);

  /// Prints matrix to stream out.
  void Print(ostream & out = cout, int width = 4);


  void MPrint(ostream & out);

  /// Prints matrix to stream out in hypre_CSRMatrix format.
  void PrintCSR(ostream & out);

  /// Destroys sparse matrix.  
  virtual ~SparseMatrix();  
};

#endif
