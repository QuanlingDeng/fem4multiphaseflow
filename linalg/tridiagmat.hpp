#ifndef FILE_TRIDIAGMAT
#define FILE_TRIDIAGMAT

/*
   Data types tridiagonal matrix, inverse tridiagonal matrix

   Author: Tzanio, Joachim
   Date:   01/18/2001
*/


using namespace std;

/// Data type tridiagonal matrix 
class TriDiagonalMatrix : public Matrix
{
private:
  double *l, *d, *u;

  friend class TriDiagonalMatrixInverse;

public:
  /// Creates tridiagonal matrix of size s. 
  TriDiagonalMatrix (int s);

  /// Returns reference to a_{ij}.  Index i, j = 0 .. size-1 
  inline double& operator()(int i, int j);

  /// Returns constant reference to a_{ij}.  Index i, j = 0 .. size-1 
  inline const double& operator() (int i, int j) const;

  /// Returns reference to a_{ij}.  Index i, j = 0 .. size-1 
  virtual double& Elem (int i, int j);

  /// Returns constant reference to a_{ij}.  Index i, j = 0 .. size-1 
  virtual const double& Elem (int i, int j) const;

  /// Matrix vector multiplication.   
  virtual void Mult (const Vector & x, Vector & y) const;

  /// Returns a pointer to the inverse matrix.
  virtual MatrixInverse * Inverse() const;

  /// Destroys tridiagonal matrix.  
  virtual ~TriDiagonalMatrix();  
};


/** Data type for inverse of tridiagonal matrix.
    Stores LU factors */
class TriDiagonalMatrixInverse : public MatrixInverse
{
private:
  double *l, *d, *u;

public:
  /// Creates tridiagonal inverse matrix. Computes factorization of mat and stores LU factors. 
  TriDiagonalMatrixInverse (const TriDiagonalMatrix & mat);
  
  /// Matrix vector multiplication with TD-inverse Matrix. 
  virtual void Mult (const Vector & x, Vector & y) const;

  /// Makes LU factorization of a new TriDiagonalMatrix
  void Rebuild(const TriDiagonalMatrix & mat);

  /// Destroys tridiagonal inverse matrix.  
  virtual ~TriDiagonalMatrixInverse ();
};



//----------------------------------------------------------------



inline double & TriDiagonalMatrix::operator()(int i, int j)
{
#ifdef DEBUG
  if ( i < 0 || i >= size || j < 0 || j >= size )
    {
      cerr << "TriDiagonalMatrix::operator()";
      // make the compiler happy
      return d[0];
    }
#endif

  switch (i-j) 
    {
    case  0: return d[i];
    case -1: return u[i];
    case  1: return l[j];
    default:
#ifdef DEBUG
      cerr << "TriDiagonalMatrix::operator()"; 
#endif
      return d[0];
    }
}

inline const double & TriDiagonalMatrix::operator()(int i, int j) const
{
#ifdef DEBUG
  if ( i < 0 || i >= size || j < 0 || j >= size )
    {
      cerr << "TriDiagonalMatrix::operator()";
      // make the compiler happy
      return d[0];
    }
#endif

  switch (i-j) 
    {
    case  0: return d[i];
    case -1: return u[i];
    case  1: return l[j];
    default:
#ifdef DEBUG
      cerr << "TriDiagonalMatrix::operator()"; 
#endif
      return d[0];
    }
}
	
#endif
