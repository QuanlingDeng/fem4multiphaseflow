#ifndef FILE_RECTMAT
#define FILE_RECTMAT

/*
  Data types rectangular matrix, inverse dense matrix

  Author: Kristi
  Date:   01/24/2001
*/

using namespace std;

/// Data type rectangular matrix
class RectangularMatrix : public Matrix
{
private:
  int height;
  double *data;

  friend class DenseMatrixInverse;
  friend void Mult (const RectangularMatrix & b, 
		    const RectangularMatrix & c, 
		    RectangularMatrix & a);

public:
  /** Default constructor for RectangularMatrix.
      Sets data = NULL size = height = 0 */
  RectangularMatrix();

  /// Creates square matrix of size s.
  RectangularMatrix (int s);

  /// Creates rectangular matrix of size n and height m.
  RectangularMatrix (int m, int n);
 
  /// Creates rectangular matrix equal to the transpose of mat.
  RectangularMatrix (const RectangularMatrix & mat, char ch);

  /// If the matrix is not a square matrix of size s then recreate it
  void SetSize(int s);

  /// If the matrix is not a matrix of size (h x w) then recreate it
  void SetSize(int h, int w);

  /// Returns the height of the matrix.
  inline int Height() const {return height;};

  /// Returns reference to a_{ij}.  Index i, j = 0 .. size-1 
  inline double& operator() (int i, int j);

  /// Returns constant reference to a_{ij}.  Index i, j = 0 .. size-1 
  inline const double& operator() (int i, int j) const;

  /// Returns reference to a_{ij}.  Index i, j = 0 .. size-1 
  virtual double& Elem (int i, int j);

  /// Returns constant reference to a_{ij}.  Index i, j = 0 .. size-1 
  virtual const double& Elem (int i, int j) const;

  /// Matrix vector multiplication.
  void Mult (const Vector & x, Vector & y) const;

  /// Multiply a vector with the transpose matrix.
  void MultTranspose (const Vector & x, Vector & y) const;

  /// Returns a pointer to the inverse matrix.
  virtual MatrixInverse * Inverse() const;

  /// Replaces the current matrix with its inverse
  void Invert();

  /// Calculates the determinant of the matrix (for 2x2 or 3x3 matrices)
  double Det();

  /// Adds the matrix A multiplied by the number c to the matrix
  void Add(const double c, RectangularMatrix &A);

  /// Sets the matrix elements equal to constant c
  RectangularMatrix &operator=(double c);

  RectangularMatrix &operator+=(RectangularMatrix &m);

  RectangularMatrix &operator-=(RectangularMatrix &m);

  /// (*this) = -(*this)
  void Neg();

  void Print (ostream & out=cout, int width=4);

  /// Destroys rectangular matrix.
  virtual ~RectangularMatrix();
};


/// Matrix matrix multiplication.  A = B * C.
void Mult (const RectangularMatrix & b, 
	   const RectangularMatrix & c, 
	   RectangularMatrix & a);

/// Calculate the inverse of a matrix (for 2x2 or 3x3 matrices)
void CalcInverse( RectangularMatrix & a, RectangularMatrix & inva );

/// Calculate the matrix A.At
void MultAAt( RectangularMatrix & a, RectangularMatrix & aat);

/// Multiply a matrix A with the transpose of a matrix B:   A*Bt
void MultABt( RectangularMatrix & A, RectangularMatrix & B,
              RectangularMatrix & ABt );

/// Multiply the transpose of a matrix A with a matrix B:   At*B
void MultAtB( RectangularMatrix & A, RectangularMatrix & B,
              RectangularMatrix & AtB );

/// Make a matrix from a vector V.Vt
void MultVVt( Vector & v, RectangularMatrix & vvt);

void MultVWt(Vector &v, Vector &w, RectangularMatrix &VWt);


/** Data type for inverse of square dense matrix.
    Stores LU factors */
class DenseMatrixInverse : public MatrixInverse
{
private:
  double *data;

public:
  /// Creates square dense matrix. Computes factorization of mat and stores LU factors. 
  DenseMatrixInverse (const RectangularMatrix & mat);

  /// Matrix vector multiplication with inverse of dense matrix.
  virtual void Mult (const Vector & x, Vector & y) const;

  /// Destroys dense inverse matrix.
  virtual ~DenseMatrixInverse ();
};




//----------------------------------------------------------------



inline double & RectangularMatrix::operator() (int i, int j)
{
#ifdef DEBUG
  if ( data == 0 || i < 0 || i >= height || j < 0 || j >= size )
    {
      cerr << "RectangularMatrix::operator()\n";
      return data[0];
    }
#endif
  
  return data[i+j*height];
}

inline const double & RectangularMatrix::operator() (int i, int j) const
{
#ifdef DEBUG
  if ( data == 0 || i < 0 || i >= height || j < 0 || j >= size )
    {
      cerr << "RectangularMatrix::operator()\n";
      return data[0];
    }
#endif
  
  return data[i+j*height];
}

#endif
