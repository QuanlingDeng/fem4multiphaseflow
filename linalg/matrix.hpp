#ifndef FILE_MATRIX
#define FILE_MATRIX

/*
  Abstract data types matrix, inverse matrix

  Author: Tzanio
  Date:   01/27/2001
*/
using namespace std;

#include "../general/array.hpp"

class  MatrixInverse;
class  CGMatrixInverse;
class  PCGMatrixInverse;
class  RectangularMatrix;

/// Abstract data type matrix
class Matrix
{
protected:
  int size;

  friend class MatrixInverse;

public:
  /// Creates matrix of width s.
  Matrix (int s) {size=s;};

  /// Returns the width of a matrix.
  inline int Size() const {return size;};

  /// Returns reference to a_{ij}.  Index i, j = 0 .. size-1 
  virtual double& Elem (int i, int j) = 0;

  /// Returns constant reference to a_{ij}.  Index i, j = 0 .. size-1 
  virtual const double& Elem (int i, int j) const = 0;

  /// Matrix vector multiplication.
  virtual void Mult (const Vector & x, Vector & y) const = 0;

  /// Add (element) submatrix to the matrix.
  virtual void AddElementMatrix (const Array<int> & dofs, const Matrix & elemmat);

  /// Add (element) submatrix to the matrix.
  virtual void AddElementMatrix (const Array<int> & dofs1,
                                 const Array<int> & dofs2,
                                 const RectangularMatrix & elemmat);

  /// Returns a pointer to (approximation) of the matrix inverse.
  virtual MatrixInverse * Inverse() const = 0;

  /** Returns a pointer to CGMatrixInverse class. When the default arguments 
      are used CG doesn't print current residual and number of iterations, 
      maximum number of iterations is 1000 and the tolerance is 10e-12. **/  
  CGMatrixInverse * CG(int printiter = 0, int maxnumiter = 1000, 
	               double rtol = 10e-12, double atol = 10e-24) const;

  /** Returns a pointer to PCGMatrixInverse class. When the default arguments 
      are used PCG doesn't print current residual and number of iterations, 
      maximum number of iterations is 1000 and the tolerance is 10e-12. **/  
  PCGMatrixInverse * PCG(const MatrixInverse & b, int printiter = 0, 
                         int maxnumiter = 1000, 
                         double rtol = 10e-12, double atol = 10e-24) const;

  /// Finalizes the matrix initialization.
  virtual void Finalize();

  /// Prints matrix to stream out.
  virtual void Print (ostream & out = cout, int width =4);

  /// Destroys matrix.
  virtual ~Matrix();
};


/// Abstract data type for matrix inverse
class MatrixInverse
{
protected:
  const Matrix * a;
  int size;

public:
  /// Creates approximation of the inverse of square matrix
  MatrixInverse (const Matrix & mat) { size = mat.size; a = &mat;}

  /// Returns the size of the matrix.
  inline int Size() const {return size;}

  /// Matrix vector multiplication with the inverse.
  virtual void Mult (const Vector & x, Vector & y) const = 0;

  /// Destroys inverse matrix.
  virtual ~MatrixInverse();
};

class IdentityInverse : public MatrixInverse
{
   public:
      IdentityInverse(const Matrix &mat) : MatrixInverse (mat) { };

      virtual void Mult (const Vector & x, Vector & y) const { y = x; };
};

#endif
