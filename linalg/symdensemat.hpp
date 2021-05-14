#ifndef FILE_SYMDENSEMAT
#define FILE_SYMDENSEMAT

/*
   Data types symmetric dense matrix, invesre symmetric dense matrix

   Authors: Joanna, Taejong, Veselin
   Date:   01/29/2001
*/

using namespace std;

/// Data type symmetric dense matrix
class SymDenseMatrix : public Matrix
{
private:
  double *data;

  friend class SymDenseMatrixInverse; // Need to declare this
                                      // class above.

public:
  /// Creates symmetric dense matrix of size s.
  SymDenseMatrix (int s);

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

  /// Returns a pointer to (approximation) of the matrix inverse.
  virtual MatrixInverse * Inverse() const;

  /// Destroys symmetric dense matrix.
  virtual ~SymDenseMatrix();
};


/** Data type for inverse of symmetric dense matrix.
    Stores L and D factors of the L D Lt decomposition */
class SymDenseMatrixInverse : public MatrixInverse
{
private:
   double *data;

public:
  /** Creates symmetric dense inverse matrix.
      Computes the L D Lt factorization of mat and stores it. */
  SymDenseMatrixInverse (const SymDenseMatrix & mat);
  
   /** Matrix vector multiplication with the inverse of symmetric
      dense matrix. */
  virtual void Mult (const Vector & x, Vector & y) const;

  /// Destroys inverse symmetric dense matrix.
  virtual ~SymDenseMatrixInverse ();
};


//--------------------------------------------
//   Implementation of the inline functions
//--------------------------------------------

// Returns reference to a_{ij}.  Index i, j = 0 .. size-1
inline double& SymDenseMatrix::operator()(int i, int j)
{
   int k;

#ifdef DEBUG
   if ( i < 0 || i >= size || j < 0 || j >= size )
   {
      cerr << "SymDenseMatrix::Elem";
      // make the compiler happy
      return data[0];
   }
#endif

   if ( i < j )
      k = j-i, j = i;
   else
      k = i-j;
   return data[k*(2*size-k+1)/2+j];
}

// Returns constant reference to a_{ij}.  Index i, j = 0 .. size-1
inline const double& SymDenseMatrix::operator() (int i, int j) const
{
   int k;

#ifdef DEBUG
   if ( i < 0 || i >= size || j < 0 || j >= size )
   {
      cerr << "SymDenseMatrix::Elem";
      // make the compiler happy
      return data[0];
   }
#endif

   if ( i < j )
      k = j-i, j = i;
   else
      k = i-j;
   return data[k*(2*size-k+1)/2+j];
}

#endif
