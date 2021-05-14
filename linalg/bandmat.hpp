#ifndef FILE_BANDMAT
#define FILE_BANDMAT

using namespace std;

/// Data type band matrix 
class BandMatrix : public Matrix
{
private:
  int hbw; // Half-band-width
  double *data;

  friend class BandMatrixInverse;

public:
  /// Creates band matrix of size s and halfbandwidth h.
  BandMatrix (int s, int h);

  /// Returns the half-band-width of the band matrix.
  int HalfBandWidth() const { return hbw; };

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

  /// Destroys band matrix.
  virtual ~BandMatrix();
};


/** Data type for inverse of band matrix.
    Stores LU factors */
class BandMatrixInverse : public MatrixInverse
{
private:
  int hbw; // Half-band-width
  double *data;

public:
  /** Creates band inverse matrix.
      Computes the LU factorization of mat and stores it. */
  BandMatrixInverse (const BandMatrix & mat);

  /// Returns the half-band-width of the original band matrix.
  int HalfBandWidth() const { return hbw; };

  /// Matrix vector multiplication with band-inverse Matrix.
  virtual void Mult (const Vector & x, Vector & y) const;

  /// Destroys band inverse matrix.
  virtual ~BandMatrixInverse ();
};


//--------------------------------------------
//   Implementation of the inline functions
//--------------------------------------------

// Returns reference to a_{ij}.  Index i, j = 0 .. size-1
inline double& BandMatrix::operator()(int i, int j)
{
#ifdef DEBUG
   if ( i < 0 || i >= size || j < 0 || j >= size || abs(i-j) > hbw )
   {
      cerr << "BandMatrix::Elem";
      // make the compiler happy
      return data[0];
   }
#endif

   if ( i < j )
      return data[(hbw+j-i)*size+i];
   else
      return data[(i-j)*size+j];
}

// Returns constant reference to a_{ij}.  Index i, j = 0 .. size-1
inline const double& BandMatrix::operator() (int i, int j) const
{
#ifdef DEBUG
   if ( i < 0 || i >= size || j < 0 || j >= size || abs(i-j) > hbw )
   {
      cerr << "BandMatrix::Elem";
      // make the compiler happy
      return data[0];
   }
#endif

   if ( i < j )
      return data[(hbw+j-i)*size+i];
   else
      return data[(i-j)*size+j];
}

#endif
