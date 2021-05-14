#ifndef FILE_SPARSEMATSMOOTHERS
#define FILE_SPARSEMATSMOOTHERS

/*
   Data types for sparse matrix smoothers

   Author: Constantin, Tony, Victor
   Date:   03/04/2001
*/

using namespace std;

/// Data type for Gauss-Seidel smoother of sparse matrix
class GSSmoother : public MatrixInverse
{	
public:

  /// Create GSSmoother.
  GSSmoother (const SparseMatrix & a);

  /// Matrix vector multiplication with GS Smoother.
  virtual void Mult ( const Vector & x, Vector & y) const;

  /// Destroys the GS Smoother.
  virtual ~GSSmoother();
};

//============================================================================


/// Data type for scaled Diagonal smoother of sparse matrix
class DSmoother : public MatrixInverse
{
private:
  /// Scale for the Diagonal smooter. 
  double scale;
  
public:

  /// Create the diagonal smoother.
  DSmoother (const SparseMatrix & a, double scale = 1.);

  /// Matrix vector multiplication with Diagonal smoother.
  virtual void Mult ( const Vector & x, Vector & y) const;

  /// Destroys the Diagonal smoother.
  virtual ~DSmoother();
};

//============================================================================

#endif
