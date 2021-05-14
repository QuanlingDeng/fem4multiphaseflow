//==========================================================================
// The class "multifrontal" defined here provides an interface to UMFPACK2.2
//==========================================================================
#ifndef FILE_MULTIFRONTAL
#define FILE_MULTIFRONTAL

using namespace std;

class multifrontal
{
 private:
  // input parameters (these are not changed even after factorization)
  int N;      // Matrix is N x N
  int Nelem;  // Number of nonzero elements in the matrix
  /*                  
		      The k-th nonzero entry (k=0..Nelem-1) in the matrix A,
		      say A[i,j], is defined by
  		         A[Index[k],Index[k+Nelem]] = Value[k], 
		      where "Index" and "Value" are arrays declared below. 
  */
  int *Index;      // Nonzero entry indices: Array of length 2*Nelem
  double  *Value;  // The matrix: an array of length Nelem
  
  // Variables needed for factorization and solution calls
  int n;
  int nelem;
  int *index;
  int lindex;
  double  *value;
  int lvalue; 
  double  *work; 
  int job;
  int transa;
  int keep[20];
  double cntl[10];
  int icntl[20];
  int info[40];
  double rinfo[20];

  int size_factor; // number of doubles given as scratch space for 
                   // factorization will be = size_factor * Nelem.

 public:
  // Creates sparse matrix with dimension dim and number of nonzero
  // elements nelem, ind is the array containing nonzero entries 
  // indices, val is the array of nonzero entries, mem_mult = size_factor
  multifrontal( int dim, int nel, int *ind, double *val, int mem_mult );

  // LU factorization, calls umd2fa subroutine from the package
  int LU_fact();

  // Solves the linear system, must be done after factorization,
  // calls umd2so subroutine from the package
  void Solve( const double *b, double *x );

  // Destroys the multifrontal
  ~multifrontal();
 
};

#endif
