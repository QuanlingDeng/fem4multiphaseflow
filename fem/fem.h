#ifndef FILE_FEM
#define FILE_FEM

 /*
   Data type FEM

   Author: V. Ginting
*/

/// FEM data type. This is a base class from which various classes are derived. 
class FEM
{
protected:

  // number of degree of freedoms, i.e., how many unknowns to be solved.
  int NumGlobalDOF;

  // a vector containing the "discretized" representation of linear functional.
  Vector rhs;

  // a matrix containing the "discretized" representation of bilinear form.
  SparseMatrix *mat;

  // contains the nodal values of coefficient in the bilinear form and in the linear functional.
  Data *data;

  /* Flag for type of finite element, i.e., LINEAR, QUADRATIC, BILINEAR, BIQUADRATIC.
     This is set in the constructor of each derived class. */
  char Type[256];

  /* Number of unknowns in an element that is unique for each finite element type and dimension.
     This is set in the constructor of each derived class. */
  int NumLocalDOF;

  // contains the length of each subintervals as a discretization of [a,b].
  Mesh *mesh;

  /* Gives back the global index of unknowns associated with element i. This function is
     unique for each finite element type and is implemented in each derived class. */
  virtual void GetElementDOF(int i, Array<int> &ind) = 0;

  /* Gives back the global index of unknowns associated with boundary element i. This function is
     unique for each finite element type and is implemented in each derived class. */
  virtual void GetBdrElementDOF(int i, Array<int> &ind) = 0;

  // Gives the vertex index of an element i.
  void GetElementVertices(int i, Array<int> &ind);

  // Gives the vertex index of a boundary element i.
  void GetBdrElementVertices(int i, Array<int> &ind);


  /* Functions for calculation of local matrix (locmat) and local vector (locrhs)
     for element i. This calculation is unique for each finite element type and is
     implemented in each derived class. This function is called from Assemble(...). */
  virtual void ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat) = 0;

  virtual void ComputeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat) = 0;

  virtual void ComputeConservativeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat) = 0;

  virtual void ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat) = 0;

  virtual void ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs) = 0;

  /* Compute local Neumann system */
  virtual void ComputeNeumannLocalSystem(int i, Array<int> &ind, 
  					 Array<double> &locrhs) = 0;

  virtual void ComputeNeumannLocalSystem(int b, int i, Array<int> &ind, 
  					 Array<double> &locrhs) = 0;

  /* Compute local Robin system */
  virtual void ComputeRobinLocalSystem(int i, Array<int> &ind, 
  				       Array<double *> &locmat, 
				       Array<double> &locrhs) = 0;

  virtual void ComputeRobinLocalSystem(int i, int b, Array<int> &ind, 
  				       Array<double *> &locmat, 
				       Array<double> &locrhs) = 0;

  /* Stablization for advection dominated problems */
  virtual void ComputeStableLocalSystem(int i, Array<int> &ind, Array<double *> &locmat) = 0;

  virtual void ComputeStableRHSLocalSystem(int i, Array<int> &ind, Array<double> &locrhs) = 0;

public:

  // Constructor
  FEM(Mesh *_mesh, Data *_data);

  // Destructor
  virtual ~FEM();

  // A function to create and fill mat and rhs.
  void Assemble();

  // A function to create and fill mat and rhs.
  void AssembleMatrix();

  // A function to update rhs.
  void UpdateRHS();

  void FinalizeMatrix(){ mat->Finalize();}

  // Solve the linear system.
  //  void Solve(Vector &sol, int maxiter=1000, int printopt=0);
  void Solve(Vector &sol, char *linalgsolv, int maxit=1000, int printit=0, 
	     double rtol=1e-14, double atol=1e-16);

  inline void GetEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat){ ComputeEllipticLocalSystem(i, ind, locmat); }

  inline void GetAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat){ ComputeAdvectionLocalSystem(i, ind, locmat); }

  inline void GetReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat){ ComputeReactionLocalSystem(i, ind, locmat); }

  inline void GetForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs){ ComputeForceLocalSystem(i, ind, locrhs); }

  inline void GetNeumannLocalSystem(int b, int i, Array<int> &ind, Array<double> &locrhs){ ComputeNeumannLocalSystem(b, i, ind, locrhs); }

  inline void GetNeumannLocalSystem(int i, Array<int> &ind, Array<double> &locrhs){ ComputeNeumannLocalSystem(i, ind, locrhs); }

  // Include Dirichlet condition into the global system
  //  void SetDirichletBoundary(Array<bool> &bdry_is_dirichlet, const Vector &Dval);
  virtual void SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
			    const Array<double *> &Dval) = 0;


  inline int GetNumGlobalDOF() { return NumGlobalDOF; };

  inline int GetNumLocalDOF() { return NumLocalDOF; };

  inline char *GetType() { return Type; };

  Element *GetElement(int i) const { return mesh->GetElement(i); }

  Element *GetBdrElement(int i) const { return mesh->GetBdrElement(i); }

  SparseMatrix *getMatrix(){ return mat; }

  inline void getRHS(Vector &_rhs){ _rhs = rhs; }

  inline int GetNBE() const { return mesh->GetNBE(); }

  inline int GetNE() const { return mesh->GetNE(); }

  inline void getElementDOF(int i, Array<int> &ind) { GetElementDOF(i, ind); };

  inline void getBdrElementDOF(int i, Array<int> &ind) { GetBdrElementDOF(i, ind); };

  /// Given a function func(Array<double> &), get the associated projection to the finite element space.
  virtual void ProjectAFunction(double (*func)(Array<double> &), Vector &pval) = 0;
 
  //Compute the L2 error between approximate and exact solutions
  virtual double ComputeL2Error(double (*func)(Array<double> &), Vector &approxsol) = 0;

  virtual double ComputeH1Error(double (*deronefunc)(Array<double> &), double (*dertwofunc)(Array<double> &), Vector &approxsol) = 0;

  //Compute the L2 error between approximate and exact solutions
  virtual double ComputeL2Error(Function *exactsol, Vector &femsol, double time=1.0) = 0;

  virtual double ComputeH1Error(Array<Function *> &exactderivative, Vector &femsol) = 0;

  virtual void ComputePPL2Error(Array<double> &errs, Function *exactsol, Vector &femsol, Array<double *> &ppsol, double time=1.0) = 0;

  virtual void ComputePPH1Error(Array<double> &errs, Array<Function *> &exactderivative, Vector &femsol, Array<double *> &ppsol) = 0;

};

#endif
