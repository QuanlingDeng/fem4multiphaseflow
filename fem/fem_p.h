#ifndef FILE_FEM_P
#define FILE_FEM_P

 /*
   Data type FEM

   Author: V. Ginting
*/

/// FEM data type. This is a base class from which various classes are derived. 
class FEM_p
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
  Mesh_p *mesh;

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

  virtual void ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat) = 0;

  virtual void ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs) = 0;

  /* Compute local Neumann system */
  virtual void ComputeNeumannLocalSystem(int i, Array<int> &ind, 
  					 Array<double> &locrhs) = 0;

public:

  // Constructor
  FEM_p(Mesh_p *_mesh, Data *_data);

  // Destructor
  virtual ~FEM_p();

  // A function to create and fill mat and rhs.
  void Assemble();

  // Solve the linear system.
  void Solve(Vector &sol, int maxiter=1000, int printopt=0);

  // Include Dirichlet condition into the global system
  void SetDirichletBoundary(Array<bool> &bdry_is_dirichlet, const Vector &Dval);

  inline int GetNumGlobalDOF() { return NumGlobalDOF; };

  inline int GetNumLocalDOF() { return NumLocalDOF; };

  inline char *GetType() { return Type; };

  Element *GetElement(int i) const { return mesh->GetElement(i); }

  Element *GetBdrElement(int i) const { return mesh->GetBdrElement(i); }

  inline int GetNBE() const { return mesh->GetNBE(); }

  /// Given a function func(Array<double> &), get the associated projection to the finite element space.
  virtual void ProjectAFunction(double (*func)(Array<double> &), Vector &pval) = 0;

};

#endif
