#ifndef FILE_LINEARFEM3D
#define FILE_LINEARFEM3D
/*
Data type LinearFEM3D
Author: R. Johnson and B. McCaskill
*/

/// LinearFEM1D data type. 
class LinearFEM3D : public FEM
{
protected:
                   
  virtual void ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat) { ; }

  virtual void ComputeConservativeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat) { ; }

  virtual void ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs);


  virtual void GetElementDOF(int i, Array<int> &ind);

  virtual void GetBdrElementDOF(int i, Array<int> &ind);


  virtual void ComputeNeumannLocalSystem(int i, Array<int> &ind, Array<double> &locrhs);

  virtual void ComputeRobinLocalSystem(int i, Array<int> &ind, Array<double *> &locmat, 
				       Array<double> &locrhs);


  virtual void ComputeStableLocalSystem(int i, Array<int> &ind, Array<double *> &locmat) { ; }

  virtual void ComputeStableRHSLocalSystem(int i, Array<int> &ind, Array<double> &locrhs) { ; }

public:

  LinearFEM3D(Mesh *_mesh, Data *_data);

  virtual void SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
				    const Array<double *> &Dval);

  virtual void ProjectAFunction(double (*func)(Array<double> &), Vector &pval);

  virtual double ComputeL2Error(double (*func)(Array<double> &), Vector &approxsol){ return -1.0; }

  virtual double ComputeRobinLocalSystem(int i, Array<int> &ind, 
				       Array<double> &locrhs) { return -1.0; }

};

#endif
