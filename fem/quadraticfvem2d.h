#ifndef FILE_QUADRATICFVEM2D
#define FILE_QUADRATICFVEM2D
/*
Data type QuadraticFVEM2D
Author: Q. Deng
Date: 06/12/2015
*/

/// QuadraticFVEM2D data type. 
class QuadraticFVEM2D : public FEM
{
protected:
             
  virtual void ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat) { ; }

  virtual void ComputeConservativeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat) { ; }

  virtual void ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs);

  virtual void GetElementDOF(int i, Array<int> &ind);

  virtual void GetBdrElementDOF(int i, Array<int> &ind);

  virtual void ComputeNeumannLocalSystem(int i, Array<int> &ind, Array<double> &locrhs) {;}

  virtual void ComputeRobinLocalSystem(int i, Array<int> &ind, Array<double *> &locmat, 
				       Array<double> &locrhs) {;}

  virtual void ComputeStableLocalSystem(int i, Array<int> &ind, Array<double *> &locmat) {;}

  virtual void ComputeStableRHSLocalSystem(int i, Array<int> &ind, Array<double> &locrhs) {;}  

public:

  QuadraticFVEM2D(Mesh *_mesh, Data *_data);

  virtual void SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
				    const Array<double *> &Dval);

  virtual void ProjectAFunction(double (*func)(Array<double> &), Vector &pval);
};

#endif
