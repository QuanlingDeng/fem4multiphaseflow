#ifndef FILE_LINEARFVEM2D
#define FILE_LINEARFVEM2D
/*
Data type LinearFVEM2D
Author: Q. Deng
Date: 06/09/2015
*/

/// LinearFVEM2D data type. 
class LinearFVEM2D : public FEM
{
protected:
                   
  virtual void ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeConservativeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

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

  LinearFVEM2D(Mesh *_mesh, Data *_data);

  virtual void SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
				    const Array<double *> &Dval);

  virtual void ProjectAFunction(double (*func)(Array<double> &), Vector &pval);

};

#endif
