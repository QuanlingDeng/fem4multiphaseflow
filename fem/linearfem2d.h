#ifndef FILE_LINEARFEM2D
#define FILE_LINEARFEM2D
/*
Data type LinearFEM2D
Author: Q. Deng and P. Torsu
*/

/// LinearFEM2D data type. 
class LinearFEM2D : public FEM
{
protected:
                   
  virtual void ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeConservativeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs);


  virtual void GetElementDOF(int i, Array<int> &ind);

  virtual void GetBdrElementDOF(int i, Array<int> &ind);


  virtual void ComputeNeumannLocalSystem(int i, Array<int> &ind, Array<double> &locrhs);

  virtual void ComputeNeumannLocalSystem(int b, int i, Array<int> &ind, Array<double> &locrhs);

  virtual void ComputeRobinLocalSystem(int i, Array<int> &ind, Array<double *> &locmat, 
				       Array<double> &locrhs);

  virtual void ComputeRobinLocalSystem(int i, int b, Array<int> &ind, Array<double *> &locmat, 
				       Array<double> &locrhs);

  virtual void ComputeStableLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeStableRHSLocalSystem(int i, Array<int> &ind, Array<double> &locrhs);

public:

  LinearFEM2D(Mesh *_mesh, Data *_data);

  virtual void SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
				    const Array<double *> &Dval);

  virtual void ProjectAFunction(double (*func)(Array<double> &), Vector &pval);

  virtual double ComputeL2Error(double (*solfunc)(Array<double> &), Vector &approxsol);

  virtual double ComputeH1Error(double (*deronefunc)(Array<double> &), double (*dertwofunc)(Array<double> &), Vector &approxsol);

  //Compute the L2 error between approximate and exact solutions
  virtual double ComputeL2Error(Function *exactsol, Vector &femsol, double time=1.0);

  virtual double ComputeH1Error(Array<Function *> &exactderivative, Vector &femsol);

  virtual void ComputePPL2Error(Array<double> &errs, Function *exactsol, Vector &femsol, Array<double *> &ppsol, double time=1.0);

  virtual void ComputePPH1Error(Array<double> &errs, Array<Function *> &exactderivative, Vector &femsol, Array<double *> &ppsol);

};

#endif
