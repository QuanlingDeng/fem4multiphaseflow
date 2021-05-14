#ifndef FILE_LINEARFEM2D_p
#define FILE_LINEARFEM2D_p
/*
Data type LinearFEM2D
Author: Q. Deng and P. Torsu
*/

/// LinearFEM2D data type. 
class LinearFEM2D_p : public FEM_p
{
protected:
                   
  virtual void ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs);

  virtual void GetElementDOF(int i, Array<int> &ind);

  virtual void GetBdrElementDOF(int i, Array<int> &ind);

  virtual void ProjectAFunction(double (*func)(Array<double> &), Vector &pval);

  virtual void ComputeNeumannLocalSystem(int i, Array<int> &ind, Array<double> &locrhs);

public:

  LinearFEM2D_p(Mesh_p *_mesh, Data *_data);
};

#endif
