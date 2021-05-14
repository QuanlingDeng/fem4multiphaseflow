#ifndef FILE_QUADRATICFEM2D_P
#define FILE_QUADRATICFEM2D_P
/*
Data type QuadraticFEM2D
Author: Q. Deng and P. Torsu
*/

/// QuadraticFEM2D data type. 
class QuadraticFEM2D_p : public FEM_p
{
protected:
                   
  virtual void ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs);

  virtual void GetElementDOF(int i, Array<int> &ind);

  virtual void GetBdrElementDOF(int i, Array<int> &ind);

  virtual void ProjectAFunction(double (*func)(Array<double> &), Vector &pval);

public:

  QuadraticFEM2D_p(Mesh_p *_mesh, Data *_data);
};

#endif
