#ifndef FILE_LINEARFEM3D_P
#define FILE_LINEARFEM3D_P
/*
Data type LinearFEM3D_p
Author: R. Johnson and B. McCaskill
*/

/// LinearFEM1D data type. 
class LinearFEM3D_p : public FEM_p
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

  LinearFEM3D_p(Mesh_p *_mesh, Data *_data);
};

#endif
