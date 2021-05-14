#ifndef FILE_LINEARFEM1D_P
#define FILE_LINEARFEM1D_P
/*
Data type LinearFEM1D
Author: M. Seo and X. Wu
*/

/// LinearFEM1D data type. 
class LinearFEM1D_p : public FEM_p
{
protected:
                   
  virtual void ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs);

  virtual void GetElementDOF(int i, Array<int> &ind);

  virtual void GetBdrElementDOF(int i, Array<int> &ind);

public:

  LinearFEM1D_p(Mesh_p *_mesh, Data *_data);
};

#endif
