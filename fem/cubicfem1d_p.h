#ifndef FILE_CUBICFEM1D_P
#define FILE_CUBICFEM1D_P

/* 
   Data type cubicFEM1D
Author:  R. Johnson, B. McCaskill
*/

//CubicFEM1D data type
class CubicFEM1D_p : public FEM_p
{
 protected:

  virtual void ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat);

  virtual void ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs);

  virtual void GetElementDOF(int i, Array<int> &ind);

  virtual void GetBdrElementDOF(int i, Array<int> &ind) {;} 

 public:

  CubicFEM1D_p(Mesh_p *_mesh, Data *_data);

  virtual void SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
				    const Array<double *> &Dval);

  virtual void ProjectAFunction(double (*func)(Array<double> &), Vector &pval) { ; }


};

#endif
 



