#ifndef FILE_PATCH
#define FILE_PATCH
#include "geometry_header.h"

class Patch {
 protected:
  Array<Container *> patch;
  Vector *face;
  
 public:

  Patch(ifstream *file, int pdim, int sdim);

  inline double GetWeight(int i) const { return patch[2]->GetVal(0,i); };

  inline int GetDegs(int i) const { return patch[3]->GetVal(i,0); };

  inline double  GetCtrlPt(int i, int j) const { return patch[1]->GetVal(i,j); };

  inline double  GetKnot(int i, int j) const { return patch[0]->GetVal(i,j); };

  inline Vector*  GetKnot(int i) const { return patch[0]->ReturnData(i); };

  inline int GetKnotSize(int i) const{ return patch[0]->GetSize(i); };

  inline int GetCtrlPtSize() const { return patch[1]->GetSize(0); }; // all same number of control points in x_i for all i  

  inline double GetFaceVal(int i) const {return face->Elem(i);};

  inline Vector* GetFaceVec() const {return face;};

  inline int GetFaceSize() const {return face->Size();};

  void Refine(Array<Vector *> &knots, Array<Vector *> &ctrlpts);

  virtual ~Patch();
};
#endif
