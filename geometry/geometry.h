#ifndef FILE_GEOMETRY
#define FILE_GEOMETRY
#include "geometry_header.h"

class Geometry {
 protected:
  Array<Patch *> patches;
  Array<Vector *> faces;
  Array<Vector *> connect;

  int NumOfBdrs;
  int NumOfPatches;
  int NumOfFaces;
  int SpaceDim;
  int PatchDim;
  
 public:
  inline int GetNumOfPatches() const { return NumOfPatches;};

  inline int GetSpaceDim() const { return SpaceDim; };

  inline int GetPatchDim() const { return PatchDim; };

  inline int GetPatchDegs(int i, int j) const { return patches[i]->GetDegs(j); };

  inline double GetCtrlPt(int i, int j, int k) const { return patches[i]->GetCtrlPt(j,k); };

  inline double GetWeight(int i, int j) const {return patches[i]->GetWeight(j); };

  inline double GetKnot(int i, int j, int k) const { return patches[i]->GetKnot(j,k); };

  inline int GetKnotSize(int i, int j) const { return patches[i]->GetKnotSize(j); } ;

  inline int GetFaceVal(int i, int j) const { return faces[i]->Elem(j); } ;

  inline int GetNumOfFaces() const { return NumOfFaces; } ;

  inline int GetFaceSize() const { return faces[0]->Size(); } ;

  inline int GetConnectVal(int i, int j) const { return connect[i]->Elem(j); } ;

  Geometry(ifstream *files);

  void PrintNURBS(int *n);

  void PrintBSSurface(int *n);

  //refines solution mesh while preserving the NURBS map 
  void Refine(Array<int *> &factor);

  void alphafill(int w, int i, int k, Vector *knots, double newknot, Vector *alpha);
  
  virtual ~Geometry();
};
#endif
