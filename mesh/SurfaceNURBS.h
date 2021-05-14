#ifndef FILE_SurfaceNURBS
#define FILE_SurfaceNURBS
#include "mesh_header.h"

class SurfaceNURBS : public Element
{
 protected:
  int indices[4];
  
 public:
 SurfaceNURBS() : Element() {};

  SurfaceNURBS(const int *ind, int bg);
  virtual void SetVertices(const int *ind);

  /// Get the element's geometry
  virtual int GetGeometryType() const { return Element::SurfaceNURBS; }

  virtual void GetVertices(Array<int> &v) const;
  virtual int GetNVertices() const {return 4; };

  virtual void GetEdges(Array<int> &e) const { ; }
  
  virtual ~SurfaceNURBS();
};

#endif
