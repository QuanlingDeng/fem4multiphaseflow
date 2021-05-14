#ifndef FILE_ThreeDNURBS
#define FILE_ThreeDNURBS
#include "mesh_header.h"

class ThreeDNURBS : public Element
{
 protected:
  int indices[4];
  
 public:
 ThreeDNURBS() : Element() {};

  ThreeDNURBS(const int *ind, int bg);
  virtual void SetVertices(const int *ind);

  /// Get the element's geometry
  virtual int GetGeometryType() const { return Element::ThreeDNURBS; }

  virtual void GetVertices(Array<int> &v) const;
  virtual int GetNVertices() const {return 4; };

  virtual void GetEdges(Array<int> &e) const { ; }
  
  virtual ~ThreeDNURBS();
};

#endif
