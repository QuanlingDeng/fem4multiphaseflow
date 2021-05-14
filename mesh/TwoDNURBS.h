#ifndef FILE_TwoDNURBS
#define FILE_TwoDNURBS
#include "mesh_header.h"

class TwoDNURBS : public Element
{
 protected:
  int indices[4];
  
 public:
 TwoDNURBS() : Element() {};

  TwoDNURBS(const int *ind, int bg);
  virtual void SetVertices(const int *ind);

  /// Get the element's geometry
  virtual int GetGeometryType() const { return Element::TwoDNURBS; }

  virtual void GetVertices(Array<int> &v) const;
  virtual int GetNVertices() const {return 4; };

  virtual void GetEdges(Array<int> &e) const { ; }
  
  virtual ~TwoDNURBS();
};

#endif
