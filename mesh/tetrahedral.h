#ifndef FILE_TETRAHEDRAL
#define FILE_TETRAHEDRAL
#include "mesh_header.h"

class Tetrahedral : public Element
{
 protected:
  int indices[4];
  
 public:
 Tetrahedral() : Element() {};

  Tetrahedral(const int *ind, int bg);
  virtual void SetVertices(const int *ind);

  /// Get the element's geometry
  virtual int GetGeometryType() const { return Element::TETRAHEDRAL; }

  virtual void GetVertices(Array<int> &v) const;
  virtual int GetNVertices() const {return 4; };

  virtual void GetEdges(Array<int> &e) const { ; }
  
  virtual ~Tetrahedral();
};

#endif
