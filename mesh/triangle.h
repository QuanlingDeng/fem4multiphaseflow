#ifndef FILE_TRIANGLE
#define FILE_TRIANGLE
#include "mesh_header.h"

class Triangle : public Element
{
 protected:
  int indices[3];
  int edges[3];

 public:
  Triangle() : Element(){};
  Triangle(const int *ind, int bg);
  Triangle(const int *ind, int *edge, int bg);

  virtual void SetVertices(const int *ind);

  /// Get the element's geometry
  virtual int GetGeometryType() const { return Element::TRIANGLE; }

  virtual void GetVertices(Array<int> &v) const ;
  virtual int GetNVertices() const {return 3; };

  virtual void GetEdges(Array<int> &e) const;

  virtual ~Triangle();
};

#endif
