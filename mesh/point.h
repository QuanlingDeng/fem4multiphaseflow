#include "mesh_header.h"

#ifndef FILE_POINT
#define FILE_POINT

/*
  Data type point element
*/

/// Data type point element
class Point : public Element
{
protected:
  int indices[1];

public:

  Point() : Element() {}

  /// Constructs point by specifying the vertex and the attribute.
  Point(const int *ind, int attr);

  /// Set the index of vertex of the element according to the input.
  virtual void SetVertices(const int *ind);

  /// Get the element's geometry
  virtual int GetGeometryType() const { return Element::POINT; }

  /// Returns the index of the element's  vertex.
  virtual void GetVertices( Array<int> &v ) const;

  /// Returns number of vertices.
  virtual int GetNVertices() const { return 1; }

  virtual void GetEdges(Array<int> &e) const { ; }

  virtual ~Point();
};

#endif
