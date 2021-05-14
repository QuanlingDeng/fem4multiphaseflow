#ifndef FILE_SEGMENT
#define FILE_SEGMENT
#include "mesh_header.h"
/*
  Data type line segment element
*/

/// Data type line segment element
class Segment : public Element
{
protected:
  int indices[2];
  int edge;

public:

 Segment() : Element() {}

  /// Constructs a segment by specifying the index of vertices and the attribute.
  Segment(const int *ind, int attr);
  Segment(const int *ind, int e, int attr);

  /// Set the index of vertices of the element according to the input.
  virtual void SetVertices(const int *ind);

  /// Get the element's geometry
  virtual int GetGeometryType() const { return Element::SEGMENT; }

  /// Returns the index of the element's  vertices.
  virtual void GetVertices(Array<int> &v) const;

  /// Returns number of vertices of the element
  virtual int GetNVertices() const { return 2; };

  virtual void GetEdges(Array<int> &e) const;

  virtual ~Segment();
};

#endif
