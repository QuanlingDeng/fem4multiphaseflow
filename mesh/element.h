#ifndef FILE_ELEMENT
#define FILE_ELEMENT

/*
  Abstract data type element
*/

/// Abstract data type element
class Element
{
protected:

  /* Element's attribute (specifying material property, etc).       */
  int attribute;

public:

  /// Constants for the classes derived from Element.
  enum Type { POINT, SEGMENT, TRIANGLE, QUADRILATERAL, TETRAHEDRAL, TwoDNURBS, ThreeDNURBS, SurfaceNURBS };

  Element() { attribute = -1; }

  /// Set the index of the element's vertices according to the input.
  virtual void SetVertices(const int *ind) = 0;

  /// Get the element's geometry
  virtual int GetGeometryType() const = 0;

  /// Returns index of the element's vertices.
  virtual void GetVertices(Array<int> &v) const = 0;

  /// Returns index of the element's edges.
  virtual void GetEdges(Array<int> &e) const = 0;

  /// Return number of vertices of the element.
  virtual int GetNVertices() const = 0;

  /// Return element's attribute.
  inline int GetAttribute() const { return attribute; }

  /// Set element's attribute.
  inline void SetAttribute (const int attr) { attribute = attr; }

  /// Destroys element.
  virtual ~Element() {};
};

#endif
