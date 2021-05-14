#ifndef FILE_VERTEX
#define FILE_VERTEX

/// Data type for vertex
template <int Dim>
class Vertex
{
protected:

  double coord[Dim];

public:

  Vertex(double x);

  Vertex(double x, double y);

  Vertex(double x, double y, double z);

  /// Returns pointer to the coordinates of the vertex.
  inline double * operator() () const { return (double*)coord; };
  
  /// Returns the i'th coordinate of the vertex.
  inline double & operator() (int i) { return coord[i]; };

  /// Returns the i'th coordinate of the vertex.
  inline const double & operator() (int i) const { return coord[i]; };

  /// Get the dimension
  inline int GetDim() const { return Dim; }

  inline void UpdateVertex(int i, double x){ coord[i] = x; }

  ~Vertex();
};

#endif
