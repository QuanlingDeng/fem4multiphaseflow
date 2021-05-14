#ifndef FILE_MESH_p
#define FILE_MESH_p
#include "mpi_mesh_header.h"
/*
  Data type mesh
*/

class Mesh_p {

protected:
  Array<Element *> element;
  Array<Element *> bdrelement;

  int NumOfVertices;
  int NumOfElements;
  int NumOfBdrElements;
  int NumOfBdrs;

public:

  /// Returns number of vertices.
  inline int GetNV() const { return NumOfVertices; }

  /// Returns number of elements.
  inline int GetNE() const { return NumOfElements; }

  /// Returns number of boundary elements.
  inline int GetNBE() const { return NumOfBdrElements; }

  /// Returns number of boundarys.
  inline int GetNBdrs() const { return NumOfBdrs; }

  Element *GetElement(int i) const { return element[i]; }

  Element *GetBdrElement(int i) const { return bdrelement[i]; }

  /// Returns the indices of the vertices of element i.
  void GetElementVertices(int i, Array<int> &ind) const { element[i]->GetVertices(ind); }

  /// Returns the indices of the vertices of boundary element i.
  void GetBdrElementVertices(int i, Array<int> &ind) const { bdrelement[i]->GetVertices(ind); }

  /// Returns the type of element i.
  int GetElementType(int i) const { return element[i]->GetGeometryType(); }

  /// Returns the type of boundary element i.
  int GetBdrElementType(int i) const { return bdrelement[i]->GetGeometryType(); }

  /// Return element i's attribute.
  int GetAttribute(int i) const { return element[i]->GetAttribute(); }; 

  /// Return boundary element i's attribute.
  int GetBdrAttribute(int i) const { return bdrelement[i]->GetAttribute(); };

  int GetElementNVertices(int i) const {return element[i]->GetNVertices(); };

 int GetBdrElementNVertices(int i) const {return bdrelement[i]->GetNVertices(); };

  /// Give the coordinates of every vertex of element i
  virtual void GetElementVerticesCoord(int i, Array<double *> &coord) const = 0;

  /// Give the coordinates of every vertex of boundary element i
  virtual void GetBdrElementVerticesCoord(int i, Array<double *> &coord) const = 0;

  /// Print the mesh information for glvis plot
  void Print(ofstream &out, int type) const;

  /// Print the mesh information for paraview plot
  void ParaviewPrint(int type) const;

  /// Print the solution to a file for glvis plot
  void Print(ofstream &out, const Vector &solution) const;

 /// Print the solution information for paraview plot
  void ParaviewPrint(const Vector &solution, int type) const;

  /// Return the spatial dimension.
  virtual int GetDim() const = 0;

  /// Return pointer to vertex i's coordinates
  virtual double *GetVertex(int i) const = 0;

  virtual double GetVertex(int j, int i) const = 0;

  virtual ~Mesh_p();
};



/// Data type for mesh The vertex is given by template mechanism
template <int Dim>
class TMesh_p : virtual public Mesh_p {

protected:
  Array<Vertex<Dim> > vertices;

  void ConstructRectangularElements(Array<int> &n, Array<double> &L);

  void ConstructTriangularElements(Array<int> &n, Array<double> &L);

  void ConstructTetrahedralElements(Array<int> &n, Array<double> &L);

  void ConstructRectangularElements(Array<ifstream *> &files);

  void ConstructTriangularElements(Array<ifstream *> &files);

  void ConstructTetrahedralElements(Array<ifstream *> &files);

public:

  /// Return the spatial dimension.
  virtual int GetDim() const { return Dim; }

  /// Return pointer to vertex i's coordinates
  virtual double *GetVertex(int i) const { return vertices[i](); }

  /// Return pointer to cordinates coordinates
  virtual double GetVertex(int j, int i) const { return vertices[j](i); }

  /// Give the coordinates of every vertex of element i
  virtual void GetElementVerticesCoord(int i, Array<double *> &coord) const;

  /// Give the coordinates of every vertex of boundary element i
  virtual void GetBdrElementVerticesCoord(int i, Array<double *> &coord) const;


  /// Given n number of elements, create a 1D mesh on (xl, xr)
  TMesh_p<Dim>(int n, double xl, double xr);

  /// Create a 1D mesh given an array of x-coord.
  TMesh_p<Dim>(Array<double> &xcoord);

 /** Creates mesh for a rectangular domain (0,L[0])x(0,L[1]) or
     a brick domain (0,L[0])x(0,L[1])x(0,L[2]).
     For a rectangular domain, we divide it into n[0]*n[1] rectangles
     if type = QUADRILATERAL or into 2*n[0]*n[1] triangles if type = TRIANGLE.
     For a brick domain, we divide it into 2*n[0]*n[1]*n[2] TETRAHEDRALs
 **/
  TMesh_p<Dim>(Array<int> &n, Array<double> &L, Element::Type type);

  TMesh_p<Dim>(Array<ifstream *> &files, Element::Type type);

  virtual ~TMesh_p();
};


#endif
