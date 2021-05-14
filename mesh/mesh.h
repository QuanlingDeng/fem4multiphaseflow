#ifndef FILE_MESH
#define FILE_MESH
#include "mesh_header.h"
/*
  Data type mesh
*/

class Mesh {

protected:
  Array<Element *> element;
  Array<Element *> bdrelement;
  Array<int*> vert_to_ele;
  Array<int*> BdrytoGlobalIndex;
  Array<int> be_to_edge;
  Array<int> NumOfBdryElem;
  Array<int> Degs;
  
  int Nx;
  int NumOfVertices;
  int NumOfElements;
  int NumOfBdrElements;
  int NumOfBdrs;
  int NumOfEdges;
  
public:

  /// Returns number of vertices.
  inline int GetNV() const { return NumOfVertices; }

  ///Returns number of elements in x direction from rectangular constructor. 
  inline int GetNx() const { return Nx; }

  //Returns number of elements in y direction from rectangular constructor. 
  inline int GetNy() const { return NumOfElements/Nx; }

  /// Returns number of elements.
  inline int GetNE() const { return NumOfElements; }

  /// Returns number of boundary elements.
  inline int GetNBE() const { return NumOfBdrElements; }

  /// Returns number of boundarys.
  inline int GetNBdrs() const { return NumOfBdrs; }

  /// Returns number of Edges.
  inline int GetNEdges() const { return NumOfEdges; }

  /// Returns number of boundary elements.
  inline int GetNBdryElem(int b) const { return NumOfBdryElem[b]; }

  /// Returns global index for boundary vertex
  inline int GetBdryGindex(int b, int i) const { return BdrytoGlobalIndex[b][i]; }

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

  int GetBdrElementEdgeIndex(int i) const {  return be_to_edge[i]; };

  /// Give the coordinates of every vertex of element i
  virtual void GetElementVerticesCoord(int i, Array<double *> &coord) const = 0;

  /// Give the coordinates of every vertex of boundary element i
  virtual void GetBdrElementVerticesCoord(int i, Array<double *> &coord) const = 0;

  /// Give the index of element's edges
  virtual void GetElementEdges(int el, Array<int> &edge_index) const = 0;

  virtual void GetEdgeElements(int ed, Array<int> &element_index) const = 0;

  virtual int GetElementToEdgeTable(Table &el_2_edge, Array<int> &be_to_edge) const = 0;

  virtual Table *GetEdgeToElementTable() const = 0;

  virtual void GetVertexToVertexTable(STable &v_to_v) const = 0;

  virtual void UpdateMesh(const Vector &displacement) = 0;

  /// Print the mesh information for glvis plot
  void Print(ofstream &out, int type) const;

  /// Print the mesh information for paraview plot
  void ParaviewPrint(int type) const;

  void ParaviewPrint(ofstream &out, const Vector &solution, double scalex=1.0, double scaley=1.0) const;

  void ParaviewPrintDeformedMesh(ofstream &out, const Vector &displacement, double scalex=1.0, double scaley=1.0) const; 

  /// Print the solution to a file for glvis plot
  void Print(ofstream &out, const Vector &solution) const;

 /// Print the solution information for paraview plot
  void ParaviewPrint(const Vector &solution, int type) const;

  /// Return the mesh type.
  virtual int GetMType() const = 0;

  /// Return the spatial dimension.
  virtual int GetDim() const = 0;

  /// Return pointer to vertex i's coordinates
  virtual double *GetVertex(int i) const = 0;

  virtual double GetVertex(int j, int i) const = 0;

  virtual ~Mesh();

  /// Return the index in vert_to_ele[i][j]
  int GetVertToEle(int i, int j);
};

/// Data type for mesh The vertex is given by template mechanism
template <int Dim>
class TMesh : virtual public Mesh {

protected:

  int mtype;

  Array<Vertex<Dim> > vertices;
  Table *el_to_edge;

  mutable Table *edge_to_el;
  
  void ConstructRectangularElements(Array<int> &n, Array<double> &L, bool generateedges=true);

  void ConstructTriangularElements(Array<int> &n, Array<double> &L, bool generateedges=true);

  void ConstructRectangularElements(Array<int> &n, Array<double*> &L, bool generateedges=true);

  void ConstructTriangularElements(Array<int> &n, Array<double*> &L, bool generateedges=true);

  void ConstructTetrahedralElements(Array<int> &n, Array<double> &L, bool generateedges=true);

  void ConstructRectangularElements(Array<ifstream *> &files, bool generateedges=true);

  void ConstructTriangularElements(Array<ifstream *> &files, bool generateedges=true);

  void ConstructTetrahedralElements(Array<ifstream *> &files, bool generateedges=true);

  //void Construct2DNURBSElements(Geometry &geom, bool generateedges=true);

  //void ConstructSurfaceNURBSElements(Geometry &geom, bool generateedges=true);

  //void Construct3DNURBSElements(Geometry &geom, bool generateedges=true);

public:

  /// Return the spatial dimension.
  virtual int GetMType() const { return mtype; }

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

  /// Give the index of element's edges
  virtual void GetElementEdges(int el, Array<int> &edge_index) const;

  virtual void GetEdgeElements(int ed, Array<int> &element_index) const;
  
  virtual int GetElementToEdgeTable(Table &el_2_edge, Array<int> &be_to_edge) const;

  virtual Table *GetEdgeToElementTable() const;
  
  virtual void GetVertexToVertexTable(STable &v_to_v) const;
  
  virtual void UpdateMesh(const Vector &displacement);

  /// Given n number of elements, create a 1D mesh on (xl, xr)
  TMesh<Dim>(int n, double xl, double xr);

  /// Create a 1D mesh given an array of x-coord.
  TMesh<Dim>(Array<double> &xcoord);

 /** Creates mesh for a rectangular domain (0,L[0])x(0,L[1]) or
     a brick domain (0,L[0])x(0,L[1])x(0,L[2]).
     For a rectangular domain, we divide it into n[0]*n[1] rectangles
     if type = QUADRILATERAL or into 2*n[0]*n[1] triangles if type = TRIANGLE.
     For a brick domain, we divide it into 2*n[0]*n[1]*n[2] TETRAHEDRALs
 **/
  TMesh<Dim>(Array<int> &n, Array<double> &L, Element::Type type, bool generateedges=true);

  TMesh<Dim>(Array<int> &n, Array<double*> &L, Element::Type type, bool generateedges=true);

  TMesh<Dim>(Array<ifstream *> &files, Element::Type type, bool generateedges=true);

  //TMesh<Dim>(Geometry &Geom, Element::Type type, bool generateedges=true);


  virtual ~TMesh();
};

#endif
