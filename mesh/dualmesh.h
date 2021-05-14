#ifndef FILE_DUALMESH
#define FILE_DUALMESH

/*
  Data type for DualMesh

  Author: Q. Deng
  Date: 09/05/2015
*/


class DualMesh
{
protected:

  Mesh *mesh;
  int dof;
  int dualmeshorder;
  int dualmeshtype;
  
public:

  /* Constructor for dualmesh. _dualmeshorder specifies dualmesh 
     with order 1 for linear, 2 for quadratic, and 3 for cubic. 
     _dualmeshtype specifies a type for a certain order of dualmesh. */
  DualMesh(Mesh *_mesh, int _dualmeshorder=1, int _dualmeshtype=0);

  // return the original mesh
  inline Mesh *GetMesh() { return mesh; }

  inline int GetDualMeshNumDOF() { return dof; }

  inline int GetDualMeshOrder() { return dualmeshorder; }

  inline int GetDualMeshType() { return dualmeshtype; }

  // return the dual info in an element
  void GetElemDualInfo(int i, Array<double *> &coord, Array<double *> &normals,
		       Array<double> &elengths, Array<double> &areas);
  
  void GetElemDualCoord(int i, Array<double *> &coord);
    
  void GetElemDualNormals(int i, Array<double *> &normals);

  void GetElemDualEdgeLengths(int i, Array<double> &elengths);
  
  void GetElemDualAreas(int i, Array<double> &areas);

  //void PrintDualMesh();

  //void PrintSaturation(ofstream &out, Vector &Saturation, double scalex = 1.0, double scaley = 1.0);

  //void PrintGLESaturation(ofstream &out, Vector &Saturation, double scalex = 0.1, double scaley = 10.0);

  ~DualMesh();

};

#endif
