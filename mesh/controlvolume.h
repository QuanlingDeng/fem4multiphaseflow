#ifndef FILE_CONTROLVOLUME
#define FILE_CONTROLVOLUME


/*
  Author: Q. Deng
  Data: 10/14/2015
*/

class ControlVolume
{
 protected:

  int numdof;
  Array<int> vindex;
  Array<double> elengths;
  Array<double *> coord;
  Array<double *> vcoord;
  Array<double *> normals;
  Array<double> areas;

 public:

  ControlVolume(const Array<double *> &v, const Array<int> &dof);

  void GetEdgeLengths(Array<double> &_elenghts);

  void GetCoordinates(Array<double *> &_coord);

  void GetElementVerticesCoord(Array<double *> &_vcoord);
  
  void GetNormals(Array<double *> &_normals);

  void GetAreas(Array<double> &_areas);

  void GetDOF(Array<int> &_dof);

  double GetArea(int k){ return areas[k]; }

  ~ControlVolume();
};

#endif
