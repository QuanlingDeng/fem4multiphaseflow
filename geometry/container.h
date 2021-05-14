#ifndef FILE_CONTAINER
#define FILE_CONTAINER
#include "geometry_header.h"

class Container  
{
 protected:
  Array<Vector *> data;
  
 public:


  ///Constructs an array of Vectors to facilitate the base Array in an Array of Arrays
  Container(ifstream *file, int m);

  Container(Array<Vector *> &input, int altsize);  

  inline int Size() {return data.Size();}
  
  inline Vector* ReturnData(int i) {return data[i]; };
  
  inline double GetVal(int i, int j) const{ return data[i]->Elem(j); };
  
  inline int GetSize(int i) const{ return data[i]->Size(); };

  void RD(Array<Vector *> &input);
  
  ~Container();
};
#endif
