#ifndef FILE_FDATA
#define FILE_FDATA

/*
  Author: Q. Deng
  Data: 10/01/2015
*/

class FData
{
 protected:
  
  Array<double *> fdata;
  Array<double > fbdrydata;
  
 public:
  
  FData(Array<double*> &_fdata, Array<double > &_fbdrydata);

  double GetFVal(int k, int el){ return fdata[k][el]; }
  double GetBdryFVal(int i){ return fbdrydata[i]; }
  
  ~FData();
};

#endif
