#include "timedependent_header.h"

//=============================================================

FData::FData(Array<double*> &_fdata, Array<double > &_fbdrydata)
{
  fdata.SetSize(_fdata.Size());
  for(int i=0; i<fdata.Size(); i++){ fdata[i] = _fdata[i]; }

  fbdrydata.SetSize(_fbdrydata.Size());
  for(int i=0; i<fbdrydata.Size(); i++){ fbdrydata[i] = _fbdrydata[i]; }
}

//=============================================================

FData::~FData()
{
  // for(int i=0; i<fdata.Size(); i++){ delete []fdata[i]; }
}
