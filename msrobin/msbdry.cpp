#include "msrobin_header.h"

//===========================================================================

msbdry::msbdry(int _numseg, int _n, int _locindex)
{
  n = _n;
  numseg = _numseg;
  locindex = _locindex;
  neps = n/numseg; 

  bdryfcts.SetSize(numseg+1);
  for(int k=0; k<=numseg; k++)
    {
      bdryfcts[k] = new double[n+1];
      for(int p=0; p<=n; p++)
	{
	  bdryfcts[k][p] = 0.0;
	} 
    }

  double m = 1.0/((double) neps);

  int ind = 0;
  for(int k=1; k<numseg; k++)
    {
      for(int p=0; p<=neps; p++)
	{
	  int ii = p + ind;
	  double val = p*m;
	  bdryfcts[k][ii] = val; 
	  bdryfcts[k][ii+neps] = 1.0-val;
	}
      ind+=neps;
    }
  
  ind = n-neps;
  for(int p=0; p<=neps; p++)
    {
      double val = p*m;
      bdryfcts[0][p] = 1.0-val;
      bdryfcts[numseg][p+ind] = val;
    }
}

//=============================================================================

msbdry::msbdry(int _n, int _locindex)
{
  n = _n;
  numseg = 0;
  locindex = _locindex;
  neps = n; 

  gind = -1;
}

//=============================================================================

msbdry::~msbdry()
{
  for(int k=0; k<bdryfcts.Size(); k++)
    {
      delete []bdryfcts[k];
    }
}
