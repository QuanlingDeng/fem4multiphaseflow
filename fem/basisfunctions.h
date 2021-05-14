#ifndef FILE_BASISFUNCTIONS
#define FILE_BASISFUNCTIONS

/*
   Data type for basisfunctions

   Author: Quanling Deng
   Date: 09/04/2015
*/

/* This is a container for all the basis functions of finite element methods */
class BasisFunctions
{
protected:

  int order;
  
public:

  BasisFunctions(int _order);

  double BF1D(int ind, double x);
  double DerBF1D(int ind, double x);

  double BF2D(int ind, double x, double y);
  double GradBF2D(int ind, int dimflag, double x, double y);

  double BF3D(int ind, double x, double y, double z) { return 0.0; }
  double GradBF3D(int ind, int dimflag, double x, double y, double z) { return 0.0; }

  ~BasisFunctions();
};

#endif
