#ifndef FILE_TRIGAUSSQUAD
#define FILE_TRIGAUSSQUAD

/*
  Data Type for Triangular Element Gaussian Quadrature
  
  Author: Q. Deng
  Date: 05/20/2015
*/

#include "array.hpp"
#include "../linalg/linalg_header.h"

using namespace std;

void GetStandTriQuadPW(int n, Array<double *> &pw);

void GetStandLineQuadPW(int n, Array<double *> &pw);

double GetElementInt(int n, double detJ, Array<double> &fv);

int ConvertN2PN(int n);

#endif
