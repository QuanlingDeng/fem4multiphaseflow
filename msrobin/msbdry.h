#ifndef FILE_MSBDRY
#define FILE_MSBDRY
 
/*
   Data type Subdomain
   A container of multiscale boundary data
   Author: Bradley McCaskill
*/

class msbdry
 {
 protected:
   
   int n;
   int locindex;
   int numseg;
   int gind;
   int neps;

   Array<double *> bdryfcts;

 public: 

   msbdry(int _numseg, int _n, int _locindex);
   msbdry(int _n, int _locindex);

   double getphifctval(int k, int p) { return bdryfcts[k][p]; }
   int getnumseg() { return numseg; }
   int getnumelem() { return n; }
   int getgbindex() { return gind; }
   int getlocindex() { return locindex; }
   int getneps(){ return neps; }

   void setgbindex(int gbindex){ gind = gbindex; }
   ~msbdry();
};

#endif
