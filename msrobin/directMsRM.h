#ifndef FILE_DIRECTMSRM
#define FILE_DIRECTMSRM
 
/*
   Data type MSRM
   An implementation of the non iterative multiscale Robin method.
   Author: Bradley McCaskill
*/

class directMsRM 
 {
 protected:

   Array<Vector *> gsol;
   Array<double *> alpha;
   Array<msd *> subdomains;

   SparseMatrix *mat;
   Vector rhs;
   Vector solalpha;

   int Nsd;
   int dim;
   int nbdrcon;

 public: 
   
   directMsRM(Array<msd*> &_subdomains, int &_dim, int _numbdrcon=1);
   double intbdryfcts(Array<double *> &bdryfcts, int n, Array<double> &h);
   double getgsol(int k, int i){ return gsol[k]->Elem(i); }
   double computel2error(Function *exactsol);
   double computel2error(Function *exactsolx, Function *exactsoly);
   double computeh1error(Array<Function *> &exactderivative);

   void solve(); 
   void addlocalsystem(int b, int bb, int k);
   void buildsol(Vector &sol, int Nx, int Ny, int nx, int ny);
   void paraviewprintmesh(); 
   void paraviewprintdeformedmesh(); 
   void paraviewprintsol(); 
   void paraviewprintellipticdata(); 
   void solve(Vector &sol, char *linalgsolv, int maxit=1000, int printit=0, 
	     double rtol=1e-14, double atol=1e-14);

   ~directMsRM();
   
};

#endif
