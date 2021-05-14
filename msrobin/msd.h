#ifndef FILE_MSD
#define FILE_MSD
 
/*
   Data type Subdomain
   A container for the subdomains
   Author: Bradley McCaskill
*/

class msd
 {
 protected:

   Mesh *locmesh;
   Data *locdata;

   Array<double *> locellipticdata;
   Array<double *> locforcedata;

   Array<double *> dbdryval;
   Array<double *> nbdryval;
   Array<double *> rbdryval;
   Array<double *> rcoeff;

   Array<bool> BdrNeumann;
   Array<bool> BdrRobin;
   Array<bool> BdrDirichlet;
   Array<bool> GBdry;

   Array<double> locpointsourcevalues;
   Array<int> locpointsourceindex;

   Array<msbdry *> boundary;
   Array<Vector *> basisfcts;

   Array<int*> nbinfo;
   Array<int> gbvindex;

   int totvert;
   int dof;
   int tbf;
   int mtype;
   int nbdrcon;

 public: 

   msd(Array<int> &_n, Array<double> &_l, int _mtype, int _nbdrcon=1);

   Mesh *getmesh() { return locmesh; }
   Data *getdata() { return locdata; }

   int getnbdrs(){ return locmesh->GetNBdrs(); }
   int getmeshtype(){ return mtype; }
   int getnbdryelem(int b){ return locmesh->GetNBdryElem(b); }
   int getnumseg(int b){return boundary[b]->getnumseg(); }
   int getgbindex(int b){return boundary[b]->getgbindex(); }
   int getlocindex(int b) { return boundary[b]->getlocindex(); }
   int gettotvert(){return totvert; }
   int getne(){return locmesh->GetNE(); }
   int getnumbasis(){return tbf; }
   int getgbvindex(int i){return gbvindex[i]; }
   int getelementnvertices(int i){return locmesh->GetElementNVertices(i); }
   int getnbinfo(int j, int b){ return nbinfo[j][b];}
   int getbdrelemgindex(int b, int p);

   bool getGBdry(int b){return GBdry[b]; }

   double getbdryfctval(int b, int k, int p) { return boundary[b]->getphifctval(k, p); }
   double getbasisval(int k, int index){ return basisfcts[k]->Elem(index); }
   double getbdryelemlength(int b, int p);
   double getvertex(int j, int i){ return locmesh->GetVertex(j, i); }
   double getellipticdata(int k, int i=0){ return locellipticdata[i][k]; }
   double getrobincoeff(int b, int i=0) { return locdata->GetRobinCoeff(b, i); }

   void setgbindex(int b, int gbindex){ boundary[b]->setgbindex(gbindex); }
   void setgbvindex(int startindex);
   void setnbinfo(Array<int*> &_nbinfo);
   void setlocellipticdata(Function *permfunc);
   void setmaterialdata(Function *modfunc, double v);
   void setmaterialdata(Array<double *> &lamecoeff);
   void setlocforcedata(Function *forcefunc);
   void setlocforcedata(Function *forcefuncx, Function *forcefuncy);
   void setpointsourcedata(const Array<double> &_fvalue, const Array<int> &_nodeindex);
   void buildmsbasisfcts(Array<int> &_ns);
   void Build(Vector &sol, double* alpha);
   void setlocbndrydata(Array<double *> &_dbdryval,
			Array<double *> &_nbdryval,
			Array<double *> &_rbdryval,
			Array<double *> &_rcoeff,
			Array<bool> &_BdrDirichlet,
			Array<bool> &_BdrNeumann,
			Array<bool> &_BdrRobin, 
			Array<bool> &_GBdry);
   void getelementvertices(int i, Array<int> &ind){ locmesh->GetElementVertices(i, ind); }
   void printbasis(ofstream &out, int k){ locmesh->Print(out, *basisfcts[k]); }
   void printmesh(ofstream &out){ locmesh->Print(out, Element::QUADRILATERAL); }

   ~msd();
   
};

#endif
