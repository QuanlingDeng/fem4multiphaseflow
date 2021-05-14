#include "../fem/fem_header.h"
#include "../linalg/linalg_header.h"
#include "../general/array.hpp"
#include "../geometry/geometry_header.h"
#include "../timedependent/time_header.h"
#include <fstream>
#include <cmath>
int main()
{
  char file1[256];
  sprintf(file1, "geometry.info");
  ifstream *files = new ifstream(file1);
  Geometry *uw = new Geometry(files);
  files->close();
  delete files;
   Array<int *> h(2);
  h[0] = new int[2];
  h[0][0] = 5;
  h[0][1] = 2;
  h[1] = new int[2];
  h[1][0] = 1;
  h[1][1] = 1;
  int *n = new int[2];
  n[0]=60;
  n[1]=60;
  Mesh *tmesh;
  tmesh = new TMesh<2> (*uw, Element::TwoDNURBS);
  // uw->Refine(h);
  // uw->PrintNURBS(n);
  //cout << "new knot is ";
  /*  for (int i=0; i<12; i++)
    {
      cout << uw->GetKnot(0,0,i) << "  ";
      }
  cout  << endl << "new ctrl pts are " << endl;
  for (int j=0; j<4; j++)
    {
      for (int i=0; i<3; i++)
	{
	  cout << uw->GetCtrlPt(0,i,j*2) << "  " ;
	}
      cout << uw->GetWeight(0,j*2) << endl;
      }
  cout << "break" << endl;
  for (int j=0; j<3; j++)
    {
      for (int i=0; i<3; i++)
	{
	  cout << uw->GetCtrlPt(1,i,j*2+1) << "  " ;
	}
      cout << uw->GetWeight(1,j*2+1) << endl;
      }


  //cout << uw->GetKnot(0,0,0) << "  "<< uw->GetKnot(0,0,1) << "  "<< uw->GetKnot(0,0,2) << "  "<< uw->GetKnot(0,0,3) << "  "<< uw->GetKnot(0,0,4) << "  "<< uw->GetKnot(0,0,5) << "  "<< uw->GetKnot(0,0,6) << "  " << uw->GetKnot(0,0,7) << "  "<< uw->GetKnot(0,0,8) << "  "<< uw->GetKnot(0,0,9) << "  "<< uw->GetKnot(0,0,10) << "  "<< uw->GetKnot(0,0,11) << endl;
  // cout << "===============================" << endl << uw->GetKnot(0,1,0) << "  "<< uw->GetKnot(0,1,1) << "  " << uw->GetKnot(0,1,2) << endl;*/
  delete uw;
  delete []h[0];
  delete []h[1];
  delete []n;
  //delete tmesh;
}
