#include "../fem/fem_header.h"

using namespace std;

int main(int argc, const char * argv[])
{
  /*
  int Nx = 4;
  int Ny = 2;
  int dof = ( Nx+1 ) * ( Ny+1 );
  double hx = 1.0/Nx;
  double hy = 1.0/Ny;
  int Nb = 2*Nx + 2*Ny;

  Array<int> n(2);
  n[0] = Nx;
  n[1] = Ny;

  Array<double> L(4);
  L[0] = 0.0;
  L[1] = 1.0;
  L[2] = 0.0;
  L[3] = 1.0;

  Mesh *tmesh;
  tmesh = new TMesh<2> (n, L, Element::TRIANGLE);
  */
 
  char file1[256];
  char file2[256];
  char file3[256];

  sprintf(file1, "E.2.node");
  sprintf(file2, "E.2.ele");
  sprintf(file3, "E.2.poly");

  //sprintf(file1, "car.2.node");
  //sprintf(file2, "car.2.ele");
  //sprintf(file3, "car.2.poly");

  //sprintf(file1, "cookmem.2.node");
  //sprintf(file2, "cookmem.2.ele");
  //sprintf(file3, "cookmem.2.poly");

  Array<ifstream *> files;
  files.SetSize(3);
  files[0] = new ifstream( file1 );
  files[1] = new ifstream( file2 );
  files[2] = new ifstream( file3 );

  Mesh *mesh;
  mesh = new TMesh<2>(files, Element::TRIANGLE);

  ofstream out;
  out.open("trimesh.out");
  mesh->Print(out, Element::TRIANGLE);
  out.close();


  delete mesh;
  for(int i=0; i<3; i++)
    {
      files[i]->close();
      delete files[i];
    }
  
  //delete tmesh;

  return 0;
}

