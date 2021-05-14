#include "../linalg/linalg_header.h"
#include "../fem/fem_header.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include "../linalg/linalg_header.h"
#include "../mesh/mesh_header.h"

using namespace std;

int main()
{
  Array<ifstream *> files(1);
  files[0] = new ifstream( "quadrilateral_mesh.mh2" );

  Mesh *mesh;
  mesh = new TMesh<2>(files, Element::QUADRILATERAL);

  ofstream out;
  out.open("mesh_info.txt");
  mesh->Print(out, 5);
  out.close();

  delete mesh;
  for(int i=0; i<1; i++)
    {
      files[i]->close();
      delete files[i];
    }
}

