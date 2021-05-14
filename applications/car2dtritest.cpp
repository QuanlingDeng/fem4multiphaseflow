#include "../fem/fem_header.h"
#include "../mesh/mesh_header.h"

double k(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return 1.0;} //2+sin(x+2*y); }
double u(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return sin(x-2.0*y); }
double f(Array<double> &xy) { double x=xy[0]; double y=xy[1]; return 1.0;} //5.0*x*u(xy)-cos(x-2.0*y);}
double g(Array<double> &xy) { double x=xy[0]; double y=xy[1]; if( fabs(x-1.0) < 1.0e-32 ) return 2*(y-2); else return 0.0;}
double b(Array<double> &xy) { double x=xy[0]; double y=xy[1];
  if( fabs(x-1.0) < 1.0e-32 || fabs(x-12.0) < 1.0e-32 )
    return u(xy);
  else if( fabs(y-2.0) < 1.0e-32 ) 
    return u(xy);
  else if( fabs(y-3.0) < 1.0e-32 && x<=3.0 )
    return u(xy);
  else if( fabs(y-3.5) < 1.0e-32 && x>=9.5 )
    return u(xy);
  else if( fabs(x-y) < 1.0e-7 && x>=3.0 && x<=4.0)
    return u(xy);
  else if( fabs(2.0*x/3.0 + y - 59.0/6.0 ) < 1.0e-7 && x>=8.0 && x<=9.5)
    return u(xy);
  else if( fabs(x/8.0 - y + 3.5 ) < 1.0e-7 && x>=4.0 && x<=8.0)
    return u(xy);
  else 
    return 0.0;
}


int main(int argc, const char * argv[])
{
  char file1[256];
  char file2[256];
  char file3[256];

  sprintf(file1, "car.2.node");
  sprintf(file2, "car.2.ele");
  sprintf(file3, "car.2.poly");

  Array<ifstream *> files;
  files.SetSize(3);
  files[0] = new ifstream( file1 );
  files[1] = new ifstream( file2 );
  files[2] = new ifstream( file3 );

  Mesh *tmesh;
  tmesh = new TMesh<2>(files, Element::TRIANGLE);

  ofstream out;
  out.open("mesh.out");
  tmesh->Print(out, Element::TRIANGLE);
  out.close();

  int NV = tmesh->GetNV();
  //int NE = tmesh->GetNE();
  //int NBE = tmesh->GetNBE();

  Array<double *> edata(1);
  edata[0] = new double[NV];

  Array<double *> fdata(1);
  fdata[0] = new double[NV];
  
  Array<double *> ndata(1);
  ndata[0] = new double[NV];
  
  Vector uu(NV);
  
  Array<double> coord(2);
  for(int i=0; i<NV; i++)
    {
      coord[0] = tmesh->GetVertex(i, 0);
      coord[1] = tmesh->GetVertex(i, 1);
      edata[0][i] = k(coord);
      fdata[0][i] = f(coord);
      ndata[0][i] = g(coord);
      uu(i) = u(coord);
    }

  Data *tdata;
  tdata = new Data();
  tdata->SetEllipticData(edata);
  tdata->SetForceData(fdata);
  tdata->SetNeumannData(ndata);

  FEM *UW = new LinearFUPG2D(tmesh, tdata);
  UW->Assemble();

  cout<<" 1234 "<<endl;

  int NB = 0;
  NB = tmesh->GetNBdrs();
  Array<bool> bdry_is_dirichlet(NB);
  for(int i=0; i<NB; i++)
    bdry_is_dirichlet[i] = true;  
  bdry_is_dirichlet[NB-1] = false;

  Vector Dval(NV);
  Dval = 0.0;
  //UW->ProjectAFunction(b, Dval);

  UW->SetDirichletBoundary(bdry_is_dirichlet, Dval);

  Vector sol(NV);
  sol = Dval;

  char solv[20] = "pcg";
  UW->Solve(sol, solv, 100, 1, 1e-10, 12-20);

  tmesh->ParaviewPrint(Element::TRIANGLE);
  tmesh->ParaviewPrint(sol, Element::TRIANGLE);  
  
  delete []ndata[0];
  delete []edata[0];
  delete []fdata[0];
  delete tdata;
  delete tmesh;

  return 0;
}
