#include "mpi_fem_header.h"

//==============================================================================

LinearFEM3D_p::LinearFEM3D_p(Mesh_p *_mesh, Data *_data) :
FEM_p(_mesh, _data)
{
    strcpy(Type, "LinearFEM3D_p");
    NumGlobalDOF = mesh->GetNV();
    NumLocalDOF = 4;
    rhs.SetSize(NumGlobalDOF);
    mat = new SparseMatrix(NumGlobalDOF, 45);
}

//============================================================================== 

void LinearFEM3D_p::GetElementDOF(int i, Array<int> &ind)
{
  mesh->GetElementVertices(i, ind);
}
//==============================================================================

void LinearFEM3D_p::GetBdrElementDOF(int i, Array<int> &ind)
{
  mesh->GetBdrElementVertices(i, ind);
}

//==============================================================================

void LinearFEM3D_p::ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalEllipticCoeff(ind[j]);

  double k1 = locdat[0];
  double k2 = locdat[1];
  double k3 = locdat[2];
  double k4 = locdat[3];

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[mesh->GetNV()];

  mesh->GetElementVerticesCoord(i, coord);
  double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
  x1 = coord[0][0];
  x2 = coord[0][1];
  x3 = coord[0][2];
  x4 = coord[0][3];
  y1 = coord[1][0];
  y2 = coord[1][1];
  y3 = coord[1][2];
  y4 = coord[1][3];
  z1 = coord[2][0];
  z2 = coord[2][1];
  z3 = coord[2][2];
  z4 = coord[2][3];
  
  //=====================================================================================
  //row1
  //=====================================================================================
locmat[0][0] = ((k1 + k2 + k3 + k4)*(pow(x2,2)*pow(y3,2) - 2*pow(x2,2)*y3*y4 + pow(x2,2)*pow(y4,2) + pow(y3,2)*pow(z2,2) - 2*y3*y4*pow(z2,2) + pow(y4,2)*pow(z2,2) + pow(x4,2)*(pow(y2,2) - 2*y2*y3 + pow(y3,2) + pow(z2 - z3,2)) - 2*y2*y3*z2*z3 + 2*y2*y4*z2*z3 + 2*y3*y4*z2*z3 - 2*pow(y4,2)*z2*z3 + pow(x2,2)*pow(z3,2) + pow(y2,2)*pow(z3,2) - 2*y2*y4*pow(z3,2) + pow(y4,2)*pow(z3,2) + pow(x3,2)*(pow(y2,2) - 2*y2*y4 + pow(y4,2) + pow(z2 - z4,2)) - 2*x2*x4*(pow(y3,2) - y3*y4 + y2*(-y3 + y4) - (z2 - z3)*(z3 - z4)) + 2*y2*y3*z2*z4 - 2*pow(y3,2)*z2*z4 - 2*y2*y4*z2*z4 + 2*y3*y4*z2*z4 - 2*pow(x2,2)*z3*z4 - 2*pow(y2,2)*z3*z4 + 2*y2*y3*z3*z4 + 2*y2*y4*z3*z4 - 2*y3*y4*z3*z4 + pow(x2,2)*pow(z4,2) + pow(y2,2)*pow(z4,2) - 2*y2*y3*pow(z4,2) + pow(y3,2)*pow(z4,2) - 2*x3*(x4*(pow(y2,2) + y3*y4 - y2*(y3 + y4) + (z2 - z3)*(z2 - z4)) + x2*(y2*(y3 - y4) - y3*y4 + pow(y4,2) + z2*z3 - z2*z4 - z3*z4 + pow(z4,2)))))/   (24.*(-(x2*y3*z1) + x2*y4*z1 + x1*y3*z2 - x1*y4*z2 + x2*y1*z3 - x1*y2*z3 + x1*y4*z3 - x2*y4*z3 + x4*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) - x2*y1*z4 + x1*y2*z4 - x1*y3*z4 + x2*y3*z4 + x3*(-(y4*z1) - y1*z2 + y4*z2 + y2*(z1 - z4) + y1*z4))); 

locmat[0][1] = -((k1 + k2 + k3 + k4)*(x1*x2*pow(y3,2) - 2*x1*x2*y3*y4 + x1*x2*pow(y4,2) + pow(y3,2)*z1*z2 - 2*y3*y4*z1*z2 + pow(y4,2)*z1*z2 - y2*y3*z1*z3 + y2*y4*z1*z3 + y3*y4*z1*z3 - pow(y4,2)*z1*z3 - y1*y3*z2*z3 + y1*y4*z2*z3 + y3*y4*z2*z3 - pow(y4,2)*z2*z3 + x1*x2*pow(z3,2) + y1*y2*pow(z3,2) - y1*y4*pow(z3,2) - y2*y4*pow(z3,2) + pow(y4,2)*pow(z3,2) + pow(x4,2)*(y1*(y2 - y3) - y2*y3 + pow(y3,2) + z1*z2 - z1*z3 - z2*z3 + pow(z3,2)) + x4*(x2*(-pow(y3,2) + y1*(y3 - y4) + y3*y4 + (z1 - z3)*(z3 - z4)) + x1*(-pow(y3,2) + y2*(y3 - y4) + y3*y4 + (z2 - z3)*(z3 - z4))) + y2*y3*z1*z4 - pow(y3,2)*z1*z4 - y2*y4*z1*z4 + y3*y4*z1*z4 + y1*y3*z2*z4 - pow(y3,2)*z2*z4 - y1*y4*z2*z4 + y3*y4*z2*z4 - 2*x1*x2*z3*z4 - 2*y1*y2*z3*z4 + y1*y3*z3*z4 + y2*y3*z3*z4 + y1*y4*z3*z4 + y2*y4*z3*z4 - 2*y3*y4*z3*z4 + x1*x2*pow(z4,2) + y1*y2*pow(z4,2) - y1*y3*pow(z4,2) - y2*y3*pow(z4,2) + pow(y3,2)*pow(z4,2) + pow(x3,2)*(y1*(y2 - y4) - y2*y4 + pow(y4,2) + z1*z2 - z1*z4 - z2*z4 + pow(z4,2)) + x3*(x4*(-2*y3*y4 + y2*(y3 + y4) + y1*(-2*y2 + y3 + y4) - 2*z1*z2 + z1*z3 + z2*z3 + z1*z4 + z2*z4 - 2*z3*z4) - x2*(y1*(y3 - y4) - y3*y4 + pow(y4,2) + z1*z3 - z1*z4 - z3*z4 + pow(z4,2)) - x1*(y2*(y3 - y4) - y3*y4 + pow(y4,2) + z2*z3 - z2*z4 - z3*z4 + pow(z4,2)))))/ (24.*(-(x2*y3*z1) + x2*y4*z1 + x1*y3*z2 - x1*y4*z2 + x2*y1*z3 - x1*y2*z3 + x1*y4*z3 - x2*y4*z3 + x4*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) - x2*y1*z4 + x1*y2*z4 - x1*y3*z4 + x2*y3*z4 + x3*(-(y4*z1) - y1*z2 + y4*z2 + y2*(z1 - z4) + y1*z4)));

locmat[0][2] = -((k1 + k2 + k3 + k4)*(-(pow(x4,2)*y1*y2) - x1*x4*pow(y2,2) + pow(x4,2)*pow(y2,2) + pow(x4,2)*y1*y3 + x1*x4*y2*y3 - pow(x4,2)*y2*y3 + x1*x4*y2*y4 - x1*x4*y3*y4 - pow(x4,2)*z1*z2 - y2*y3*z1*z2 + y2*y4*z1*z2 + y3*y4*z1*z2 - pow(y4,2)*z1*z2 - x1*x4*pow(z2,2) + pow(x4,2)*pow(z2,2) + y1*y3*pow(z2,2) - y1*y4*pow(z2,2) - y3*y4*pow(z2,2) + pow(y4,2)*pow(z2,2) + pow(x4,2)*z1*z3 + pow(y2,2)*z1*z3 - 2*y2*y4*z1*z3 + pow(y4,2)*z1*z3 + x1*x4*z2*z3 - pow(x4,2)*z2*z3 - y1*y2*z2*z3 + y1*y4*z2*z3 + y2*y4*z2*z3 - pow(y4,2)*z2*z3 + x3*(x4*(-pow(y2,2) + y1*(y2 - y4) + y2*y4 + (z1 - z2)*(z2 - z4)) + x1*(pow(y2,2) - 2*y2*y4 + pow(y4,2) + pow(z2 - z4,2))) - pow(y2,2)*z1*z4 + y2*y3*z1*z4 + y2*y4*z1*z4 - y3*y4*z1*z4 +  x1*x4*z2*z4 + y1*y2*z2*z4 - 2*y1*y3*z2*z4 + y2*y3*z2*z4 + y1*y4*z2*z4 - 2*y2*y4*z2*z4 + y3*y4*z2*z4 - x1*x4*z3*z4 + y1*y2*z3*z4 - pow(y2,2)*z3*z4 - y1*y4*z3*z4 + y2*y4*z3*z4 - y1*y2*pow(z4,2) + pow(y2,2)*pow(z4,2) + y1*y3*pow(z4,2) - y2*y3*pow(z4,2) + pow(x2,2)*(y1*(y3 - y4) - y3*y4 + pow(y4,2) + z1*z3 - z1*z4 - z3*z4 + pow(z4,2)) + x2*(x4*(y2*y3 - 2*y2*y4 + y3*y4 + y1*(y2 - 2*y3 + y4) + z1*z2 - 2*z1*z3 + z2*z3 + z1*z4 - 2*z2*z4 + z3*z4) - x3*(y1*(y2 - y4) - y2*y4 + pow(y4,2) + z1*z2 - z1*z4 - z2*z4 + pow(z4,2)) - x1*(y2*(y3 - y4) - y3*y4 + pow(y4,2) + z2*z3 - z2*z4 - z3*z4 + pow(z4,2)))))/   (24.*(-(x2*y3*z1) + x2*y4*z1 + x1*y3*z2 - x1*y4*z2 + x2*y1*z3 - x1*y2*z3 + x1*y4*z3 - x2*y4*z3 + x4*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) - x2*y1*z4 + x1*y2*z4 - x1*y3*z4 + x2*y3*z4 + 	 x3*(-(y4*z1) - y1*z2 + y4*z2 + y2*(z1 - z4) + y1*z4)));

 locmat[0][3] = -((k1 + k2 + k3 + k4)*(x1*x4*pow(y2,2) - 2*x1*x4*y2*y3 + x1*x4*pow(y3,2) + y2*y3*z1*z2 - pow(y3,2)*z1*z2 - y2*y4*z1*z2 +  y3*y4*z1*z2 + x1*x4*pow(z2,2) - y1*y3*pow(z2,2) + pow(y3,2)*pow(z2,2) + y1*y4*pow(z2,2) - y3*y4*pow(z2,2) - pow(y2,2)*z1*z3 + y2*y3*z1*z3 + y2*y4*z1*z3 - y3*y4*z1*z3 - 2*x1*x4*z2*z3 + y1*y2*z2*z3 + y1*y3*z2*z3 - 2*y2*y3*z2*z3 -  2*y1*y4*z2*z3 + y2*y4*z2*z3 + y3*y4*z2*z3 + x1*x4*pow(z3,2) - y1*y2*pow(z3,2) + pow(y2,2)*pow(z3,2) +  y1*y4*pow(z3,2) - y2*y4*pow(z3,2) + x3*(-(x4*(pow(y2,2) - y2*y3 + y1*(-y2 + y3) - (z1 - z2)*(z2 - z3))) -  x1*(pow(y2,2) + y3*y4 - y2*(y3 + y4) + (z2 - z3)*(z2 - z4))) + pow(x3,2)*(pow(y2,2) - y2*y4 + y1*(-y2 + y4) - (z1 - z2)*(z2 - z4)) +  pow(x2,2)*(pow(y3,2) - y3*y4 + y1*(-y3 + y4) - (z1 - z3)*(z3 - z4)) + pow(y2,2)*z1*z4 - 2*y2*y3*z1*z4 + pow(y3,2)*z1*z4 - y1*y2*z2*z4 + y1*y3*z2*z4 + y2*y3*z2*z4 - pow(y3,2)*z2*z4 + y1*y2*z3*z4 - pow(y2,2)*z3*z4 - y1*y3*z3*z4 + y2*y3*z3*z4 + x2*(-(x4*(y1*(y2 - y3) - y2*y3 + pow(y3,2) + z1*z2 - z1*z3 - z2*z3 + pow(z3,2))) + x1*(-pow(y3,2) + y2*(y3 - y4) + y3*y4 + (z2 - z3)*(z3 - z4)) + x3*(y1*(y2 + y3 - 2*y4) + y3*y4 + y2*(-2*y3 + y4) + z1*z2 + z1*z3 - 2*z2*z3 - 2*z1*z4 + z2*z4 + z3*z4))))/ (24.*(-(x2*y3*z1) + x2*y4*z1 + x1*y3*z2 - x1*y4*z2 + x2*y1*z3 - x1*y2*z3 + x1*y4*z3 - x2*y4*z3 + x4*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) - x2*y1*z4 + x1*y2*z4 - x1*y3*z4 + x2*y3*z4 + x3*(-(y4*z1) - y1*z2 + y4*z2 + y2*(z1 - z4) + y1*z4)));


  //=====================================================================================
  //row2
  //=====================================================================================
 locmat[1][0] = locmat[0][1];

locmat[1][1] = ((k1 + k2 + k3 + k4)*(pow(x1,2)*pow(y3,2) - 2*pow(x1,2)*y3*y4 + pow(x1,2)*pow(y4,2) + pow(y3,2)*pow(z1,2) - 2*y3*y4*pow(z1,2) + pow(y4,2)*pow(z1,2) + pow(x4,2)*(pow(y1 - y3,2) + pow(z1 - z3,2)) - 2*y1*y3*z1*z3 + 2*y1*y4*z1*z3 + 2*y3*y4*z1*z3 - 2*pow(y4,2)*z1*z3 + pow(x1,2)*pow(z3,2) + pow(y1,2)*pow(z3,2) - 2*y1*y4*pow(z3,2) + pow(y4,2)*pow(z3,2) - 2*x3*(x4*((y1 - y3)*(y1 - y4) + (z1 - z3)*(z1 - z4)) + x1*((y1 - y4)*(y3 - y4) + (z1 - z4)*(z3 - z4))) + pow(x3,2)*(pow(y1 - y4,2) + pow(z1 - z4,2)) + 2*x1*x4*((y1 - y3)*(y3 - y4) + (z1 - z3)*(z3 - z4)) - 2*(pow(x1,2)*z3 + (y1 - y3)*(-(y3*z1) + y4*z1 + y1*z3 - y4*z3))*z4 + (pow(x1,2) + pow(y1 - y3,2))*pow(z4,2)))/(24.*(-(x2*y3*z1) + x2*y4*z1 + x1*y3*z2 - x1*y4*z2 + x2*y1*z3 - x1*y2*z3 + x1*y4*z3 - x2*y4*z3 + x4*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) + (-(x2*y1) + x1*y2 - x1*y3 + x2*y3)*z4 + x3*(-(y4*z1) - y1*z2 + y4*z2 + y2*(z1 - z4) + y1*z4)));

locmat[1][2] = -((k1 + k2 + k3 + k4)*(pow(x4,2)*pow(y1,2) + x1*x4*y1*y2 - pow(x4,2)*y1*y2 + x1*x4*y1*y3 - pow(x4,2)*y1*y3 + pow(x1,2)*y2*y3 - 2*x1*x4*y2*y3 + pow(x4,2)*y2*y3 - 2*x1*x4*y1*y4 - pow(x1,2)*y2*y4 + x1*x4*y2*y4 - pow(x1,2)*y3*y4 + x1*x4*y3*y4 + pow(x1,2)*pow(y4,2) + pow(x4,2)*pow(z1,2) + y2*y3*pow(z1,2) - y2*y4*pow(z1,2) - y3*y4*pow(z1,2) + pow(y4,2)*pow(z1,2) + x1*x4*z1*z2 - pow(x4,2)*z1*z2 - y1*y3*z1*z2 + y1*y4*z1*z2 + y3*y4*z1*z2 - pow(y4,2)*z1*z2 + x1*x4*z1*z3 - pow(x4,2)*z1*z3 - y1*y2*z1*z3 + y1*y4*z1*z3 + y2*y4*z1*z3 - pow(y4,2)*z1*z3 + pow(x1,2)*z2*z3 - 2*x1*x4*z2*z3 + pow(x4,2)*z2*z3 + pow(y1,2)*z2*z3 - 2*y1*y4*z2*z3 + pow(y4,2)*z2*z3 - 2*x1*x4*z1*z4 + y1*y2*z1*z4 + y1*y3*z1*z4 - 2*y2*y3*z1*z4 - 2*y1*y4*z1*z4 + y2*y4*z1*z4 + y3*y4*z1*z4 - pow(x1,2)*z2*z4 + x1*x4*z2*z4 - pow(y1,2)*z2*z4 + y1*y3*z2*z4 + y1*y4*z2*z4 - y3*y4*z2*z4 - pow(x1,2)*z3*z4 + x1*x4*z3*z4 - pow(y1,2)*z3*z4 + y1*y2*z3*z4 + y1*y4*z3*z4 - y2*y4*z3*z4 + pow(x1,2)*pow(z4,2) + pow(y1,2)*pow(z4,2) - y1*y2*pow(z4,2) - y1*y3*pow(z4,2) + y2*y3*pow(z4,2) + x3*(-(x4*(pow(y1,2) + y2*y4 - y1*(y2 + y4) + (z1 - z2)*(z1 - z4))) - x1*(y1*(y2 - y4) - y2*y4 + pow(y4,2) + z1*z2 - z1*z4 - z2*z4 + pow(z4,2))) + x2*(-(x4*(pow(y1,2) + y3*y4 - y1*(y3 + y4) + (z1 - z3)*(z1 - z4))) + x3*(pow(y1,2) - 2*y1*y4 + pow(y4,2) + pow(z1 - z4,2)) - x1*(y1*(y3 - y4) - y3*y4 + pow(y4,2) + z1*z3 - z1*z4 - z3*z4 + pow(z4,2)))))/(24.*(-(x2*y3*z1) + x2*y4*z1 + x1*y3*z2 - x1*y4*z2 + x2*y1*z3 - x1*y2*z3 + x1*y4*z3 - x2*y4*z3 +  x4*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) - x2*y1*z4 + x1*y2*z4 - x1*y3*z4 + x2*y3*z4 + x3*(-(y4*z1) - y1*z2 + y4*z2 + y2*(z1 - z4) + y1*z4)));

 locmat[1][3] = ((-k1 - k2 - k3 - k4)*(-(x1*x4*y1*y2) + x1*x4*y1*y3 - pow(x1,2)*y2*y3 + x1*x4*y2*y3 + pow(x1,2)*pow(y3,2) - x1*x4*pow(y3,2) + pow(x1,2)*y2*y4 - pow(x1,2)*y3*y4 - y2*y3*pow(z1,2) + pow(y3,2)*pow(z1,2) + y2*y4*pow(z1,2) - y3*y4*pow(z1,2) - x1*x4*z1*z2 + y1*y3*z1*z2 - pow(y3,2)*z1*z2 - y1*y4*z1*z2 + y3*y4*z1*z2 + x1*x4*z1*z3 + y1*y2*z1*z3 - 2*y1*y3*z1*z3 + y2*y3*z1*z3 + y1*y4*z1*z3 - 2*y2*y4*z1*z3 + y3*y4*z1*z3 - pow(x1,2)*z2*z3 + x1*x4*z2*z3 - pow(y1,2)*z2*z3 + y1*y3*z2*z3 + y1*y4*z2*z3 - y3*y4*z2*z3 + pow(x1,2)*pow(z3,2) - x1*x4*pow(z3,2) + pow(y1,2)*pow(z3,2) - y1*y2*pow(z3,2) - y1*y4*pow(z3,2) + y2*y4*pow(z3,2) + x2*(x4*(pow(y1,2) - 2*y1*y3 + pow(y3,2) + pow(z1 - z3,2)) - x3*(pow(y1,2) + y3*y4 - y1*(y3 + y4) + (z1 - z3)*(z1 - z4)) + x1*(-pow(y3,2) + y1*(y3 - y4) + y3*y4 + (z1 - z3)*(z3 - z4))) + pow(x3,2)*(pow(y1,2) + y2*y4 - y1*(y2 + y4) + (z1 - z2)*(z1 - z4)) - y1*y2*z1*z4 + y1*y3*z1*z4 + y2*y3*z1*z4 -  pow(y3,2)*z1*z4 + pow(x1,2)*z2*z4 + pow(y1,2)*z2*z4 - 2*y1*y3*z2*z4 + pow(y3,2)*z2*z4 - pow(x1,2)*z3*z4 - pow(y1,2)*z3*z4 + y1*y2*z3*z4 + y1*y3*z3*z4 - y2*y3*z3*z4 + x3*(-(x4*(pow(y1,2) + y2*y3 - y1*(y2 + y3) + (z1 - z2)*(z1 - z3))) + x1*(y2*(y3 - 2*y4) + y3*y4 + y1*(y2 - 2*y3 + y4) + z1*z2 - 2*z1*z3 + z2*z3 + z1*z4 - 2*z2*z4 + z3*z4))))/ (24.*(-(x2*y3*z1) + x2*y4*z1 + x1*y3*z2 - x1*y4*z2 + x2*y1*z3 - x1*y2*z3 + x1*y4*z3 - x2*y4*z3 +  x4*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) - x2*y1*z4 + x1*y2*z4 - x1*y3*z4 + x2*y3*z4 + x3*(-(y4*z1) - y1*z2 + y4*z2 + y2*(z1 - z4) + y1*z4))); 

  //=====================================================================================
  //row3
  //=====================================================================================

 locmat[2][0] = locmat[0][2];

 locmat[2][1] = locmat[1][2];

locmat[2][2] = ((k1 + k2 + k3 + k4)*(pow(x1,2)*pow(y2,2) - 2*pow(x1,2)*y2*y4 + pow(x1,2)*pow(y4,2) + pow(y2,2)*pow(z1,2) - 2*y2*y4*pow(z1,2) + pow(y4,2)*pow(z1,2) + pow(x4,2)*(pow(y1 - y2,2) + pow(z1 - z2,2)) - 2*y1*y2*z1*z2 + 2*y1*y4*z1*z2 + 2*y2*y4*z1*z2 - 2*pow(y4,2)*z1*z2 + pow(x1,2)*pow(z2,2) + pow(y1,2)*pow(z2,2) - 2*y1*y4*pow(z2,2) + pow(y4,2)*pow(z2,2) - 2*x2*(x4*((y1 - y2)*(y1 - y4) + (z1 - z2)*(z1 - z4)) + x1*((y1 - y4)*(y2 - y4) + (z1 - z4)*(z2 - z4))) + pow(x2,2)*(pow(y1 - y4,2) + pow(z1 - z4,2)) + 2*x1*x4*((y1 - y2)*(y2 - y4) + (z1 - z2)*(z2 - z4)) - 2*(pow(x1,2)*z2 + (y1 - y2)*(-(y2*z1) + y4*z1 + y1*z2 - y4*z2))*z4 + (pow(x1,2) + pow(y1 - y2,2))*pow(z4,2)))/(24.*(-(x2*y3*z1) + x2*y4*z1 + x1*y3*z2 - x1*y4*z2 + x2*y1*z3 - x1*y2*z3 + x1*y4*z3 - x2*y4*z3 + x4*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) + (-(x2*y1) + x1*y2 - x1*y3 + x2*y3)*z4 + x3*(-(y4*z1) - y1*z2 + y4*z2 + y2*(z1 - z4) + y1*z4)));

locmat[2][3] = -((k1 + k2 + k3 + k4)*(x1*x4*y1*y2 + pow(x1,2)*pow(y2,2) - x1*x4*pow(y2,2) - x1*x4*y1*y3 - pow(x1,2)*y2*y3 + x1*x4*y2*y3 - pow(x1,2)*y2*y4 + pow(x1,2)*y3*y4 + pow(y2,2)*pow(z1,2) - y2*y3*pow(z1,2) - y2*y4*pow(z1,2) + y3*y4*pow(z1,2) + x1*x4*z1*z2 - 2*y1*y2*z1*z2 + y1*y3*z1*z2 + y2*y3*z1*z2 + y1*y4*z1*z2 + y2*y4*z1*z2 - 2*y3*y4*z1*z2 + pow(x1,2)*pow(z2,2) - x1*x4*pow(z2,2) + pow(y1,2)*pow(z2,2) - y1*y3*pow(z2,2) - y1*y4*pow(z2,2) + y3*y4*pow(z2,2) - x1*x4*z1*z3 + y1*y2*z1*z3 - pow(y2,2)*z1*z3 - y1*y4*z1*z3 + y2*y4*z1*z3 - pow(x1,2)*z2*z3 + x1*x4*z2*z3 - pow(y1,2)*z2*z3 + y1*y2*z2*z3 + y1*y4*z2*z3 - y2*y4*z2*z3 + x3*(x4*(pow(y1,2) - 2*y1*y2 + pow(y2,2) + pow(z1 - z2,2)) + x1*(-pow(y2,2) + y1*(y2 - y4) + y2*y4 + (z1 - z2)*(z2 - z4))) + pow(x2,2)*(pow(y1,2) + y3*y4 - y1*(y3 + y4) + (z1 - z3)*(z1 - z4)) + y1*y2*z1*z4 - pow(y2,2)*z1*z4 - y1*y3*z1*z4 + y2*y3*z1*z4 - pow(x1,2)*z2*z4 - pow(y1,2)*z2*z4 + y1*y2*z2*z4 + y1*y3*z2*z4 - y2*y3*z2*z4 + pow(x1,2)*z3*z4 + pow(y1,2)*z3*z4 - 2*y1*y2*z3*z4 + pow(y2,2)*z3*z4 + x2*(-(x4*(pow(y1,2) + y2*y3 - y1*(y2 + y3) + (z1 - z2)*(z1 - z3))) - x3*(pow(y1,2) + y2*y4 - y1*(y2 + y4) + (z1 - z2)*(z1 - z4)) + x1*(-2*y3*y4 + y2*(y3 + y4) + y1*(-2*y2 + y3 + y4) - 2*z1*z2 + z1*z3 + z2*z3 + z1*z4 + z2*z4 - 2*z3*z4))))/(24.*(-(x2*y3*z1) + x2*y4*z1 + x1*y3*z2 - x1*y4*z2 + x2*y1*z3 - x1*y2*z3 + x1*y4*z3 - x2*y4*z3 + x4*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) - x2*y1*z4 + x1*y2*z4 - x1*y3*z4 + x2*y3*z4 + x3*(-(y4*z1) - y1*z2 + y4*z2 + y2*(z1 - z4) + y1*z4)));

  //=====================================================================================
  //row4
  //=====================================================================================

 locmat[3][0] = locmat[0][3];

 locmat[3][1] = locmat[1][3];

 locmat[3][2] = locmat[2][3];

locmat[3][3] = ((k1 + k2 + k3 + k4)*(pow(x1,2)*pow(y2,2) - 2*pow(x1,2)*y2*y3 + pow(x1,2)*pow(y3,2) + pow(y2,2)*pow(z1,2) - 2*y2*y3*pow(z1,2) + pow(y3,2)*pow(z1,2) + pow(x3,2)*(pow(y1 - y2,2) + pow(z1 - z2,2)) - 2*y1*y2*z1*z2 + 2*y1*y3*z1*z2 + 2*y2*y3*z1*z2 - 2*pow(y3,2)*z1*z2 + pow(x1,2)*pow(z2,2) + pow(y1,2)*pow(z2,2) - 2*y1*y3*pow(z2,2) + pow(y3,2)*pow(z2,2) - 2*x2*(x3*((y1 - y2)*(y1 - y3) + (z1 - z2)*(z1 - z3)) + x1*((y1 - y3)*(y2 - y3) + (z1 - z3)*(z2 - z3))) + pow(x2,2)*(pow(y1 - y3,2) + pow(z1 - z3,2)) + 2*x1*x3*((y1 - y2)*(y2 - y3) + (z1 - z2)*(z2 - z3)) - 2*(pow(x1,2)*z2 + (y1 - y2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2))*z3 + (pow(x1,2) + pow(y1 - y2,2))*pow(z3,2)))/(24.*(-(x2*y3*z1) + x2*y4*z1 + x1*y3*z2 - x1*y4*z2 + x2*y1*z3 - x1*y2*z3 + x1*y4*z3 - x2*y4*z3 + x4*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) + (-(x2*y1) + x1*y2 - x1*y3 + x2*y3)*z4 + x3*(-(y4*z1) - y1*z2 + y4*z2 + y2*(z1 - z4) + y1*z4)));

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
  
}

//==============================================================================

void LinearFEM3D_p::ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  ;
}

//==============================================================================

void LinearFEM3D_p::ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalForceCoeff(ind[j]);
 
  double f1 = locdat[0];
  double f2 = locdat[1];
  double f3 = locdat[2];
  double f4 = locdat[3];

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
   double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
  x1 = coord[0][0];
  x2 = coord[0][1];
  x3 = coord[0][2];
  x4 = coord[0][3];
  y1 = coord[1][0];
  y2 = coord[1][1];
  y3 = coord[1][2];
  y4 = coord[1][3];
  z1 = coord[2][0];
  z2 = coord[2][1];
  z3 = coord[2][2];
  z4 = coord[2][3];

locrhs[0] = ((2*f1 + f2 + f3 + f4)*(-(x2*y3*z1) + x2*y4*z1 + x1*y3*z2 - x1*y4*z2 + x2*y1*z3 - x1*y2*z3 + x1*y4*z3 - x2*y4*z3 + x4*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) - x2*y1*z4 + x1*y2*z4 - x1*y3*z4 + x2*y3*z4 + x3*(-(y4*z1) - y1*z2 + y4*z2 + y2*(z1 - z4) + y1*z4)))/120.;

locrhs[1] = ((f1 + 2*f2 + f3 + f4)*(-(x2*y3*z1) + x2*y4*z1 + x1*y3*z2 - x1*y4*z2 + x2*y1*z3 - x1*y2*z3 + x1*y4*z3 - x2*y4*z3 + x4*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) - x2*y1*z4 + x1*y2*z4 - x1*y3*z4 + x2*y3*z4 + x3*(-(y4*z1) - y1*z2 + y4*z2 + y2*(z1 - z4) + y1*z4)))/120.;

locrhs[2] = ((f1 + f2 + 2*f3 + f4)*(-(x2*y3*z1) + x2*y4*z1 + x1*y3*z2 - x1*y4*z2 + x2*y1*z3 - x1*y2*z3 + x1*y4*z3 - x2*y4*z3 + x4*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) - x2*y1*z4 + x1*y2*z4 - x1*y3*z4 + x2*y3*z4 + x3*(-(y4*z1) - y1*z2 + y4*z2 + y2*(z1 - z4) + y1*z4)))/120.;

locrhs[3] = ((f1 + f2 + f3 + 2*f4)*(-(x2*y3*z1) + x2*y4*z1 + x1*y3*z2 - x1*y4*z2 + x2*y1*z3 - x1*y2*z3 + x1*y4*z3 - x2*y4*z3 + x4*(y3*(z1 - z2) + y1*(z2 - z3) + y2*(-z1 + z3)) + (x1*(y2 - y3) + x2*(-y1 + y3))*z4 + x3*(y4*(-z1 + z2) + y2*(z1 - z4) + y1*(-z2 + z4))))/120.;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}


//==============================================================================

void LinearFEM3D_p::ProjectAFunction(double (*func)(Array<double> &), Vector &pval)
{
  int numvert = mesh->GetNV();
  pval.SetSize(numvert);
  for(int j=0; j<numvert; j++)
    {
      Array<double> coord(3);
      for(int i=0; i<3; i++)
	{
	  coord[i] = mesh->GetVertex(j, i);
	}
      pval(j) = func(coord);  
    }
}

//=============================================================================
void  LinearFEM3D_p::ComputeNeumannLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    {
      locdat[j] = data->GetNodalNeumannCoeff(ind[j]); 
    }
 
  double g0 = locdat[0];
  double g1 = locdat[1];
  double g2 = locdat[2];
  
  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    {
      coord[k] = new double[GetBdrElement(i)->GetNVertices()];
    }

  mesh->GetBdrElementVerticesCoord(i, coord);
  double x0,x1,x2,y0,y1,y2,z1,z2,z0;
  x0 = coord[0][0];
  x1 = coord[0][1];
  x2 = coord[0][2];
  y0 = coord[1][0];
  y1 = coord[1][1];
  y2 = coord[1][2];
  z0 = coord[2][0];
  z1 = coord[2][1];
  z2 = coord[2][2];

  //compute norm of cross product of vectors which define plane containing the triangle
  double CPi = (y1-y0)*(z2-z0) - (y2-y0)*(z1-z0);
  double CPj = -((x1-x0)*(z2-z0) - (z1-z0)*(x2-x0));
  double CPk = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0);
  double coeff = sqrt(CPi*CPi + CPj*CPj + CPk*CPk);

  locrhs[0] = coeff * (g0/12.0 + g1/24.0 + g2/24.0);

  locrhs[1] = coeff * (g0/24.0 + g1/12.0 + g2/24.0);

  locrhs[2] = coeff * (g0/24.0 + g1/24.0 + g2/12.0); 

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }

}
