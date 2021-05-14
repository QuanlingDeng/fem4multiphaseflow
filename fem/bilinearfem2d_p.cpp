#include "mpi_fem_header.h"

//==============================================================================

BilinearFEM2D_p::BilinearFEM2D_p(Mesh_p *_mesh, Data *_data) :
FEM_p(_mesh, _data)
{
    strcpy(Type, "LinearFEM1D");
    NumGlobalDOF = mesh->GetNV();
    NumLocalDOF = 4;
    rhs.SetSize(NumGlobalDOF);
    mat = new SparseMatrix(NumGlobalDOF, 10);
}

//============================================================================== 

void BilinearFEM2D_p::GetElementDOF(int i, Array<int> &ind)
{
  GetElementVertices(i, ind);
}

//============================================================================== 

void BilinearFEM2D_p::GetBdrElementDOF(int i, Array<int> &ind)
{
  ind.SetSize(2);
  GetBdrElementVertices(i, ind);
}

//==============================================================================
void BilinearFEM2D_p::ProjectAFunction(double (*func)(Array<double> &), Vector &pval)
{
    pval.SetSize(NumGlobalDOF);
    Array<double> coord(2);
    for(int i=0; i<NumGlobalDOF; i++)
    {
        coord[0] = mesh->GetVertex(i, 0);
        coord[1] = mesh->GetVertex(i, 1);
        
        pval(i) = func(coord);
    }
}



//==============================================================================
void BilinearFEM2D_p::ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalEllipticCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
    double x1 = coord[0][0];
    double y1 = coord[1][0];
    
    double x2 = coord[0][1];
    double y2 = coord[1][2];
    
    double lenx = x2 - x1;
    double leny = y2 - y1;
    
    double k1 = locdat[0];
    double k2 = locdat[1];
    double k3 = locdat[2];
    double k4 = locdat[3];
    
    
    locmat[0][0] = 1.0/(24*lenx*leny)*( (-k4)*(-3*lenx*lenx-leny*leny) + 3*k1*(lenx*lenx+leny*leny) + k3*(lenx*lenx+leny*leny) + k2*(lenx*lenx+3*leny*leny) );
    locmat[0][1] = 1.0/(24*lenx*leny)*( k1*(lenx*lenx-3*leny*leny) + k2*(lenx*lenx-3*leny*leny) + (k3+k4)*(-lenx-leny)*(-lenx+leny)   );
    locmat[0][2] = -1.0/(24*lenx*leny)*(k1+k2+k3+k4)*(lenx*lenx+leny*leny);
    locmat[0][3] = -1.0/(24*lenx*leny)*( k1*(3*lenx*lenx-leny*leny) + k4*(3*lenx*lenx-leny*leny) - (k2+k3)*(-lenx*lenx+leny*leny)    );
    locmat[1][0] = locmat[0][1];
    locmat[1][1] = 1.0/(24*lenx*leny)*(  3*(k2+k3)*lenx*lenx + k4*(lenx*lenx+leny*leny) + k1*(lenx*lenx+3*leny*leny) + (3*k2+k3)*leny*leny  );
    locmat[1][2] = 1.0/(24*lenx*leny)*(  -4*(k1+k4)*lenx*lenx + 3*(k1-k2-k3+k4)*lenx*lenx + (-k1+k2+k3-k4)*leny*leny + 2*(k1+k4)*leny*leny   );
    locmat[1][3] = -1.0/(24*lenx*leny)*(k1+k2+k3+k4)*(lenx*lenx+leny*leny);
    locmat[2][0] = locmat[0][2];
    locmat[2][1] = locmat[1][2];
    locmat[2][2] = 1.0/(24*lenx*leny)*(  3*(k2+k3)*lenx*lenx + k1*(lenx*lenx+leny*leny) + k4*(lenx*lenx+3*leny*leny) + (k2+3*k3)*leny*leny   );
    locmat[2][3] = 1.0/(24*lenx*leny)*(  (k3+k4)*(lenx*lenx-3*leny*leny) + k1*(-lenx-leny)*(-lenx+leny) + k2*(-lenx-leny)*(-lenx+leny)   );
    locmat[3][0] = locmat[0][3];
    locmat[3][1] = locmat[1][3];
    locmat[3][2] = locmat[2][3];
    locmat[3][3] = 1.0/(24*lenx*leny)*(  -k1*(-3*lenx*lenx-leny*leny) + k2*(lenx*lenx+leny*leny) + 3*k4*(lenx*lenx+leny*leny) + k3*(lenx*lenx+3*leny*leny)   );
    
    

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void BilinearFEM2D_p::ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
    Array<double> locdat(ind.Size());
    for (int j=0; j<ind.Size(); j++)
        locdat[j] = data->GetNodalReactionCoeff(ind[j]);
    
    Array<double *> coord(mesh->GetDim());
    for (int k=0; k<coord.Size(); k++)
        coord[k] = new double[GetElement(i)->GetNVertices()];
    
    mesh->GetElementVerticesCoord(i, coord);
    
    double x1 = coord[0][0];
    double y1 = coord[1][0];
    
    double x2 = coord[0][1];
    double y2 = coord[1][2];
    
    double lenx = x2 - x1;
    double leny = y2 - y1;
    
    double r1 = locdat[0];
    double r2 = locdat[1];
    double r3 = locdat[2];
    double r4 = locdat[3];
    
    
    locmat[0][0] = 1.0/144.0 *(9 *r1 + 3 *r2 + r3 + 3 *r4) *lenx*leny;
    locmat[0][1] = 1.0/144.0 *(3 *r1 + 3 *r2 + r3 + r4) *lenx*leny;
    locmat[0][2] = 1.0/144.0 *(r1 + r2 + r3 + r4) *lenx*leny;
    locmat[0][3] = 1.0/144.0 *(3 *r1 + r2 + r3 + 3* r4) *lenx*leny;
    locmat[1][0] = locmat[0][1];
    locmat[1][1] = 1.0/144.0 *(3 *r1 + 9* r2 + 3 *r3 + r4) *lenx*leny;
    locmat[1][2] = 1.0/144.0* (r1 + 3* r2 + 3* r3 + r4) *lenx*leny;
    locmat[1][3] = 1.0/144.0 *(r1 + r2 + r3 + r4) *lenx*leny;
    locmat[2][0] = locmat[0][2];
    locmat[2][1] = locmat[1][2];
    locmat[2][2] = 1.0/144.0 *(r1 + 3* (r2 + 3* r3 + r4)) *lenx*leny;
    locmat[2][3] = 1.0/144.0 *(r1 + r2 + 3 *(r3 + r4)) *lenx*leny;
    locmat[3][0] = locmat[0][3];
    locmat[3][1] = locmat[1][3];
    locmat[3][2] = locmat[2][3];
    locmat[3][3] = 1.0/144.0 *(3* r1 + r2 + 3* r3 + 9* r4)*lenx*leny;
    
    
    
    for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
  
}

//==============================================================================

void BilinearFEM2D_p::ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalForceCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  
    double x1 = coord[0][0];
    double y1 = coord[1][0];
    
    double x2 = coord[0][1];
    double y2 = coord[1][2];
    
    double lenx = x2 - x1;
    double leny = y2 - y1;
    
    double f1 = locdat[0];
    double f2 = locdat[1];
    double f3 = locdat[2];
    double f4 = locdat[3];

    
  locrhs[0] = 1.0/36.0*(4*f1+2*f2+f3+2*f4)*lenx*leny;
  locrhs[1] = 1.0/36.0*(2*f1+4*f2+2*f3+f4)*lenx*leny;
  locrhs[2] = 1.0/36.0*(f1+2*f2+4*f3+2*f4)*lenx*leny;
  locrhs[3] = 1.0/36.0*(2*f1+f2+2*f3+4*f4)*lenx*leny;

  
  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}


