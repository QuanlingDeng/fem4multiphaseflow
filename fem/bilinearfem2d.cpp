#include "fem_header.h"

//==============================================================================

BilinearFEM2D::BilinearFEM2D(Mesh *_mesh, Data *_data) :
FEM(_mesh, _data)
{
    strcpy(Type, "BilinearFEM2D");
    NumGlobalDOF = mesh->GetNV();
    NumLocalDOF = 4;  
    rhs.SetSize(NumGlobalDOF);
    mat = new SparseMatrix(NumGlobalDOF, 10);
}

//============================================================================== 

void BilinearFEM2D::GetElementDOF(int i, Array<int> &ind)
{
  GetElementVertices(i, ind);
}

//============================================================================== 

void BilinearFEM2D::GetBdrElementDOF(int i, Array<int> &ind)
{
  ind.SetSize(2);
  GetBdrElementVertices(i, ind);
}

//==============================================================================
void BilinearFEM2D::ProjectAFunction(double (*func)(Array<double> &), Vector &pval)
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
void BilinearFEM2D::ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    {
      locdat[j] = data->GetNodalEllipticCoeff(ind[j]);
    }

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    {
      coord[k] = new double[GetElement(i)->GetNVertices()];
    }
  mesh->GetElementVerticesCoord(i, coord);
 
  double x1 = coord[0][0];
  double y1 = coord[1][0];
    
  double x2 = coord[0][1];
  double y2 = coord[1][1];

  double x3 = coord[0][3];
  double y3 = coord[1][3];

  double lenx = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
  double leny = sqrt((x3 - x1)*(x3 - x1) + (y3 - y1)*(y3 - y1));
    
  double k1 = locdat[0];
  double k2 = locdat[1];
  double k3 = locdat[2];
  double k4 = locdat[3];
     
  locmat[0][0] = 1.0/(24.0*lenx*leny)*( (-k4)*(-3.0*lenx*lenx-leny*leny) + 3.0*k1*(lenx*lenx+leny*leny) + k3*(lenx*lenx+leny*leny) + k2*(lenx*lenx+3.0*leny*leny) );
  locmat[0][1] = 1.0/(24.0*lenx*leny)*( k1*(lenx*lenx-3.0*leny*leny) + k2*(lenx*lenx-3.0*leny*leny) + (k3+k4)*(-lenx-leny)*(-lenx+leny)   );
  locmat[0][2] = -1.0/(24.0*lenx*leny)*(k1+k2+k3+k4)*(lenx*lenx+leny*leny);
  locmat[0][3] = -1.0/(24.0*lenx*leny)*( k1*(3.0*lenx*lenx-leny*leny) + k4*(3.0*lenx*lenx-leny*leny) - (k2+k3)*(-lenx*lenx+leny*leny)    );
  locmat[1][0] = locmat[0][1];
  locmat[1][1] = 1.0/(24.0*lenx*leny)*(  3.0*(k2+k3)*lenx*lenx + k4*(lenx*lenx+leny*leny) + k1*(lenx*lenx+3.0*leny*leny) + (3.0*k2+k3)*leny*leny  );
  locmat[1][2] = 1.0/(24.0*lenx*leny)*(  -4.0*(k1+k4)*lenx*lenx + 3.0*(k1-k2-k3+k4)*lenx*lenx + (-k1+k2+k3-k4)*leny*leny + 2.0*(k1+k4)*leny*leny   );
  locmat[1][3] = -1.0/(24.0*lenx*leny)*(k1+k2+k3+k4)*(lenx*lenx+leny*leny);
  locmat[2][0] = locmat[0][2];
  locmat[2][1] = locmat[1][2];
  locmat[2][2] = 1.0/(24.0*lenx*leny)*(  3*(k2+k3)*lenx*lenx + k1*(lenx*lenx+leny*leny) + k4*(lenx*lenx+3*leny*leny) + (k2+3*k3)*leny*leny   );
  locmat[2][3] = 1.0/(24.0*lenx*leny)*(  (k3+k4)*(lenx*lenx-3.0*leny*leny) + k1*(-lenx-leny)*(-lenx+leny) + k2*(-lenx-leny)*(-lenx+leny)   );
  locmat[3][0] = locmat[0][3];
  locmat[3][1] = locmat[1][3];
  locmat[3][2] = locmat[2][3];
  locmat[3][3] = 1.0/(24.0*lenx*leny)*(  -k1*(-3.0*lenx*lenx-leny*leny) + k2*(lenx*lenx+leny*leny) + 3.0*k4*(lenx*lenx+leny*leny) + k3*(lenx*lenx+3.0*leny*leny) );
    
  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void BilinearFEM2D::ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
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
        
    locmat[0][0] = 1.0/144.0 *(9.0*r1 + 3*r2 + r3 + 3 *r4)*lenx*leny;
    locmat[0][1] = 1.0/144.0 *(3.0*r1 + 3*r2 + r3 + r4)*lenx*leny;
    locmat[0][2] = 1.0/144.0 *(r1 + r2 + r3 + r4)*lenx*leny;
    locmat[0][3] = 1.0/144.0 *(3.0*r1 + r2 + r3 + 3* r4)*lenx*leny;
    locmat[1][0] = locmat[0][1];
    locmat[1][1] = 1.0/144.0 *(3.0*r1 + 9.0* r2 + 3 *r3 + r4)*lenx*leny;
    locmat[1][2] = 1.0/144.0* (r1 + 3.0*r2 + 3.0*r3 + r4)*lenx*leny;
    locmat[1][3] = 1.0/144.0 *(r1 + r2 + r3 + r4)*lenx*leny;
    locmat[2][0] = locmat[0][2];
    locmat[2][1] = locmat[1][2];
    locmat[2][2] = 1.0/144.0 *(r1 + 3.0* (r2 + 3.0* r3 + r4))*lenx*leny;
    locmat[2][3] = 1.0/144.0 *(r1 + r2 + 3.0 *(r3 + r4))*lenx*leny;
    locmat[3][0] = locmat[0][3];
    locmat[3][1] = locmat[1][3];
    locmat[3][2] = locmat[2][3];
    locmat[3][3] = 1.0/144.0 *(3.0* r1 + r2 + 3.0*r3 + 9.0*r4)*lenx*leny;
    
    for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void BilinearFEM2D::ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
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

//==============================================================================
void BilinearFEM2D::SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
			       const Array<double *> &Dval)
{
  Element *el;  
  for(int j=0; j<bdry_is_dirichlet.Size(); j++)
    {
      for (int i=0; i<GetNBE(); i++)
	{
	  el = GetBdrElement(i);
	  int attr = el->GetAttribute();
	  int nv = el->GetNVertices();

	  if (bdry_is_dirichlet[j][attr])
	    {
	      Array<int> row;
	      GetBdrElementDOF(i, row);

	      Array<int> index(nv);
	      for(int kk=0; kk<nv; kk++) { index[kk] = row[kk]; }
	      int ind = 0;
	      for (int ii=nv*j; ii<nv*j+nv; ii++)
		{
		  mat->EliminateRowCol(row[ii], Dval[j][index[ind]], rhs);
		  ind++;
		}
	    }
	}
    }

}

//==============================================================================
void BilinearFEM2D::ComputeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  Array<double *> locdat(mesh->GetDim());
  for (int k=0; k<locdat.Size(); k++)
    locdat[k] = new double[GetElement(i)->GetNVertices()];

  //set v nodal values
  for (int j=0; j<ind.Size(); j++)
    {
      locdat[0][j] = data->GetNodalAdvectionCoeff(ind[j], 0);
      locdat[1][j] = data->GetNodalAdvectionCoeff(ind[j], 1);
    }
  
  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    {
      coord[k] = new double[GetElement(i)->GetNVertices()];
    }
  mesh->GetElementVerticesCoord(i, coord);
 
  double x1 = coord[0][0];
  double y1 = coord[1][0];
    
  double x2 = coord[0][1];
  double y2 = coord[1][1];

  double x3 = coord[0][3];
  double y3 = coord[1][3];

  double dx = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
  double dy = sqrt((x3 - x1)*(x3 - x1) + (y3 - y1)*(y3 - y1));

  locmat[0][0] = (-((6*locdat[1][0] + 2*locdat[1][1] + locdat[1][2] + 3*locdat[1][3])*dx) - (6*locdat[0][0] + 3*locdat[0][1] + locdat[0][2] + 2*locdat[0][3])*dy)/72.;
  locmat[0][1] = (-((2*locdat[1][0] + 2*locdat[1][1] + locdat[1][2] + locdat[1][3])*dx) + (6*locdat[0][0] + 3*locdat[0][1] + locdat[0][2] + 2*locdat[0][3])*dy)/72.;
  locmat[0][2] = ((2*locdat[1][0] + 2*locdat[1][1] + locdat[1][2] + locdat[1][3])*dx + (2*locdat[0][0] + locdat[0][1] + locdat[0][2] + 2*locdat[0][3])*dy)/72.;
  locmat[0][3] = ((6*locdat[1][0] + 2*locdat[1][1] + locdat[1][2] + 3*locdat[1][3])*dx - (2*locdat[0][0] + locdat[0][1] + locdat[0][2] + 2*locdat[0][3])*dy)/72.;
 
  locmat[1][0] = (-((2*locdat[1][0] + 2*locdat[1][1] + locdat[1][2] + locdat[1][3])*dx) - (3*locdat[0][0] + 6*locdat[0][1] + 2*locdat[0][2] + locdat[0][3])*dy)/72.;
  locmat[1][1] = (-((2*locdat[1][0] + 6*locdat[1][1] + 3*locdat[1][2] + locdat[1][3])*dx) + (3*locdat[0][0] + 6*locdat[0][1] + 2*locdat[0][2] + locdat[0][3])*dy)/72.; 
  locmat[1][2] = ((2*locdat[1][0] + 6*locdat[1][1] + 3*locdat[1][2] + locdat[1][3])*dx + (locdat[0][0] + 2*(locdat[0][1] + locdat[0][2]) + locdat[0][3])*dy)/72.;
  locmat[1][3] = ((2*locdat[1][0] + 2*locdat[1][1] + locdat[1][2] + locdat[1][3])*dx - (locdat[0][0] + 2*(locdat[0][1] + locdat[0][2]) + locdat[0][3])*dy)/72.;

  locmat[2][0] = (-((locdat[1][0] + locdat[1][1] + 2*(locdat[1][2] + locdat[1][3]))*dx) - (locdat[0][0] + 2*(locdat[0][1] + locdat[0][2]) + locdat[0][3])*dy)/72.;
  locmat[2][1] = (-((locdat[1][0] + 3*locdat[1][1] + 6*locdat[1][2] + 2*locdat[1][3])*dx) + (locdat[0][0] + 2*(locdat[0][1] + locdat[0][2]) + locdat[0][3])*dy)/72.;
  locmat[2][2] = ((locdat[1][0] + 3*locdat[1][1] + 6*locdat[1][2] + 2*locdat[1][3])*dx + (locdat[0][0] + 2*locdat[0][1] + 6*locdat[0][2] + 3*locdat[0][3])*dy)/72.;
  locmat[2][3] = ((locdat[1][0] + locdat[1][1] + 2*(locdat[1][2] + locdat[1][3]))*dx - (locdat[0][0] + 2*locdat[0][1] + 6*locdat[0][2] + 3*locdat[0][3])*dy)/72.;

  locmat[3][0] = (-((3*locdat[1][0] + locdat[1][1] + 2*locdat[1][2] + 6*locdat[1][3])*dx) - (2*locdat[0][0] + locdat[0][1] + locdat[0][2] + 2*locdat[0][3])*dy)/72.;
  locmat[3][1] = (-((locdat[1][0] + locdat[1][1] + 2*(locdat[1][2] + locdat[1][3]))*dx) + (2*locdat[0][0] + locdat[0][1] + locdat[0][2] + 2*locdat[0][3])*dy)/72.;
  locmat[3][2] = ((locdat[1][0] + locdat[1][1] + 2*(locdat[1][2] + locdat[1][3]))*dx + (2*locdat[0][0] + locdat[0][1] + 3*locdat[0][2] + 6*locdat[0][3])*dy)/72.;
  locmat[3][3] = ((3*locdat[1][0] + locdat[1][1] + 2*locdat[1][2] + 6*locdat[1][3])*dx - (2*locdat[0][0] + locdat[0][1] + 3*locdat[0][2] + 6*locdat[0][3])*dy)/72.;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
  for (int k=0; k<locdat.Size(); k++) { delete []locdat[k]; }
}

//==============================================================================

void BilinearFEM2D::ComputeConservativeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
}

//==============================================================================

void BilinearFEM2D::ComputeNeumannLocalSystem(int i, Array<int> &ind, 
					    Array<double> &locrhs)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    {
     locdat[j] = data->GetNodalNeumannCoeff(ind[j]);
    }

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetBdrElement(i)->GetNVertices()];

  mesh->GetBdrElementVerticesCoord(i, coord);

  double h = (coord[0][1] - coord[0][0]) * (coord[0][1] - coord[0][0]) 
                + (coord[1][1] - coord[1][0]) * (coord[1][1] - coord[1][0]);

  h = sqrt(h);
  locrhs[0] = -h * ( 2.0*locdat[0] + locdat[1] ) / 6.0;
  locrhs[1] = -h * ( locdat[0] + 2.0*locdat[1] ) / 6.0;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void BilinearFEM2D::ComputeNeumannLocalSystem(int b, int i, Array<int> &ind, 
					    Array<double> &locrhs)
{

  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    {
      locdat[j] = data->GetNeumannBdryVal(b, ind[j]);
    }

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetBdrElement(i)->GetNVertices()];

  mesh->GetBdrElementVerticesCoord(i, coord);

  double h = (coord[0][1] - coord[0][0]) * (coord[0][1] - coord[0][0]) 
                + (coord[1][1] - coord[1][0]) * (coord[1][1] - coord[1][0]);

  h = sqrt(h);

  //cout << h << endl;
  locrhs[0] = -h * ( 2.0*locdat[0] + locdat[1] ) / 6.0;
  locrhs[1] = -h * ( locdat[0] + 2.0*locdat[1] ) / 6.0;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//============================================================================

void BilinearFEM2D::ComputeRobinLocalSystem(int i, Array<int> &ind, 
					  Array<double *> &locmat, 
					  Array<double> &locrhs)
{
  double g0,g1;
  g0 = data->GetNodalRobinCoeff(ind[0], 0);
  g1 = data->GetNodalRobinCoeff(ind[1], 0);

  double n0,n1;
  n0 = data->GetNodalRobinCoeff(ind[0], 1);
  n1 = data->GetNodalRobinCoeff(ind[1], 1);

  Array<double *> coord(2);
  coord[0] = new double[2];
  coord[1] = new double[2];
  
  mesh->GetBdrElementVerticesCoord(i, coord);

  double dx, dy, leng;
  dx = coord[0][0] - coord[0][1];
  dy = coord[1][0] - coord[1][1];

  leng = sqrt(pow(dx,2) + pow(dy,2));

  locmat[0][0] = leng * (g1 + 3.0*g0)/12.0;
  locmat[0][1] = leng * (g0 + g1)/12.0;

  locmat[1][0] = leng * (g0 + g1)/12.0;
  locmat[1][1] = leng * (3.0*g1 + g0)/12.0;
  
  locrhs[0] = -leng * ( 1.0/6.0 * n1 + 1.0/3.0 * n0 );
  locrhs[1] = -leng * ( 1.0/3.0 * n1 + 1.0/6.0 * n0 );

  delete []coord[0];
  delete []coord[1];
}

//==============================================================================

void BilinearFEM2D::ComputeRobinLocalSystem(int i, int b, Array<int> &ind, 
					    Array<double *> &locmat, Array<double> &locrhs)
{
  double g0,g1;
  g0 = data->GetRobinCoeff(b, ind[0]);
  g1 = data->GetRobinCoeff(b, ind[1]);

  double n0,n1;
  n0 = data->GetRobinBdryVal(b, ind[0]);
  n1 = data->GetRobinBdryVal(b, ind[1]);

  Array<double *> coord(2);
  coord[0] = new double[2];
  coord[1] = new double[2];
  
  mesh->GetBdrElementVerticesCoord(i, coord);

  double dx, dy, leng;
  dx = coord[0][0] - coord[0][1];
  dy = coord[1][0] - coord[1][1];

  leng = sqrt(pow(dx,2) + pow(dy,2));

  //cout << leng << endl;

  locmat[0][0] = leng * (g1 + 3.0*g0)/12.0;
  locmat[0][1] = leng * (g0 + g1)/12.0;

  locmat[1][0] = leng * (g0 + g1)/12.0;
  locmat[1][1] = leng * (3.0*g1 + g0)/12.0;
  
  locrhs[0] = -leng * ( 1.0/6.0 * n1 + 1.0/3.0 * n0 );
  locrhs[1] = -leng * ( 1.0/3.0 * n1 + 1.0/6.0 * n0 );

  delete []coord[0];
  delete []coord[1];
}

//==============================================================================

double BilinearFEM2D::ComputeL2Error(double (*solfunc)(Array<double> &), Vector &approxsol)
{ 
  Array<int> v(2);
  double w[3];

  w[0] = 0.27777777777777777778; w[1] = 0.44444444444444444444; w[2] = w[0];

  int numquadpoints = 9;
  Array<double> quadweight(numquadpoints);
  int ind = 0;
  for (int j=0; j<3; j++) {
    for (int i=0; i<3; i++) {
      quadweight[ind] = w[i] * w[j];
      ind++;
    }
  }

  Array<double> xi(3), eta(3);
  eta[0] = xi[0] = 0.11270166537925831148;
  eta[1] = xi[1] = 0.5;
  eta[2] = xi[2] = 1.0-xi[0];

  Array<double *> quadbasisvalues(numquadpoints);
  for (int i=0; i<numquadpoints; i++) { quadbasisvalues[i] = new double[4]; }

  ind = 0;
  for (int j=0; j<3; j++)
    {
      for (int i=0; i<3; i++) 
	{
	  double valx = 1.0 - xi[i];
	  double valy = 1.0 - eta[j];
	  quadbasisvalues[ind][0] = valx * valy;
	  quadbasisvalues[ind][1] = xi[i] * valy;
	  quadbasisvalues[ind][2] = xi[i] * eta[j];
	  quadbasisvalues[ind][3] = valx * eta[j];
	  ind++;
	}
    }
  
  double val = 0.0;
  for (int i=0; i<GetNE(); i++)
    {
      double tempval = 0.0;
      GetElementVertices(i, v);
      
      Array<double *> coords(mesh->GetDim());
      for (int k=0; k<coords.Size(); k++){ coords[k] = new double[GetElement(i)->GetNVertices()]; }
      mesh->GetElementVerticesCoord(i, coords);
	 
      double Hx = coords[0][1] - coords[0][0];
      double Hy = coords[1][3] - coords[1][0];

      ind = 0;
      Array<double> coord(2);
      for(int k=0; k<3; k++)
	{
	  coord[1] = coords[1][0] + eta[k]*Hy ; // (j + eta[k]) * Hy;

	  for (int l=0; l<3; l++)
	    {
	      coord[0] = coords[0][0] + xi[l]*Hx; // (i + xi[l]) * Hx;

	      double diff = 0.0; 
	      for (int ii=0; ii<4; ii++) 
		{ 
		  diff += approxsol(v[ii]) * quadbasisvalues[ind][ii]; 
		}

	      diff -= solfunc(coord);

	      tempval += diff*diff*quadweight[ind];
	      ind++;
	    }
	}	      

      tempval *= Hx*Hy;
	  
      val += tempval;

      for (int k=0; k<coords.Size(); k++){ delete []coords[k]; }
    }

  for (int i=0; i<9; i++) { delete []quadbasisvalues[i]; }

  return sqrt(val);
  //cout << "The L2 norm error is: " << sqrt(val) << endl;
}

//==============================================================================

double BilinearFEM2D::ComputeH1Error(double (*deronefunc)(Array<double> &), double (*dertwofunc)(Array<double> &), Vector &approxsol)
{ 
  Array<int> v(4);   
  double w[3];

  w[0] = 0.27777777777777777778; w[1] = 0.44444444444444444444; w[2] = w[0];

  int numquadpoints = 9;
  Array<double> quadweight(numquadpoints);
  int ind = 0;
  for (int j=0; j<3; j++) {
    for (int i=0; i<3; i++) {
      quadweight[ind] = w[i] * w[j];
      ind++;
    }
  }

  Array<double> xi(3), eta(3);
  eta[0] = xi[0] = 0.11270166537925831148;
  eta[1] = xi[1] = 0.5;
  eta[2] = xi[2] = 1-xi[0];
  
  Array<RectangularMatrix *> quadbasisvalues(numquadpoints);
  for (int i=0; i<numquadpoints; i++) { quadbasisvalues[i] = new RectangularMatrix(4,2); }

  ind = 0;
  for (int j=0; j<3; j++)
    for (int i=0; i<3; i++) {

      quadbasisvalues[ind]->Elem(0,0) = -1.0 + eta[j]; quadbasisvalues[ind]->Elem(0,1) = -1.0 + xi[i]; 
      quadbasisvalues[ind]->Elem(1,0) =  1.0 - eta[j]; quadbasisvalues[ind]->Elem(1,1) = -xi[i];
      quadbasisvalues[ind]->Elem(2,0) =  eta[j];       quadbasisvalues[ind]->Elem(2,1) = xi[i];
      quadbasisvalues[ind]->Elem(3,0) = -eta[j];       quadbasisvalues[ind]->Elem(3,1) = 1.0 - xi[i];
      ind++;
    }
  
  Array<double> diff(2);
  double val = 0.0;
  
  for (int i=0; i<GetNE(); i++)
    {     
      double tempval = 0.0;
      GetElementVertices(i, v);

      Array<double *> coords(mesh->GetDim());
      for (int k=0; k<coords.Size(); k++)
	coords[k] = new double[GetElement(i)->GetNVertices()];
	 
      mesh->GetElementVerticesCoord(i, coords);
	 
      double Hx = coords[0][1] - coords[0][0];
      double Hy = coords[1][3] - coords[1][0];

      ind = 0;
      Array<double> coord(2);
      for (int k=0; k<3; k++)
	{
	  coord[1] = coords[1][0] + eta[k]*Hy ; // (j + eta[k]) * Hy;

	  for (int l=0; l<3; l++)
	    {
	      coord[0] = coords[0][0] + xi[l]*Hx; // (i + xi[l]) * Hx;

	      diff[0] = diff[1] = 0.0;
	      for (int ii=0; ii<4; ii++)
		for (int jj=0; jj<2; jj++)
		  diff[jj] += approxsol(v[ii]) * quadbasisvalues[ind]->Elem(ii,jj);
	      diff[0] /= Hx; diff[1] /= Hy;
	      diff[0] -= deronefunc(coord);
	      diff[1] -= dertwofunc(coord);

	      tempval += (diff[0]*diff[0] + diff[1]*diff[1])*quadweight[ind];
	      ind++;
	    }
	}	      
     val += tempval*Hx*Hy;	
     for (int k=0; k<coords.Size(); k++){ delete []coords[k]; }
    }

  for (int i=0; i<numquadpoints; i++) { delete quadbasisvalues[i]; }

  return  sqrt(val);
  //cout << "H1 Error: " << sqrt(val) << endl;
}

//==============================================================================

double BilinearFEM2D::ComputeL2Error(Function *exactsol, Vector &approxsol, double time)
{ 
  Array<int> v(2);
  double w[3];

  w[0] = 0.27777777777777777778; w[1] = 0.44444444444444444444; w[2] = w[0];

  int numquadpoints = 9;
  Array<double> quadweight(numquadpoints);
  int ind = 0;
  for (int j=0; j<3; j++) {
    for (int i=0; i<3; i++) {
      quadweight[ind] = w[i] * w[j];
      ind++;
    }
  }

  Array<double> xi(3), eta(3);
  eta[0] = xi[0] = 0.11270166537925831148;
  eta[1] = xi[1] = 0.5;
  eta[2] = xi[2] = 1.0-xi[0];

  Array<double *> quadbasisvalues(numquadpoints);
  for (int i=0; i<numquadpoints; i++) { quadbasisvalues[i] = new double[4]; }

  ind = 0;
  for (int j=0; j<3; j++)
    for (int i=0; i<3; i++) {
      double valx = 1.0 - xi[i];
      double valy = 1.0 - eta[j];
      quadbasisvalues[ind][0] = valx * valy;
      quadbasisvalues[ind][1] = xi[i] * valy;
      quadbasisvalues[ind][2] = xi[i] * eta[j];
      quadbasisvalues[ind][3] = valx * eta[j];
      ind++;
    }
  
  double val = 0.0;
  for (int i=0; i<GetNE(); i++)
    {
     
      double tempval = 0.0;
      GetElementVertices(i, v);

      Array<double *> coords(mesh->GetDim());
      for (int k=0; k<coords.Size(); k++)
	coords[k] = new double[GetElement(i)->GetNVertices()];
	 
      mesh->GetElementVerticesCoord(i, coords);
	 
      double Hx = coords[0][1] - coords[0][0];
      double Hy = coords[1][3] - coords[1][0];

      ind = 0;
      double coord[3];
      coord[2] = time;
      for(int k=0; k<3; k++)
	{
	  coord[1] = coords[1][0] + eta[k]*Hy ; // (j + eta[k]) * Hy;

	  for (int l=0; l<3; l++)
	    {
	      coord[0] = coords[0][0] + xi[l]*Hx; // (i + xi[l]) * Hx;

	      double diff = 0.0; 
	      for (int ii=0; ii<4; ii++) 
		{ 
		  diff += approxsol(v[ii]) * quadbasisvalues[ind][ii]; 
		}

	      diff -= exactsol->Eval(coord);

	      tempval += diff*diff*quadweight[ind];
	      ind++;
	    }
	}	      

      tempval *= Hx*Hy;
	  
      val += tempval;

      for (int k=0; k<coords.Size(); k++){ delete []coords[k]; }
    }

  for (int i=0; i<9; i++) { delete []quadbasisvalues[i]; }

  return sqrt(val);
  //cout << "The L2 error is: " << sqrt(val) << endl;
}

//==============================================================================

double BilinearFEM2D::ComputeH1Error(Array<Function *> &exactderivative, Vector &approxsol)
{
  Array<int> v(4);   
  double w[3];

  w[0] = 0.27777777777777777778; w[1] = 0.44444444444444444444; w[2] = w[0];

  int numquadpoints = 9;
  Array<double> quadweight(numquadpoints);
  int ind = 0;
  for (int j=0; j<3; j++) {
    for (int i=0; i<3; i++) {
      quadweight[ind] = w[i] * w[j];
      ind++;
    }
  }

  Array<double> xi(3), eta(3);
  eta[0] = xi[0] = 0.11270166537925831148;
  eta[1] = xi[1] = 0.5;
  eta[2] = xi[2] = 1-xi[0];
  
  Array<RectangularMatrix *> quadbasisvalues(numquadpoints);
  for (int i=0; i<numquadpoints; i++) { quadbasisvalues[i] = new RectangularMatrix(4,2); }

  ind = 0;
  for (int j=0; j<3; j++)
    for (int i=0; i<3; i++) {

      quadbasisvalues[ind]->Elem(0,0) = -1.0 + eta[j]; quadbasisvalues[ind]->Elem(0,1) = -1.0 + xi[i]; 
      quadbasisvalues[ind]->Elem(1,0) =  1.0 - eta[j]; quadbasisvalues[ind]->Elem(1,1) = -xi[i];
      quadbasisvalues[ind]->Elem(2,0) =  eta[j];       quadbasisvalues[ind]->Elem(2,1) = xi[i];
      quadbasisvalues[ind]->Elem(3,0) = -eta[j];       quadbasisvalues[ind]->Elem(3,1) = 1.0 - xi[i];
      ind++;
    }
  
  Array<double> diff(2);
  double val = 0.0;
  
  for (int i=0; i<GetNE(); i++)
    {     
      double tempval = 0.0;
      GetElementVertices(i, v);

      Array<double *> coords(mesh->GetDim());
      for (int k=0; k<coords.Size(); k++)
	coords[k] = new double[GetElement(i)->GetNVertices()];
	 
      mesh->GetElementVerticesCoord(i, coords);
	 
      double Hx = coords[0][1] - coords[0][0];
      double Hy = coords[1][3] - coords[1][0];

      ind = 0;
      double coord[2];
      for (int k=0; k<3; k++)
	{
	  coord[1] = coords[1][0] + eta[k]*Hy ; // (j + eta[k]) * Hy;

	  for (int l=0; l<3; l++)
	    {
	      coord[0] = coords[0][0] + xi[l]*Hx; // (i + xi[l]) * Hx;

	      diff[0] = diff[1] = 0.0;
	      for (int ii=0; ii<4; ii++)
		for (int jj=0; jj<2; jj++)
		  diff[jj] += approxsol(v[ii]) * quadbasisvalues[ind]->Elem(ii,jj);
	      
	      diff[0] /= Hx; diff[1] /= Hy;
	      diff[0] -= exactderivative[0]->Eval(coord);
	      diff[1] -= exactderivative[1]->Eval(coord);

	      tempval += (diff[0]*diff[0] + diff[1]*diff[1])*quadweight[ind];
	      ind++;
	    }
	}	      
     val += tempval*Hx*Hy;	
     for (int k=0; k<coords.Size(); k++){ delete []coords[k]; }
    }

  for (int i=0; i<numquadpoints; i++) { delete quadbasisvalues[i]; }

  return sqrt(val);
  //cout << "H1 Error: " << sqrt(val) << endl;
}

//==============================================================================

void BilinearFEM2D::ComputeStableLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  ;
}

//==============================================================================

void BilinearFEM2D::ComputeStableRHSLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
{
  ;
}
