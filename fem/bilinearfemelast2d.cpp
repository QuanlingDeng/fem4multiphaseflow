#include "fem_header.h"
//==============================================================================

BilinearFEMELAST2D::BilinearFEMELAST2D(Mesh *_mesh, Data *_data) :
FEM(_mesh, _data)
{
    strcpy(Type, "BilinearFEMELAST2D");
    NumGlobalDOF = 2*mesh->GetNV();
    NumLocalDOF = 8;
    rhs.SetSize(NumGlobalDOF);
    mat = new SparseMatrix(NumGlobalDOF, 25); 
}

//============================================================================= 

void BilinearFEMELAST2D::GetElementDOF(int i, Array<int> &ind)
{
  Array<int> tempind;
  GetElementVertices(i, tempind);
  ind.SetSize(8);  

  for(int i=0; i<4; i++)
    {
      ind[i] = tempind[i];
      ind[i+4] = tempind[i] + mesh->GetNV();
    }
}

//==============================================================================

void BilinearFEMELAST2D::GetBdrElementDOF(int i, Array<int> &ind)
{
  Array<int> tempind;
  GetBdrElementVertices(i, tempind);
  ind.SetSize(4);
  for(int i=0; i<2; i++)
    {
      ind[i] = tempind[i];
      ind[i+2] = tempind[i] + mesh->GetNV();
    }
}


//==============================================================================

void BilinearFEMELAST2D::ComputeEllipticLocalSystem(int i, Array<int> &ind, 
						  Array<double *> &locmat)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<4; j++)
    {
       locdat[j] = data->GetNodalEllipticCoeff(ind[j], 0);
       locdat[j+4] = data->GetNodalEllipticCoeff(ind[j], 1);
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
    
  double l0 = locdat[0];
  double l1 = locdat[1];
  double l2 = locdat[2];
  double l3 = locdat[3];

  double m0 = locdat[4];
  double m1 = locdat[5];
  double m2 = locdat[6];
  double m3 = locdat[7];
     
  locmat[0][0] = (2*pow(dy,2)*(3*l0 + 3*l1 + l2 + l3 + 3*m0 + 3*m1 + m2 + m3) + pow(dx,2)*(3*m0 + m1 + m2 + 3*m3))/(48.*dx*dy);
  locmat[0][1] = (pow(dx,2)*(m0 + m1 + m2 + m3) - 2*pow(dy,2)*(3*l0 + 3*l1 + l2 + l3 + 3*m0 + 3*m1 + m2 + m3))/(48.*dx*dy);
  locmat[0][2] = -(pow(dx,2)*(m0 + m1 + m2 + m3) + 2*pow(dy,2)*(l0 + l1 + l2 + l3 + m0 + m1 + m2 + m3))/(48.*dx*dy);
  locmat[0][3] = (2*pow(dy,2)*(l0 + l1 + l2 + l3 + m0 + m1 + m2 + m3) - pow(dx,2)*(3*m0 + m1 + m2 + 3*m3))/(48.*dx*dy);

  locmat[0][4] = (8*l0 + 4*l1 + 2*l2 + 4*l3 + 4*m0 + 2*m1 + m2 + 2*m3)/72.;
  locmat[0][5] = (4*l0 + 8*l1 + 4*l2 + 2*l3 - 4*m0 - 2*m1 - m2 - 2*m3)/72.;
  locmat[0][6] = (-4*l0 - 8*l1 - 4*l2 - 2*l3 - 2*m0 - m1 - 2*m2 - 4*m3)/72.;
  locmat[0][7] = (-8*l0 - 4*l1 - 2*l2 - 4*l3 + 2*m0 + m1 + 2*m2 + 4*m3)/72.;

  locmat[1][0] = locmat[0][1];
  locmat[1][1] = (2*pow(dy,2)*(3*l0 + 3*l1 + l2 + l3 + 3*m0 + 3*m1 + m2 + m3) + pow(dx,2)*(m0 + 3*(m1 + m2) + m3))/(48.*dx*dy);
  locmat[1][2] = (2*pow(dy,2)*(l0 + l1 + l2 + l3 + m0 + m1 + m2 + m3) - pow(dx,2)*(m0 + 3*(m1 + m2) + m3))/(48.*dx*dy);
  locmat[1][3] = -(pow(dx,2)*(m0 + m1 + m2 + m3) + 2*pow(dy,2)*(l0 + l1 + l2 + l3 + m0 + m1 + m2 + m3))/(48.*dx*dy);

  locmat[1][4] = (-8*l0 - 4*l1 - 2*l2 - 4*l3 + 2*(m0 + 2*m1 + m2) + m3)/72.;
  locmat[1][5] = (-2*(2*l0 + 4*l1 + 2*l2 + l3 + m0 + 2*m1 + m2) - m3)/72.;
  locmat[1][6] = (4*l0 + 8*l1 + 4*l2 + 2*l3 - m0 - 2*(m1 + 2*m2 + m3))/72.;
  locmat[1][7] = (8*l0 + 4*l1 + 2*l2 + 4*l3 + m0 + 2*(m1 + 2*m2 + m3))/72.;

  locmat[2][0] = locmat[0][2];
  locmat[2][1] = locmat[1][2];
  locmat[2][2] = (pow(dx,2)*(m0 + 3*(m1 + m2) + m3) + 2*pow(dy,2)*(l0 + l1 + 3*l2 + 3*l3 + m0 + m1 + 3*(m2 + m3)))/(48.*dx*dy);
  locmat[2][3] = (pow(dx,2)*(m0 + m1 + m2 + m3) - 2*pow(dy,2)*(l0 + l1 + 3*l2 + 3*l3 + m0 + m1 + 3*(m2 + m3)))/(48.*dx*dy);

  locmat[2][4] = (-2*(2*l0 + l1 + 2*l2 + 4*l3 + m0 + 2*m1 + m2) - m3)/72.;
  locmat[2][5] = (-2*l0 - 4*l1 - 8*l2 - 4*l3 + 2*(m0 + 2*m1 + m2) + m3)/72.;
  locmat[2][6] = (2*l0 + 4*l1 + 8*l2 + 4*l3 + m0 + 2*(m1 + 2*m2 + m3))/72.;
  locmat[2][7] = (4*l0 + 2*l1 + 4*l2 + 8*l3 - m0 - 2*(m1 + 2*m2 + m3))/72.;
  
  locmat[3][0] = locmat[0][3];
  locmat[3][1] = locmat[1][3];
  locmat[3][2] = locmat[2][3];
  locmat[3][3] = (pow(dx,2)*(3*m0 + m1 + m2 + 3*m3) + 2*pow(dy,2)*(l0 + l1 + 3*l2 + 3*l3 + m0 + m1 + 3*(m2 + m3)))/(48.*dx*dy);

  locmat[3][4] = (4*l0 + 2*l1 + 4*l2 + 8*l3 - 4*m0 - 2*m1 - m2 - 2*m3)/72.;
  locmat[3][5] = (2*l0 + 4*l1 + 8*l2 + 4*l3 + 4*m0 + 2*m1 + m2 + 2*m3)/72.;
  locmat[3][6] = (-2*l0 - 4*l1 - 8*l2 - 4*l3 + 2*m0 + m1 + 2*m2 + 4*m3)/72.;
  locmat[3][7] = (-4*l0 - 2*l1 - 4*l2 - 8*l3 - 2*m0 - m1 - 2*m2 - 4*m3)/72.;
  
  //===============================================================================================================================
  
  locmat[4][0] = (8*l0 + 4*l1 + 2*l2 + 4*l3 + 4*m0 + 2*m1 + m2 + 2*m3)/72.;
  locmat[4][1] = (-8*l0 - 4*l1 - 2*l2 - 4*l3 + 2*(m0 + 2*m1 + m2) + m3)/72.;
  locmat[4][2] = (-2*(2*l0 + l1 + 2*l2 + 4*l3 + m0 + 2*m1 + m2) - m3)/72.;
  locmat[4][3] = (4*l0 + 2*l1 + 4*l2 + 8*l3 - 4*m0 - 2*m1 - m2 - 2*m3)/72.;

  locmat[4][4] = (pow(dy,2)*(3*m0 + 3*m1 + m2 + m3) + 2*pow(dx,2)*(3*l0 + l1 + l2 + 3*l3 + 3*m0 + m1 + m2 + 3*m3))/(48.*dx*dy);
  locmat[4][5] = (2*pow(dx,2)*(l0 + l1 + l2 + l3 + m0 + m1 + m2 + m3) - pow(dy,2)*(3*m0 + 3*m1 + m2 + m3))/(48.*dx*dy);
  locmat[4][6] = -(pow(dy,2)*(m0 + m1 + m2 + m3) + 2*pow(dx,2)*(l0 + l1 + l2 + l3 + m0 + m1 + m2 + m3))/(48.*dx*dy);
  locmat[4][7] = (pow(dy,2)*(m0 + m1 + m2 + m3) - 2*pow(dx,2)*(3*l0 + l1 + l2 + 3*l3 + 3*m0 + m1 + m2 + 3*m3))/(48.*dx*dy);

  locmat[5][0] = (4*l0 + 8*l1 + 4*l2 + 2*l3 - 4*m0 - 2*m1 - m2 - 2*m3)/72.;
  locmat[5][1] = (-2*(2*l0 + 4*l1 + 2*l2 + l3 + m0 + 2*m1 + m2) - m3)/72.;
  locmat[5][2] = (-2*l0 - 4*l1 - 8*l2 - 4*l3 + 2*(m0 + 2*m1 + m2) + m3)/72.;
  locmat[5][3] = (2*l0 + 4*l1 + 8*l2 + 4*l3 + 4*m0 + 2*m1 + m2 + 2*m3)/72.;

  locmat[5][4] = locmat[4][5]; 
  locmat[5][5] = (pow(dy,2)*(3*m0 + 3*m1 + m2 + m3) + 2*pow(dx,2)*(l0 + 3*l1 + 3*l2 + l3 + m0 + 3*(m1 + m2) + m3))/(48.*dx*dy);
  locmat[5][6] = (pow(dy,2)*(m0 + m1 + m2 + m3) - 2*pow(dx,2)*(l0 + 3*l1 + 3*l2 + l3 + m0 + 3*(m1 + m2) + m3))/(48.*dx*dy);
  locmat[5][7] = -(pow(dy,2)*(m0 + m1 + m2 + m3) + 2*pow(dx,2)*(l0 + l1 + l2 + l3 + m0 + m1 + m2 + m3))/(48.*dx*dy);

  locmat[6][0] = (-4*l0 - 8*l1 - 4*l2 - 2*l3 - 2*m0 - m1 - 2*m2 - 4*m3)/72.;
  locmat[6][1] = (4*l0 + 8*l1 + 4*l2 + 2*l3 - m0 - 2*(m1 + 2*m2 + m3))/72.;
  locmat[6][2] = (2*l0 + 4*l1 + 8*l2 + 4*l3 + m0 + 2*(m1 + 2*m2 + m3))/72.;
  locmat[6][3] = (-2*l0 - 4*l1 - 8*l2 - 4*l3 + 2*m0 + m1 + 2*m2 + 4*m3)/72.;

  locmat[6][4] = locmat[4][6];
  locmat[6][5] = locmat[5][6];
  locmat[6][6] = (2*pow(dx,2)*(l0 + 3*l1 + 3*l2 + l3 + m0 + 3*(m1 + m2) + m3) + pow(dy,2)*(m0 + m1 + 3*(m2 + m3)))/(48.*dx*dy);
  locmat[6][7] = (2*pow(dx,2)*(l0 + l1 + l2 + l3 + m0 + m1 + m2 + m3) - pow(dy,2)*(m0 + m1 + 3*(m2 + m3)))/(48.*dx*dy);

  locmat[7][0] = (-8*l0 - 4*l1 - 2*l2 - 4*l3 + 2*m0 + m1 + 2*m2 + 4*m3)/72.;
  locmat[7][1] = (8*l0 + 4*l1 + 2*l2 + 4*l3 + m0 + 2*(m1 + 2*m2 + m3))/72.;
  locmat[7][2] = (4*l0 + 2*l1 + 4*l2 + 8*l3 - m0 - 2*(m1 + 2*m2 + m3))/72.;
  locmat[7][3] = (-4*l0 - 2*l1 - 4*l2 - 8*l3 - 2*m0 - m1 - 2*m2 - 4*m3)/72.;

  locmat[7][4] = locmat[4][7];
  locmat[7][5] = locmat[5][7];
  locmat[7][6] = locmat[6][7];
  locmat[7][7] = (2*pow(dx,2)*(3*l0 + l1 + l2 + 3*l3 + 3*m0 + m1 + m2 + 3*m3) + pow(dy,2)*(m0 + m1 + 3*(m2 + m3)))/(48.*dx*dy);
    
  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }

}

//=============================================================================

void BilinearFEMELAST2D::ComputeReactionLocalSystem(int i, Array<int> &ind, 
						  Array<double *> &locmat)
{
  ;
}

//=============================================================================

void BilinearFEMELAST2D::ComputeForceLocalSystem(int i, Array<int> &ind, 
					       Array<double> &locrhs)
{
  int indsize = ind.Size();
  Array<double> locdat(indsize);
  for (int j=0; j<4; j++)
    {
      locdat[j] = data->GetNodalForceCoeff(ind[j], 0);
      locdat[j+4] = data->GetNodalForceCoeff(ind[j],1);
    }

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

  f1 = locdat[4];
  f2 = locdat[5];
  f3 = locdat[6];
  f4 = locdat[7];
    
  locrhs[4] = 1.0/36.0*(4*f1+2*f2+f3+2*f4)*lenx*leny;
  locrhs[5] = 1.0/36.0*(2*f1+4*f2+2*f3+f4)*lenx*leny;
  locrhs[6] = 1.0/36.0*(f1+2*f2+4*f3+2*f4)*lenx*leny;
  locrhs[7] = 1.0/36.0*(2*f1+f2+2*f3+4*f4)*lenx*leny;

  
  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void BilinearFEMELAST2D::ComputeNeumannLocalSystem(int i, Array<int> &ind, 
						 Array<double> &locrhs)
{

  double g00,g01,g10,g11;
  g00 = data->GetNodalNeumannCoeff(ind[0], 0);
  g01 = data->GetNodalNeumannCoeff(ind[1], 0);
  g10 = data->GetNodalNeumannCoeff(ind[0], 1);
  g11 = data->GetNodalNeumannCoeff(ind[1], 1);
  
  Array<double *> coord(2);
  coord[0] = new double[2];
  coord[1] = new double[2];
  
  mesh->GetBdrElementVerticesCoord(i, coord);

  double dx, dy, leng;
  dx = coord[0][0] - coord[0][1];
  dy = coord[1][0] - coord[1][1];

  leng = sqrt(pow(dx,2) + pow(dy,2));
  
  locrhs[0] = leng * ( 1.0/6.0 * g01 + 1.0/3.0 * g00 );
  locrhs[1] = leng * ( 1.0/3.0 * g01 + 1.0/6.0 * g00 );
  locrhs[2] = leng * ( 1.0/6.0 * g11 + 1.0/3.0 * g10 );
  locrhs[3] = leng * ( 1.0/3.0 * g11 + 1.0/6.0 * g10 );
  
  delete []coord[0];
  delete []coord[1];

}

//==============================================================================

void BilinearFEMELAST2D::ComputeNeumannLocalSystem(int b, int i, Array<int> &ind, 
					    Array<double> &locrhs)
{
  int nb = mesh->GetNBdrs();

  double g00,g01,g10,g11;
  g00 = data->GetNeumannBdryVal(b, ind[0]);
  g01 = data->GetNeumannBdryVal(b, ind[1]);
  g10 = data->GetNeumannBdryVal(b+nb, ind[0]);
  g11 = data->GetNeumannBdryVal(b+nb, ind[1]);
  
  Array<double *> coord(2);
  coord[0] = new double[2];
  coord[1] = new double[2];
  
  mesh->GetBdrElementVerticesCoord(i, coord);

  double dx, dy, leng;
  dx = coord[0][0] - coord[0][1];
  dy = coord[1][0] - coord[1][1];

  leng = sqrt(pow(dx,2) + pow(dy,2));
  
  locrhs[0] = leng * ( 1.0/6.0 * g01 + 1.0/3.0 * g00 );
  locrhs[1] = leng * ( 1.0/3.0 * g01 + 1.0/6.0 * g00 );
  locrhs[2] = leng * ( 1.0/6.0 * g11 + 1.0/3.0 * g10 );
  locrhs[3] = leng * ( 1.0/3.0 * g11 + 1.0/6.0 * g10 );

 for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void BilinearFEMELAST2D::ComputeRobinLocalSystem(int i, int b, Array<int> &ind, 
					  Array<double *> &locmat, 
					  Array<double> &locrhs)
{
  int nb = mesh->GetNBdrs();

  double g00,g01,g10,g11;
  g00 = -data->GetRobinCoeff(b, ind[0]);
  g01 = -data->GetRobinCoeff(b, ind[1]);

  g10 = -data->GetRobinCoeff(b+nb, ind[0]);
  g11 = -data->GetRobinCoeff(b+nb, ind[1]);

  double n00,n01,n10,n11;
  n00 = data->GetRobinBdryVal(b, ind[0]);
  n01 = data->GetRobinBdryVal(b, ind[1]);

  n10 = data->GetRobinBdryVal(b+nb, ind[0]);
  n11 = data->GetRobinBdryVal(b+nb, ind[1]);

  Array<double *> coord(2);
  coord[0] = new double[2];
  coord[1] = new double[2];
  
  mesh->GetBdrElementVerticesCoord(i, coord);

  double dx, dy, leng;
  dx = coord[0][0] - coord[0][1];
  dy = coord[1][0] - coord[1][1];
  leng = sqrt(pow(dx,2) + pow(dy,2));

  locmat[0][0] = leng * (g01 + 3.0*g00)/12.0;
  locmat[0][1] = leng * (g00 + g01)/12.0;

  locmat[0][2] = 0.0;
  locmat[0][3] = 0.0; 

  locmat[1][0] = leng * (g00 + g01)/12.0;
  locmat[1][1] = leng * (3.0*g01 + g00)/12.0;

  locmat[1][2] = 0.0;
  locmat[1][3] = 0.0;

  locmat[2][0] = 0.0;
  locmat[2][1] = 0.0;


  locmat[2][2] = leng * (g11 + 3.0*g10)/12.0;
  locmat[2][3] = leng * (g10 + g11)/12.0;

  locmat[3][0] = 0.0;
  locmat[3][1] = 0.0;

  locmat[3][2] = leng * (g10 + g11)/12.0;
  locmat[3][3] = leng * (3.0*g11 + g10)/12.0;

  
  locrhs[0] = leng * ( 1.0/6.0 * n01 + 1.0/3.0 * n00 );
  locrhs[1] = leng * ( 1.0/3.0 * n01 + 1.0/6.0 * n00 );

  locrhs[2] = leng * ( 1.0/6.0 * n11 + 1.0/3.0 * n10 );
  locrhs[3] = leng * ( 1.0/3.0 * n11 + 1.0/6.0 * n10 );

  delete []coord[0];
  delete []coord[1];
}

//==============================================================================

void BilinearFEMELAST2D::ProjectAFunction(double (*func)(Array<double> &), 
					Vector &pval)
{
  
  Array<double> coord(2);
  for(int i=0; i<mesh->GetNV(); i++)
    {
      coord[0] = mesh->GetVertex(i, 0);
      coord[1] = mesh->GetVertex(i, 1);      
      pval(i) = func(coord);
    }
  
}


//==============================================================================

void BilinearFEMELAST2D::ComputeConservativeStress(Vector &sol,
						 Array<double *> &stress)
{ 

}
//=============================================================================

void BilinearFEMELAST2D::ComputeBasisStress(Vector &sol,
					  Array<double *> &basisstress)
{ 


}


//==============================================================================

void BilinearFEMELAST2D::SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
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
	      
	      // This assumes row is returned as dofs associated with first
	      // direction followed by second and so on. In 2d triangular 
	      // case it contains row[0] = dofv1, row[1] = dofv2, 
	      // row[2] = dofv1 + NV, row[3] =  dofv2 +NV. See the function in 
	      // linearfemelast2d.cpp for an example.

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
  //mat->Finalize();
  // mat->Print();
}

//==============================================================================
double BilinearFEMELAST2D::ComputeL2Error(Function *exactsol, Vector &approxsol, double time)
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
