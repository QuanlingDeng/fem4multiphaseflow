#include "fem_header.h"

//==============================================================================

QuadraticBSplineFEM2D::QuadraticBSplineFEM2D(Mesh *_mesh, Data *_data) :
FEM(_mesh, _data)
{
    strcpy(Type, "QuadraticBSplineFEM2D");
    NumGlobalDOF = (mesh->GetNx()+2)*(mesh->GetNy()+2);
    NumLocalDOF = 9;  
    rhs.SetSize(NumGlobalDOF);
    mat = new SparseMatrix(NumGlobalDOF, NumGlobalDOF);
}

//============================================================================== 

void QuadraticBSplineFEM2D::GetElementDOF(int i, Array<int> &ind)
{
  //GetElementVertices(i, ind);
  ind.SetSize(9);
  int counter = 0;
  for (int j=0; j<3; j++)
    {
      int row = i/mesh->GetNx()+j;
      ind[counter] = row * (mesh->GetNx()+2) + i % mesh->GetNx();
      counter ++;
      ind[counter] = ind[counter-1] + 1;
      counter++;
      ind[counter] = ind[counter-1] + 1;
      counter++;
    }
}

//============================================================================== 

void QuadraticBSplineFEM2D::GetBdrElementDOF(int i, Array<int> &ind)
{
  ind.SetSize(3);
  int att = mesh->GetBdrAttribute(i);
  int s1 = -((att-3)*(att-1)*(att-2))/6;
  int s2 = (att*(att-2)*(att-3))/2;
  int s3 = -(att*(att-1)*(att-3))/2; 
  int s4 =(att*(att-1)*(att-2))/6;
  for (int j=0; j<3; j++)
    {
      ind[j] = (s1*(i+j) + s2*(mesh->GetNx()-1+(i-mesh->GetNx()+j)*(mesh->GetNx()+2)) + s3*(i-mesh->GetNx()-mesh->GetNy()+j+(mesh->GetNx()+2)*(mesh->GetNy()+1)) + s4*((i-2*mesh->GetNx()-mesh->GetNy()+j)*(mesh->GetNx()+2)));
    }
}

//==============================================================================
void QuadraticBSplineFEM2D::ProjectAFunction(double (*func)(Array<double> &), Vector &pval)
{
  pval.SetSize(NumGlobalDOF);
  Array<double> coordb(2);
  Array<double> coordt(2);
  double *wts = new double[3]; // use 3-pt gq to integrate the L2 product of the "exterior" basis splines
  double *pts = new double[3];
  wts[0] = 8.0/9.0;
  wts[1] = 5.0/9.0;
  wts[2] = wts[1];
  pts[0] = 0.0;
  pts[1] = -sqrt(0.6);
  pts[2] = -pts[1];
  for(int i=0; i<NumGlobalDOF; i++)
    {
      pval(i) = 0;
    }
  Array<Vector *> alpha(2);
  SparseMatrix *proj = new SparseMatrix(mesh->GetNx()+2,mesh->GetNx()+2);
  Vector *prhsb = new Vector(mesh->GetNx()+2); 
  Vector *prhst = new Vector(mesh->GetNx()+2); 
  alpha[0] = new Vector(prhsb->Size());
  alpha[1] = new Vector(prhsb->Size());
  for (int i=0; i<prhsb->Size(); i++)
    {
      prhsb->Elem(i) = 0.0;
      alpha[0]->Elem(i) = 0.0;
      alpha[1]->Elem(i) = 0.0;
    }
  gsl_vector *B;
  gsl_bspline_workspace *bw; 
  B = gsl_vector_alloc(mesh->GetNx()+2);
  bw = gsl_bspline_alloc(3,mesh->GetNx()+1);
  gsl_bspline_knots_uniform(0.0,mesh->GetVertex(mesh->GetNx(),0),bw);
  //horizontal
  coordb[1] = 0.0;
  coordt[1] = mesh->GetVertex(mesh->GetNy()*(mesh->GetNx()+1),1);
  for (int i=0; i<mesh->GetNx(); i++)
    {
      double a = mesh->GetVertex(i,0);
      double b = mesh->GetVertex(i+1,0);
      for (int j=0; j<3; j++)
	{
	  coordb[0] = (b-a)/2.0*pts[j]+(b+a)/2.0;
	  coordt[0] = coordb[0];
	  gsl_bspline_eval(coordb[0],B,bw);
	  for (int k=0; k<3; k++)
	    {
	      prhsb->Elem(i+k) +=  wts[j] * func(coordb) * gsl_vector_get(B,i+k); 
	      prhst->Elem(i+k) +=  wts[j] * func(coordt) * gsl_vector_get(B,i+k); 
	      proj->Elem(i+k,i+k) +=  wts[j] * gsl_vector_get(B,i+k) * gsl_vector_get(B,i+k);
	      proj->Elem(i,i+k) +=   wts[j] * gsl_vector_get(B,i+k) * gsl_vector_get(B,i);
	    }
	  proj->Elem(i+1,i+2) += wts[j] * gsl_vector_get(B,i+1) * gsl_vector_get(B,i+2);
	}
    }
  for (int i=0; i<mesh->GetNx()+2; i++)
    {
      double a = mesh->GetVertex(i,0);
      double b = mesh->GetVertex(i+1,0);
      prhsb->Elem(i) *= (b-a)/2;
      prhst->Elem(i) *= (b-a)/2;
      proj->Elem(i,i) *= (b-a)/2;
    }
  for (int i=0; i<mesh->GetNx(); i++)
    {
      double a = mesh->GetVertex(i,0);
      double b = mesh->GetVertex(i+1,0);
      proj->Elem(i,i+2) *= (b-a)/2;
      proj->Elem(i,i+1) *= (b-a)/2;
      proj->Elem(i+1,i+2) *= (b-a)/2;
      proj->Elem(i+2,i) = proj->Elem(i,i+2);
      proj->Elem(i+1,i) = proj->Elem(i,i+1);
      proj->Elem(i+2,i+1) = proj->Elem(i+1,i+2);
    }
  proj->Elem(0,0) = 1.0;
  proj->Elem(mesh->GetNx()+1,mesh->GetNx()+1) = 1.0;
  for (int i=1; i<3; i++)
    {
      proj->Elem(0,i) = 0;
      proj->Elem(mesh->GetNx()+1,mesh->GetNx()+1-i) = 0;
    }
  proj->Finalize();
  coordb[0] = 0;
  coordb[1] = 0;
  coordt[0] = 0;
  prhsb->Elem(0) = func(coordb);
  prhst->Elem(0) = func(coordt);
  coordb[0] = mesh->GetVertex(mesh->GetNx(),0);
  prhsb->Elem(mesh->GetNx()+1) = func(coordb);
  coordt[0] = mesh->GetVertex(mesh->GetNx(),0);
  prhst->Elem(mesh->GetNx()+1) = func(coordt);
  MatrixInverse *invmat;
  MatrixInverse *prec = new GSSmoother(*proj);
  invmat = new PCGMatrixInverse(*proj, *prec);
  invmat->Mult(*prhsb, *alpha[0]);
  invmat->Mult(*prhst, *alpha[1]);
int ind;
  for (int i; i<mesh->GetNx()+2; i++)
    {
      pval(i) = alpha[0]->Elem(i);
      ind = (mesh->GetNx()+2)*(mesh->GetNy()+1);
      pval(i+ind) = alpha[1]->Elem(i);
    }
  delete invmat;
  delete prec;
  delete proj;
  gsl_vector_free(B);
  gsl_bspline_free(bw);

  //vertical
  prhsb->SetSize(mesh->GetNy()+2);
  prhst->SetSize(mesh->GetNy()+2);
  gsl_vector *By;
  gsl_bspline_workspace *bwy; 
  By = gsl_vector_alloc(mesh->GetNy()+2);
  bwy = gsl_bspline_alloc(3,mesh->GetNy()+1);
  gsl_bspline_knots_uniform(0.0,mesh->GetVertex(mesh->GetNy()*(mesh->GetNx()+1),1),bw);
  SparseMatrix *projv = new SparseMatrix(mesh->GetNy()+2,mesh->GetNy()+2);
  alpha[0]->SetSize(prhsb->Size());
  alpha[1]->SetSize(prhsb->Size());
  for (int i=0; i<alpha[0]->Size(); i++)
    {
      alpha[0]->Elem(i) = 0.0;
      alpha[1]->Elem(i) = 0.0;
      prhsb->Elem(i) = 0.0;
      prhst->Elem(i) = 0.0;
    }
  
  coordb[0] = mesh->GetVertex(mesh->GetNx(),0); //rhs
  coordt[0] = 0.0; //lhs
  for (int i=0; i<mesh->GetNy(); i++)
    {
      double a = mesh->GetVertex(i*(mesh->GetNx()+1),1);
      double b = mesh->GetVertex((i+1)*(mesh->GetNx()+1),1);
      for (int j=0; j<3; j++)
	{
	  coordb[1] = (b-a)/2.0*pts[j]+(b+a)/2.0;
	  coordt[1] = coordb[0];
	  gsl_bspline_eval(coordb[1],By,bwy);
	  for (int k=0; k<3; k++)
	    {
	      prhsb->Elem(i+k) +=  wts[j] * func(coordb) * gsl_vector_get(By,i+k); 
	      prhst->Elem(i+k) +=  wts[j] * func(coordt) * gsl_vector_get(By,i+k); 
	      proj->Elem(i+k,i+k) +=  wts[j] * gsl_vector_get(By,i+k) * gsl_vector_get(By,i+k);
	      proj->Elem(i,i+k) +=   wts[j] * gsl_vector_get(By,i+k) * gsl_vector_get(By,i);
	    }
	  proj->Elem(i+1,i+2) += wts[j] * gsl_vector_get(By,i+1) * gsl_vector_get(By,i+2);
	}
    }
 for (int i=0; i<mesh->GetNy()+2; i++)
    {
      double a = mesh->GetVertex(i*(mesh->GetNx()+1),1);
      double b = mesh->GetVertex((i+1)*(mesh->GetNx()+1),1);
      prhsb->Elem(i) *= (b-a)/2;
      prhst->Elem(i) *= (b-a)/2;
      projv->Elem(i,i) *= (b-a)/2;
    }
  for (int i=0; i<mesh->GetNy(); i++)
    {
      double a = mesh->GetVertex(i*(mesh->GetNx()+1),1);
      double b = mesh->GetVertex((i+1)*(mesh->GetNx()+1),1);
      projv->Elem(i,i+2) *= (b-a)/2;
      projv->Elem(i,i+1) *= (b-a)/2;
      projv->Elem(i+1,i+2) *= (b-a)/2;
      projv->Elem(i+2,i) = projv->Elem(i,i+2);
      projv->Elem(i+1,i) = projv->Elem(i,i+1);
      projv->Elem(i+2,i+1) = projv->Elem(i+1,i+2);
    }
  projv->Elem(0,0) = 1.0;
  projv->Elem(mesh->GetNy()+1,mesh->GetNy()+1) = 1.0;
  for (int i=1; i<3; i++)
    {
      projv->Elem(0,i) = 0;
      projv->Elem(mesh->GetNy()+1,mesh->GetNy()+1-i) = 0;
    }
  projv->Finalize();
  coordb[1] = 0;
  coordt[1] = 0;
  prhsb->Elem(0) = func(coordb);
  prhst->Elem(0) = func(coordt);
  coordb[1] = mesh->GetVertex((mesh->GetNx()+1)*mesh->GetNy(),1);
  prhsb->Elem(mesh->GetNy()+1) = func(coordb);
  coordt[1] = coordb[1];
  prhst->Elem(mesh->GetNy()+1) = func(coordt);

  MatrixInverse *invmaty;
  MatrixInverse *precy = new GSSmoother(*projv);
  invmaty = new PCGMatrixInverse(*projv, *precy);
  invmaty->Mult(*prhsb, *alpha[0]);
  invmaty->Mult(*prhst, *alpha[1]);


  for (int i=0; i<mesh->GetNy()+2; i++)
    {
      ind = i*(mesh->GetNx()+2);
      pval(ind) = alpha[1]->Elem(i);
      pval(ind + mesh->GetNx()+1) = alpha[0]->Elem(i);
    }

  delete invmaty;
  delete precy;
  gsl_vector_free(By);
  gsl_bspline_free(bwy);
  delete alpha[0];
  delete alpha[1];
  delete prhsb;
  delete prhst;
  delete []wts;
  delete []pts;

}

//==============================================================================
void QuadraticBSplineFEM2D::ComputeEllipticLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  double *wts = new double[3]; // use 3-pt gq to integrate the L2 product of the basis splines in both x and y directions
  double *pts = new double[3];
  wts[0] = 8.0/9.0;
  wts[1] = 5.0/9.0;
  wts[2] = wts[1];
  pts[0] = 0.0;
  pts[1] = -sqrt(0.6);
  pts[2] = -pts[1];
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
  gsl_bspline_workspace *bwx;      
  gsl_matrix *Bx;                 
  gsl_bspline_deriv_workspace *wx;
  gsl_bspline_workspace *bwy;      
  gsl_matrix *By;                 
  gsl_bspline_deriv_workspace *wy;
  bwx = gsl_bspline_alloc(3,mesh->GetNx()+1);
  wx = gsl_bspline_deriv_alloc(3);
  Bx = gsl_matrix_alloc(mesh->GetNx()+2,2);
  gsl_bspline_knots_uniform(0.0,mesh->GetVertex(mesh->GetNx(),0),bwx);
  bwy = gsl_bspline_alloc(3,mesh->GetNx()+1);
  wy = gsl_bspline_deriv_alloc(3);
  By = gsl_matrix_alloc(mesh->GetNy()+2,2);
  gsl_bspline_knots_uniform(0.0,mesh->GetVertex(mesh->GetNy()*(mesh->GetNx()+1),1),bwy);

  double ax = coord[0][0];
  double bx = coord[0][1];
  double ay = coord[1][0];
  double by = coord[1][2];
 
  double k1 = locdat[0];
  double k2 = locdat[1];
  double k3 = locdat[2];
  double k4 = locdat[3];



  for (int i=0; i<3; i++)
    {
      for (int j=0; j<3; j++)
	{
	  for (int k=0; k<9; k++)
	    {
	      for (int l=k; l<9; l++)
		{
		 
		}
	  /*  locmat[0][0] = 1.0/(24.0*lenx*leny)*( (-k4)*(-3.0*lenx*lenx-leny*leny) + 3.0*k1*(lenx*lenx+leny*leny) + k3*(lenx*lenx+leny*leny) + k2*(lenx*lenx+3.0*leny*leny) );
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
	  locmat[3][3] = 1.0/(24.0*lenx*leny)*(  -k1*(-3.0*lenx*lenx-leny*leny) + k2*(lenx*lenx+leny*leny) + 3.0*k4*(lenx*lenx+leny*leny) + k3*(lenx*lenx+3.0*leny*leny)  );*/
	    }
	}
    }
  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
  delete []wts;
  delete []pts;
  gsl_matrix_free(Bx);
  gsl_bspline_free(bwx);
  gsl_bspline_deriv_free(wx);
  gsl_matrix_free(By);
  gsl_bspline_free(bwy);
  gsl_bspline_deriv_free(wy);
}
/*
//==============================================================================

void QuadraticBSplineFEM2D::ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
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

void QuadraticBSplineFEM2D::ComputeForceLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
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
void QuadraticBSplineFEM2D::SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
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
void QuadraticBSplineFEM2D::ComputeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
}

//==============================================================================

void QuadraticBSplineFEM2D::ComputeConservativeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
}

//==============================================================================

void QuadraticBSplineFEM2D::ComputeNeumannLocalSystem(int i, Array<int> &ind, 
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
  locrhs[0] = h * ( 2.0*locdat[0] + locdat[1] ) / 6.0;
  locrhs[1] = h * ( locdat[0] + 2.0*locdat[1] ) / 6.0;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void QuadraticBSplineFEM2D::ComputeNeumannLocalSystem(int b, int i, Array<int> &ind, 
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

void QuadraticBSplineFEM2D::ComputeRobinLocalSystem(int i, Array<int> &ind, 
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

void QuadraticBSplineFEM2D::ComputeRobinLocalSystem(int i, int b, Array<int> &ind, 
					  Array<double *> &locmat, 
					  Array<double> &locrhs)
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
double QuadraticBSplineFEM2D::ComputeL2Error(double (*solfunc)(Array<double> &), Vector &approxsol)
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

double QuadraticBSplineFEM2D::ComputeH1Error(double (*deronefunc)(Array<double> &), double (*dertwofunc)(Array<double> &), Vector &approxsol)
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

double QuadraticBSplineFEM2D::ComputeL2Error(Function *exactsol, Vector &approxsol)
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
      double coord[2];
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
double QuadraticBSplineFEM2D::ComputeH1Error(Array<Function *> &exactderivative, Vector &approxsol)
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

void QuadraticBSplineFEM2D::ComputeStableLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{

}

//==============================================================================

void QuadraticBSplineFEM2D::ComputeStableRHSLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
{

}
*/
