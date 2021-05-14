#include "fem_header.h"
#include "../general/general_header.h"

LinearFEM2D::LinearFEM2D(Mesh *_mesh, Data *_data) :
FEM(_mesh, _data)
{
    strcpy(Type, "LinearFEM2D");
    NumGlobalDOF = mesh->GetNV();
    NumLocalDOF = 3;
    rhs.SetSize(NumGlobalDOF);
    mat = new SparseMatrix(NumGlobalDOF, 10);
}

//============================================================================== 

void LinearFEM2D::GetElementDOF(int i, Array<int> &ind)
{
  GetElementVertices(i, ind);
}

//============================================================================== 

void LinearFEM2D::GetBdrElementDOF(int i, Array<int> &ind)
{
  GetBdrElementVertices(i, ind);
}

//==============================================================================

void LinearFEM2D::ComputeEllipticLocalSystem(int i, Array<int> &ind, 
					     Array<double *> &locmat)
{
  
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalEllipticCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[1][2] - coord[1][0]) 
                 - (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);
  double t = ( locdat[0] + locdat[1] + locdat[2] ) / ( 6.0*detJ );

  locmat[0][0] = t * ( (coord[1][1] - coord[1][2]) * (coord[1][1] 
                 - coord[1][2]) + (coord[0][2] - coord[0][1]) * 
                 (coord[0][2] - coord[0][1]) );
  locmat[0][1] = t * ( (coord[1][1] - coord[1][2]) * (coord[1][2] 
                 - coord[1][0]) + (coord[0][2] - coord[0][1]) * 
                 (coord[0][0] - coord[0][2]) );
  locmat[0][2] = t * ( (coord[1][1] - coord[1][2]) * (coord[1][0] 
                 - coord[1][1]) + (coord[0][2] - coord[0][1]) * 
                 (coord[0][1] - coord[0][0]) );

  locmat[1][0] = locmat[0][1];
  locmat[1][1] = t * ( (coord[1][0] - coord[1][2]) * (coord[1][0] 
                 - coord[1][2]) + (coord[0][2] - coord[0][0]) * 
                 (coord[0][2] - coord[0][0]) );
  locmat[1][2] = t * ( (coord[1][2] - coord[1][0]) * (coord[1][0] 
                 - coord[1][1]) + (coord[0][0] - coord[0][2]) * 
                 (coord[0][1] - coord[0][0]) );

  locmat[2][0] = locmat[0][2];
  locmat[2][1] = locmat[1][2];
  locmat[2][2] =  t * ( (coord[1][0] - coord[1][1]) * (coord[1][0] 
                  - coord[1][1]) + (coord[0][1] - coord[0][0]) * 
                  (coord[0][1] - coord[0][0]) );
    
  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
  
  /*  
  for(int k=0; k<locmat.Size(); k++)
    for(int j=0; j<locmat.Size(); j++)
      locmat[k][j] = 0.0;
	  
  int np = ConvertN2PN(5);
  Array<double *> tpw(np);
  for(int j=0; j<np; j++) tpw[j] = new double[3];
  GetStandTriQuadPW(5, tpw);

  Array<double *> c(mesh->GetDim()); // c for coordinate
  for (int j=0; j<c.Size(); j++)
    c[j] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, c);

  double t1 = c[0][1] - c[0][0];
  double t2 = c[0][2] - c[0][0];
  double t3 = c[1][1] - c[1][0];
  double t4 = c[1][2] - c[1][0];
  double detJ = t1 * t4 - t2 * t3;
    
  // Jacobian Inverse Transpose * detJ
  Array<double *> jit(2);
  jit[0] = new double[2];
  jit[1] = new double[2];
  jit[0][0] = c[1][2] - c[1][0];
  jit[0][1] = c[1][0] - c[1][1];
  jit[1][0] = c[0][0] - c[0][2];
  jit[1][1] = c[0][1] - c[0][0];

  double gpx, gpy, xx, yy, dotprod, vm[2], vn[2];

  BasisFunctions *phi = new BasisFunctions(1);
  Function *kk = data->GetEllipticFunction();
  
  for(int m=0; m<locmat.Size(); m++)
    {
      for(int n=m; n<locmat.Size(); n++)
	{
	  for(int k=0; k<np; k++)
	    {
	      gpx = tpw[k][0];
	      gpy = tpw[k][1];
	      xx = c[0][0] + t1*gpx + t2*gpy;
	      yy = c[1][0] + t3*gpx + t4*gpy;

	      vn[0] = jit[0][0] * phi->GradBF2D(n, 0, gpx, gpy) + jit[0][1] * phi->GradBF2D(n, 1, gpx, gpy);
	      vn[1] = jit[1][0] * phi->GradBF2D(n, 0, gpx, gpy) + jit[1][1] * phi->GradBF2D(n, 1, gpx, gpy);

	      vm[0] = jit[0][0] * phi->GradBF2D(m, 0, gpx, gpy) + jit[0][1] * phi->GradBF2D(m, 1, gpx, gpy);
	      vm[1] = jit[1][0] * phi->GradBF2D(m, 0, gpx, gpy) + jit[1][1] * phi->GradBF2D(m, 1, gpx, gpy);

	      dotprod = vn[0] * vm[0] + vn[1] * vm[1];

	      double coord[2];
	      coord[0] = xx;
	      coord[1] = yy;
	      locmat[m][n] += tpw[k][2] * kk->Eval(coord) * dotprod / detJ;
	    } 	  	  
	}
    }
  
  for(int m=0; m<locmat.Size(); m++)
    {
      for(int n=m+1; n<locmat.Size(); n++)
	{
	  locmat[n][m] = locmat[m][n];
	}
    }

  for (int k=0; k<jit.Size(); k++) { delete []jit[k]; }
  for (int k=0; k<c.Size(); k++) { delete []c[k]; }
  */
}

//==============================================================================

void LinearFEM2D::ComputeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
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
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  
  //build local matrix 
  locmat[0][0] = ( coord[1][1] - coord[1][2] ) * ( 2.0*locdat[0][0] + locdat[0][1] + locdat[0][2] ) / 24.0 + ( coord[0][2] - coord[0][1] ) * ( 2.0*locdat[1][0] + locdat[1][1] + locdat[1][2] ) / 24.0;
  locmat[0][1] = ( coord[1][2] - coord[1][0] ) * ( 2.0*locdat[0][0] + locdat[0][1] + locdat[0][2] ) / 24.0 + ( coord[0][0] - coord[0][2] ) * ( 2.0*locdat[1][0] + locdat[1][1] + locdat[1][2] ) / 24.0;
  locmat[0][2] = ( coord[1][0] - coord[1][1] ) * ( 2.0*locdat[0][0] + locdat[0][1] + locdat[0][2] ) / 24.0 + ( coord[0][1] - coord[0][0] ) * ( 2.0*locdat[1][0] + locdat[1][1] + locdat[1][2] ) / 24.0;

  locmat[1][0] = ( coord[1][1] - coord[1][2] ) * ( locdat[0][0] + 2.0*locdat[0][1] + locdat[0][2] ) / 24.0 + ( coord[0][2] - coord[0][1] ) * ( locdat[1][0] + 2.0*locdat[1][1] + locdat[1][2] ) / 24.0;
  locmat[1][1] = ( coord[1][2] - coord[1][0] ) * ( locdat[0][0] + 2.0*locdat[0][1] + locdat[0][2] ) / 24.0 + ( coord[0][0] - coord[0][2] ) * ( locdat[1][0] + 2.0*locdat[1][1] + locdat[1][2] ) / 24.0;
  locmat[1][2] = ( coord[1][0] - coord[1][1] ) * ( locdat[0][0] + 2.0*locdat[0][1] + locdat[0][2] ) / 24.0 + ( coord[0][1] - coord[0][0] ) * ( locdat[1][0] + 2.0*locdat[1][1] + locdat[1][2] ) / 24.0;

  locmat[2][0] = ( coord[1][1] - coord[1][2] ) * ( locdat[0][0] + locdat[0][1] + 2.0*locdat[0][2] ) / 24.0 + ( coord[0][2] - coord[0][1] ) * ( locdat[1][0] + locdat[1][1] + 2.0*locdat[1][2] ) / 24.0;
  locmat[2][1] = ( coord[1][2] - coord[1][0] ) * ( locdat[0][0] + locdat[0][1] + 2.0*locdat[0][2] ) / 24.0 + ( coord[0][0] - coord[0][2] ) * ( locdat[1][0] + locdat[1][1] + 2.0*locdat[1][2] ) / 24.0;
  locmat[2][2] = ( coord[1][0] - coord[1][1] ) * ( locdat[0][0] + locdat[0][1] + 2.0*locdat[0][2] ) / 24.0 + ( coord[0][1] - coord[0][0] ) * ( locdat[1][0] + locdat[1][1] + 2.0*locdat[1][2] ) / 24.0;


  for (int k=0; k<locdat.Size(); k++) { delete []locdat[k]; }
  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; } 
}

//==============================================================================

void LinearFEM2D::ComputeConservativeAdvectionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  Array<double *> locdat(mesh->GetDim());
  for (int k=0; k<locdat.Size(); k++)
    locdat[k] = new double[GetElement(i)->GetNVertices()];

  for (int j=0; j<ind.Size(); j++)
    {
      locdat[0][j] = data->GetNodalConservativeAdvectionCoeff(ind[j], 0);
      locdat[1][j] = data->GetNodalConservativeAdvectionCoeff(ind[j], 1);
    }
  
  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);

  // A_{ij} = - A_{ji}, where the second A is from ComputeAdvectionLocalSystem
  locmat[0][0] = -( coord[1][1] - coord[1][2] ) * ( 2.0*locdat[0][0] + locdat[0][1] + locdat[0][2] ) / 24.0 - ( coord[0][2] - coord[0][1] ) * ( 2.0*locdat[1][0] + locdat[1][1] + locdat[1][2] ) / 24.0;
  locmat[0][1] = -( coord[1][1] - coord[1][2] ) * ( locdat[0][0] + 2.0*locdat[0][1] + locdat[0][2] ) / 24.0 - ( coord[0][2] - coord[0][1] ) * ( locdat[1][0] + 2.0*locdat[1][1] + locdat[1][2] ) / 24.0;
  locmat[0][2] = -( coord[1][1] - coord[1][2] ) * ( locdat[0][0] + locdat[0][1] + 2.0*locdat[0][2] ) / 24.0 - ( coord[0][2] - coord[0][1] ) * ( locdat[1][0] + locdat[1][1] + 2.0*locdat[1][2] ) / 24.0;

  locmat[1][0] = -( coord[1][2] - coord[1][0] ) * ( 2.0*locdat[0][0] + locdat[0][1] + locdat[0][2] ) / 24.0 - ( coord[0][0] - coord[0][2] ) * ( 2.0*locdat[1][0] + locdat[1][1] + locdat[1][2] ) / 24.0;
  locmat[1][1] = -( coord[1][2] - coord[1][0] ) * ( locdat[0][0] + 2.0*locdat[0][1] + locdat[0][2] ) / 24.0 - ( coord[0][0] - coord[0][2] ) * ( locdat[1][0] + 2.0*locdat[1][1] + locdat[1][2] ) / 24.0;
  locmat[1][2] = -( coord[1][2] - coord[1][0] ) * ( locdat[0][0] + locdat[0][1] + 2.0*locdat[0][2] ) / 24.0 - ( coord[0][0] - coord[0][2] ) * ( locdat[1][0] + locdat[1][1] + 2.0*locdat[1][2] ) / 24.0;

  locmat[2][0] = -( coord[1][0] - coord[1][1] ) * ( 2.0*locdat[0][0] + locdat[0][1] + locdat[0][2] ) / 24.0 - ( coord[0][1] - coord[0][0] ) * ( 2.0*locdat[1][0] + locdat[1][1] + locdat[1][2] ) / 24.0;
  locmat[2][1] = -( coord[1][0] - coord[1][1] ) * ( locdat[0][0] + 2.0*locdat[0][1] + locdat[0][2] ) / 24.0 - ( coord[0][1] - coord[0][0] ) * ( locdat[1][0] + 2.0*locdat[1][1] + locdat[1][2] ) / 24.0;
  locmat[2][2] = -( coord[1][0] - coord[1][1] ) * ( locdat[0][0] + locdat[0][1] + 2.0*locdat[0][2] ) / 24.0 - ( coord[0][1] - coord[0][0] ) * ( locdat[1][0] + locdat[1][1] + 2.0*locdat[1][2] ) / 24.0;

  for (int k=0; k<locdat.Size(); k++) { delete []locdat[k]; }
  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; } 
}

//==============================================================================

void LinearFEM2D::ComputeReactionLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalReactionCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[1][2] - coord[1][0]) 
                - (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);
  
  locmat[0][0] = detJ * ( 3.0*locdat[0] + locdat[1] + locdat[2] ) / 60.0;
  locmat[0][1] = detJ * ( 2.0*locdat[0] + 2.0*locdat[1] + locdat[2] ) / 120.0;
  locmat[0][2] = detJ * ( 2.0*locdat[0] + locdat[1] + 2.0*locdat[2] ) / 120.0;
  locmat[1][0] = locmat[0][1];
  locmat[1][1] = detJ * ( locdat[0] + 3.0*locdat[1] + locdat[2] ) / 60.0;
  locmat[1][2] = detJ * ( locdat[0] + 2.0*locdat[1] + 2.0*locdat[2] ) / 120.0;
  locmat[2][0] = locmat[0][2];
  locmat[2][1] = locmat[0][1];
  locmat[2][2] = detJ * ( locdat[0] + locdat[1] + 3.0*locdat[2] ) / 60.0;

  for (int k=0; k<coord.Size(); k++)
    delete []coord[k];
}

//==============================================================================

void LinearFEM2D::ComputeForceLocalSystem(int i, Array<int> &ind, 
					  Array<double> &locrhs)
{
  
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalForceCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[1][2] - coord[1][0]) 
                - (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);

  locrhs[0] = detJ * ( 2.0*locdat[0] + locdat[1] + locdat[2] ) / 24.0;
  locrhs[1] = detJ * ( locdat[0] + 2.0*locdat[1] + locdat[2] ) / 24.0;
  locrhs[2] = detJ * ( locdat[0] + locdat[1] + 2.0*locdat[2] ) / 24.0;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
  
  /*
  for(int j=0; j<locrhs.Size(); j++) locrhs[j] = 0.0;

  int np = ConvertN2PN(5);
  Array<double *> tpw(np);
  for(int j=0; j<np; j++) tpw[j] = new double[3];
  GetStandTriQuadPW(5, tpw);

  Array<double *> c(mesh->GetDim()); // c for coordinate
  for (int j=0; j<c.Size(); j++)
    c[j] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, c);

  double t1 = c[0][1] - c[0][0];
  double t2 = c[0][2] - c[0][0];
  double t3 = c[1][1] - c[1][0];
  double t4 = c[1][2] - c[1][0];
  double detJ = t1 * t4 - t2 * t3;
  
  double gpx, gpy, xx, yy;

  BasisFunctions *phi = new BasisFunctions(1);
  Function *ff = data->GetForceFunction();
  for(int m=0; m<locrhs.Size(); m++)
    {
      for(int k=0; k<np; k++)
	{
	  gpx = tpw[k][0];
	  gpy = tpw[k][1];
	  xx = c[0][0] + t1*gpx + t2*gpy;
	  yy = c[1][0] + t3*gpx + t4*gpy;

	  double coord[2];
	  coord[0] = xx;
	  coord[1] = yy;
	  locrhs[m] += tpw[k][2] * ff->Eval(coord) * phi->BF2D(m, gpx, gpy) * detJ;
	}
    }
  
  for (int k=0; k<c.Size(); k++) { delete []c[k]; }
  */
}

//==============================================================================

void LinearFEM2D::ProjectAFunction(double (*func)(Array<double> &), Vector &pval)
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

void LinearFEM2D::ComputeNeumannLocalSystem(int i, Array<int> &ind, 
					    Array<double> &locrhs)
{
  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalNeumannCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetBdrElement(i)->GetNVertices()];
  
  mesh->GetBdrElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[0][1] - coord[0][0]) 
                + (coord[1][1] - coord[1][0]) * (coord[1][1] - coord[1][0]);

  detJ = sqrt(detJ);

  locrhs[0] = -detJ * ( 2.0*locdat[0] + locdat[1] ) / 6.0;
  locrhs[1] = -detJ * ( locdat[0] + 2.0*locdat[1] ) / 6.0;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//==============================================================================

void LinearFEM2D::ComputeNeumannLocalSystem(int b, int i, Array<int> &ind, 
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
  locrhs[0] = -h * ( 2.0*locdat[0] + locdat[1] ) / 6.0;
  locrhs[1] = -h * ( locdat[0] + 2.0*locdat[1] ) / 6.0;

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//============================================================================

void LinearFEM2D::ComputeRobinLocalSystem(int i, Array<int> &ind, 
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

void LinearFEM2D::ComputeRobinLocalSystem(int i, int b, Array<int> &ind, 
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

void LinearFEM2D::ComputeStableLocalSystem(int i, Array<int> &ind, Array<double *> &locmat)
{
  double dd = data->GetElementalStableCoeff(i);

  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalAdvectionCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[1][2] - coord[1][0]) - (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);

  Vector gpp(8);
  Vector gp(8);
  Vector gw(4);

  Array<double *> tpw(4);
  for(int i=0; i<4; i++) tpw[i] = new double[3];
  GetStandTriQuadPW(3, tpw);

  double tt1 = coord[0][1] - coord[0][0];
  double tt2 = coord[0][2] - coord[0][0];
  double tt3 = coord[1][1] - coord[1][0];
  double tt4 = coord[1][2] - coord[1][0];

  for(int k=0; k<4; k++)
    {
      gpp(2*k) = tpw[k][0];
      gpp(2*k+1) = tpw[k][1];
      gp(2*k) = coord[0][0] + tt1*gpp(2*k) + tt2*gpp(2*k+1);
      gp(2*k+1) = coord[1][0] + tt3*gpp(2*k) + tt4*gpp(2*k+1);
      gw(k) = 25.0 * detJ / 96.0;	  
    }
  gw(0) = -27.0 * detJ / 96.0;
  for(int k=0; k<4; k++) delete []tpw[k];
  
  double adx = 0.0;
  double ady = 0.0;
  adx = ( locdat[0]*(coord[1][1] - coord[1][2]) + locdat[1]*(coord[1][2] - coord[1][0]) + locdat[2]*(coord[1][0] - coord[1][1]) ) / detJ;
  ady = ( locdat[0]*(coord[0][2] - coord[0][1]) + locdat[1]*(coord[0][0] - coord[0][2]) + locdat[2]*(coord[0][1] - coord[0][0]) ) / detJ;

  double t1 = coord[1][1] - coord[1][2];
  double t2 = coord[0][2] - coord[0][1];
  double sum = 0.0;
  for(int j=0; j<4; j++)
    {
      sum += gw(j) * dd * ( adx * t1 + ady * t2 ) * ( adx * t1 + ady * t2 ) / ( detJ * detJ * ( adx*adx + ady*ady ) );
    }
  locmat[0][0] = sum;

  t1 = coord[1][2] - coord[1][0];
  t2 = coord[0][0] - coord[0][2];
  sum = 0.0;
  for(int j=0; j<4; j++)
    {
      sum += gw(j) * dd * ( adx * t1 + ady * t2 ) * ( adx * t1 + ady * t2 ) / ( detJ * detJ * ( adx*adx + ady*ady ) );
    }
  locmat[1][1] = sum;
  
  t1 = coord[1][0] - coord[1][1];
  t2 = coord[0][1] - coord[0][0];
  sum = 0.0;
  for(int j=0; j<4; j++)
    {
      sum += gw(j) * dd * ( adx * t1 + ady * t2 ) * ( adx * t1 + ady * t2 ) / ( detJ * detJ * ( adx*adx + ady*ady ) );
    }
  locmat[2][2] = sum;

  t1 = coord[1][1] - coord[1][2];
  t2 = coord[0][2] - coord[0][1];
  double t3 = coord[1][2] - coord[1][0];
  double t4 = coord[0][0] - coord[0][2];
  sum = 0.0;
  for(int j=0; j<4; j++)
    {
      sum += gw(j) * dd * ( adx * t1 + ady * t2 ) * ( adx * t3 + ady * t4 ) / ( detJ * detJ * ( adx*adx + ady*ady ) );
    }
  locmat[0][1] = sum;
  locmat[1][0] = sum;

  t1 = coord[1][1] - coord[1][2];
  t2 = coord[0][2] - coord[0][1];
  t3 = coord[1][0] - coord[1][1];
  t4 = coord[0][1] - coord[0][0];
  sum = 0.0;
  for(int j=0; j<4; j++)
    {
      sum += gw(j) * dd * ( adx * t1 + ady * t2 ) * ( adx * t3 + ady * t4 ) / ( detJ * detJ * ( adx*adx + ady*ady ) );
    }
  locmat[0][2] = sum;
  locmat[2][0] = sum;
  
  t1 = coord[1][2] - coord[1][0];
  t2 = coord[0][0] - coord[0][2];
  t3 = coord[1][0] - coord[1][1];
  t4 = coord[0][1] - coord[0][0];
  sum = 0.0;
  for(int j=0; j<4; j++)
    {
      sum += gw(j) * dd * ( adx * t1 + ady * t2 ) * ( adx * t3 + ady * t4 ) / ( detJ * detJ * ( adx*adx + ady*ady ) );
    }
  locmat[1][2] = sum;
  locmat[2][1] = sum;
}

//==============================================================================

void LinearFEM2D::ComputeStableRHSLocalSystem(int i, Array<int> &ind, Array<double> &locrhs)
{
  double dd = data->GetElementalStableCoeff(i);

  Array<double> locdat(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locdat[j] = data->GetNodalForceCoeff(ind[j]);

  Array<double> locd(ind.Size());
  for (int j=0; j<ind.Size(); j++)
    locd[j] = data->GetNodalAdvectionCoeff(ind[j]);

  Array<double *> coord(mesh->GetDim());
  for (int k=0; k<coord.Size(); k++)
    coord[k] = new double[GetElement(i)->GetNVertices()];
  
  mesh->GetElementVerticesCoord(i, coord);
  double detJ = (coord[0][1] - coord[0][0]) * (coord[1][2] - coord[1][0]) - (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);

  locrhs[0] = 0.0;
  locrhs[1] = 0.0;
  locrhs[2] = 0.0;

  // ===== stabilization ==============================
  Vector gpp(8);
  Vector gp(8);
  Vector gw(4);

  Array<double *> tpw(4);
  for(int i=0; i<4; i++) tpw[i] = new double[3];
  GetStandTriQuadPW(3, tpw);

  double tt1 = coord[0][1] - coord[0][0];
  double tt2 = coord[0][2] - coord[0][0];
  double tt3 = coord[1][1] - coord[1][0];
  double tt4 = coord[1][2] - coord[1][0];

  for(int k=0; k<4; k++)
    {
      gpp(2*k) = tpw[k][0];
      gpp(2*k+1) = tpw[k][1];
      gp(2*k) = coord[0][0] + tt1*gpp(2*k) + tt2*gpp(2*k+1);
      gp(2*k+1) = coord[1][0] + tt3*gpp(2*k) + tt4*gpp(2*k+1);
      gw(k) = 25.0 * detJ / 96.0;	  
    }
  gw(0) = -27.0 * detJ / 96.0;
  for(int k=0; k<4; k++) delete []tpw[k];
    
  double adx = 0.0;
  double ady = 0.0;
  adx = ( locd[0]*(coord[1][1] - coord[1][2]) + locd[1]*(coord[1][2] - coord[1][0]) + locd[2]*(coord[1][0] - coord[1][1]) ) / detJ;
  ady = ( locd[0]*(coord[0][2] - coord[0][1]) + locd[1]*(coord[0][0] - coord[0][2]) + locd[2]*(coord[0][1] - coord[0][0]) ) / detJ;

  double add = 0.0;
  double t1 = coord[1][1] - coord[1][2];
  double t2 = coord[0][2] - coord[0][1];
  for(int j=0; j<4; j++)
    {
      add = locdat[0] * ( 1.0-gpp(2*j) - gpp(2*j+1) ) + locdat[1] * gpp(2*j) + locdat[2] * gpp(2*j+1);
      locrhs[0] += gw(j) * ( adx * t1 + ady * t2 ) * add * dd / ( ( adx*adx + ady*ady )  * detJ );
    }
  
  t1 = coord[1][2] - coord[1][0];
  t2 = coord[0][0] - coord[0][2];
  for(int j=0; j<4; j++)
    {
      add = locdat[0] * ( 1.0-gpp(2*j) - gpp(2*j+1) ) + locdat[1] * gpp(2*j) + locdat[2] * gpp(2*j+1);
      locrhs[1] += gw(j) * ( adx * t1 + ady * t2 ) * add * dd / ( ( adx*adx + ady*ady )  * detJ );
    }
  
  t1 = coord[1][0] - coord[1][1];
  t2 = coord[0][1] - coord[0][0];
  for(int j=0; j<4; j++)
    {
      add = locdat[0] * ( 1.0-gpp(2*j) - gpp(2*j+1) ) + locdat[1] * gpp(2*j) + locdat[2] * gpp(2*j+1);
      locrhs[2] += gw(j) * ( adx * t1 + ady * t2 ) * add * dd / ( ( adx*adx + ady*ady )  * detJ );
    }

  for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
}

//============================================================================

void LinearFEM2D::SetDirichletBoundary(Array<bool *> &bdry_is_dirichlet, 
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
	      // row[2] = dofv1 + NV, row[3] =  dofv2+NV. See the function in 
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
  // mat->Finalize();
  // mat->Print();
}

//==============================================================================

double LinearFEM2D::ComputeL2Error(double (*solfunc)(Array<double> &), Vector &approxsol)
{ 
  double uh, tr = 0.0;
  double gp[8] = { 1.0/3.0, 1.0/3.0, 0.2, 0.6, 0.2, 0.2, 0.6, 0.2 };
  double gw[4] = { -27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0 };
  Array<double> x(3);
  Array<double> y(3);
   
  double err = 0.0;
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      mesh->GetElementVertices(i, ind);
       
      for(int k=0; k<3; k++)
	{
	  x[k] = mesh->GetVertex(ind[k], 0);
	  y[k] = mesh->GetVertex(ind[k], 1);	  
	}
      
      double t1 = x[1] - x[0];
      double t2 = x[2] - x[0];
      double t3 = y[1] - y[0];
      double t4 = y[2] - y[0];
      double detJ = t1 * t4 - t2 * t3;

      Array<double> gpcoord(2);
      for(int k=0; k<4; k++)
	{
	  double xx = gp[2*k];
	  double yy = gp[2*k+1];
	  gpcoord[0] = x[0] + t1*xx + t2*yy;
	  gpcoord[1] = y[0] + t3*xx + t4*yy;

	  uh = approxsol(ind[0]) * (1 - xx - yy) + approxsol(ind[1]) * xx + approxsol(ind[2]) * yy;
	  tr = solfunc(gpcoord);
	  err += gw[k] * (uh - tr)*(uh - tr) * detJ;	 
	}
    }  

  return sqrt(err);
  //cout<<"The L2 norm error is: "<< sqrt(err)<<endl;
}

//==============================================================================

double LinearFEM2D::ComputeL2Error(Function *exactsol, Vector &femsol, double time)
{ 
  double uh, tr = 0.0;
  double gp[8] = { 1.0/3.0, 1.0/3.0, 0.2, 0.6, 0.2, 0.2, 0.6, 0.2 };
  double gw[4] = { -27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0 };
  Array<double> x(3);
  Array<double> y(3);
   
  double err = 0.0;
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      mesh->GetElementVertices(i, ind);
       
      for(int k=0; k<3; k++)
	{
	  x[k] = mesh->GetVertex(ind[k], 0);
	  y[k] = mesh->GetVertex(ind[k], 1);	  
	}
      
      double t1 = x[1] - x[0];
      double t2 = x[2] - x[0];
      double t3 = y[1] - y[0];
      double t4 = y[2] - y[0];
      double detJ = t1 * t4 - t2 * t3;

      double gpcoord[3];
      gpcoord[2] = time;
      for(int k=0; k<4; k++)
	{
	  double xx = gp[2*k];
	  double yy = gp[2*k+1];
	  gpcoord[0] = x[0] + t1*xx + t2*yy;
	  gpcoord[1] = y[0] + t3*xx + t4*yy;

	  uh = femsol(ind[0]) * (1 - xx - yy) + femsol(ind[1]) * xx + femsol(ind[2]) * yy;
	  tr = exactsol->Eval(gpcoord);
	  err += gw[k] * (uh - tr)*(uh - tr) * detJ;	 
	}
    }  
  return sqrt(err);
}

//==============================================================================

void LinearFEM2D::ComputePPL2Error(Array<double> &errs, Function *exactsol, Vector &femsol, Array<double *> &ppsol, double time)
{
  errs[0] = 0.0;
  errs[1] = 0.0;
  errs[2] = 0.0;

  double gp[8] = { 1.0/3.0, 1.0/3.0, 0.2, 0.6, 0.2, 0.2, 0.6, 0.2 };
  double gw[4] = { -27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0 };
  Array<double> x(3);
  Array<double> y(3);
   
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      mesh->GetElementVertices(i, ind);
       
      for(int k=0; k<3; k++)
	{
	  x[k] = mesh->GetVertex(ind[k], 0);
	  y[k] = mesh->GetVertex(ind[k], 1);	  
	}
      
      double t1 = x[1] - x[0];
      double t2 = x[2] - x[0];
      double t3 = y[1] - y[0];
      double t4 = y[2] - y[0];
      double detJ = t1 * t4 - t2 * t3; 

      double gpcoord[3];
      gpcoord[2] = time;
      for(int k=0; k<4; k++)
	{
	  double xx = gp[2*k];
	  double yy = gp[2*k+1];
	  gpcoord[0] = x[0] + t1*xx + t2*yy;
	  gpcoord[1] = y[0] + t3*xx + t4*yy;

	  double uu = exactsol->Eval(gpcoord);
	  double femuh = femsol(ind[0]) * (1 - xx - yy) + femsol(ind[1]) * xx + femsol(ind[2]) * yy;
	  double ppuh = ppsol[0][i] * (1 - xx - yy) + ppsol[1][i] * xx + ppsol[2][i] * yy;

	  double tem = uu - femuh;
	  errs[0] += gw[k] * tem * tem * detJ;	 
	  tem = uu - ppuh;
	  errs[1] += gw[k] * tem * tem * detJ;	 
	  tem = ppuh - femuh;
	  errs[2] += gw[k] * tem * tem * detJ;	 
	}
    }

  errs[0] = sqrt(errs[0]);
  errs[1] = sqrt(errs[1]);
  errs[2] = sqrt(errs[2]);
}

//==============================================================================

double LinearFEM2D::ComputeH1Error(double (*deronefunc)(Array<double> &), double (*dertwofunc)(Array<double> &), Vector &approxsol)
{ 
  double uh, tr = 0.0;
  double gp[8] = { 1.0/3.0, 1.0/3.0, 0.2, 0.6, 0.2, 0.2, 0.6, 0.2 };
  double gw[4] = { -27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0 };
  Array<double> x(3);
  Array<double> y(3);
   
  double herr = 0.0;
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      mesh->GetElementVertices(i, ind);
       
      for(int k=0; k<3; k++)
	{
	  x[k] = mesh->GetVertex(ind[k], 0);
	  y[k] = mesh->GetVertex(ind[k], 1);	  
	}
      
      double t1 = x[1] - x[0];
      double t2 = x[2] - x[0];
      double t3 = y[1] - y[0];
      double t4 = y[2] - y[0];

      double s1 = y[1] - y[2];
      double s2 = y[2] - y[0];
      double s3 = y[0] - y[1];
      double s4 = x[2] - x[1];
      double s5 = x[0] - x[2];
      double s6 = x[1] - x[0];
      double detJ = t1 * t4 - t2 * t3;

      Array<double> gpcoord(2);
      for(int k=0; k<4; k++)
	{
	  double xx = gp[2*k];
	  double yy = gp[2*k+1];
	  gpcoord[0] = x[0] + t1*xx + t2*yy;
	  gpcoord[1] = y[0] + t3*xx + t4*yy;
 
	  uh = ( approxsol(ind[0])*s1 + approxsol(ind[1])*s2 + approxsol(ind[2])*s3 ) / detJ;
	  uh -= deronefunc(gpcoord);
	  uh = uh*uh;
	  tr = ( approxsol(ind[0])*s4 + approxsol(ind[1])*s5 + approxsol(ind[2])*s6 ) / detJ;
	  tr -= dertwofunc(gpcoord);
	  tr = tr*tr;
	  herr += ( uh + tr ) * gw[k] * detJ;	
	}
    }  
  return sqrt(herr);
}

//==============================================================================

double LinearFEM2D::ComputeH1Error(Array<Function *> &exactderivative, Vector &approxsol)
{ 
  double uh, tr = 0.0;
  double gp[8] = { 1.0/3.0, 1.0/3.0, 0.2, 0.6, 0.2, 0.2, 0.6, 0.2 };
  double gw[4] = { -27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0 };
  Array<double> x(3);
  Array<double> y(3);
   
  double herr = 0.0;
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      mesh->GetElementVertices(i, ind);
       
      for(int k=0; k<3; k++)
	{
	  x[k] = mesh->GetVertex(ind[k], 0);
	  y[k] = mesh->GetVertex(ind[k], 1);	  
	}
      
      double t1 = x[1] - x[0];
      double t2 = x[2] - x[0];
      double t3 = y[1] - y[0];
      double t4 = y[2] - y[0];

      double s1 = y[1] - y[2];
      double s2 = y[2] - y[0];
      double s3 = y[0] - y[1];
      double s4 = x[2] - x[1];
      double s5 = x[0] - x[2];
      double s6 = x[1] - x[0];
      double detJ = t1 * t4 - t2 * t3;

      double gpcoord[2];
      for(int k=0; k<4; k++)
	{
	  double xx = gp[2*k];
	  double yy = gp[2*k+1];
	  gpcoord[0] = x[0] + t1*xx + t2*yy;
	  gpcoord[1] = y[0] + t3*xx + t4*yy;
 
	  uh = ( approxsol(ind[0])*s1 + approxsol(ind[1])*s2 + approxsol(ind[2])*s3 ) / detJ;
	  uh -= exactderivative[0]->Eval(gpcoord);
	  uh = uh*uh;
	  tr = ( approxsol(ind[0])*s4 + approxsol(ind[1])*s5 + approxsol(ind[2])*s6 ) / detJ;
	  tr -= exactderivative[1]->Eval(gpcoord);
	  tr = tr*tr;
	  herr += ( uh + tr ) * gw[k] * detJ;	
	}
    }  

  return sqrt(herr);
  //cout << "The H1 semi-norm error is: " << sqrt(herr) << endl;
}

//==============================================================================

void LinearFEM2D::ComputePPH1Error(Array<double> &errs, Array<Function *> &exactgrad, Vector &femsol, Array<double *> &ppsol)
{
  errs[0] = 0.0;
  errs[1] = 0.0;
  errs[2] = 0.0;

  double gp[8] = { 1.0/3.0, 1.0/3.0, 0.2, 0.6, 0.2, 0.2, 0.6, 0.2 };
  double gw[4] = { -27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0 };
  Array<double> x(3);
  Array<double> y(3);
   
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      mesh->GetElementVertices(i, ind);
       
      for(int k=0; k<3; k++)
	{
	  x[k] = mesh->GetVertex(ind[k], 0);
	  y[k] = mesh->GetVertex(ind[k], 1);	  
	}
      
      double t1 = x[1] - x[0];
      double t2 = x[2] - x[0];
      double t3 = y[1] - y[0];
      double t4 = y[2] - y[0];

      double s1 = y[1] - y[2];
      double s2 = y[2] - y[0];
      double s3 = y[0] - y[1];
      double s4 = x[2] - x[1];
      double s5 = x[0] - x[2];
      double s6 = x[1] - x[0];
      double detJ = t1 * t4 - t2 * t3;

      double gpcoord[2];
      for(int k=0; k<4; k++)
	{
	  double xx = gp[2*k];
	  double yy = gp[2*k+1];
	  gpcoord[0] = x[0] + t1*xx + t2*yy;
	  gpcoord[1] = y[0] + t3*xx + t4*yy;

	  double ux = exactgrad[0]->Eval(gpcoord);
	  double uy = exactgrad[1]->Eval(gpcoord);

	  double femuhx = 0.0;
	  double femuhy = 0.0;
	  femuhx = ( femsol(ind[0])*s1 + femsol(ind[1])*s2 + femsol(ind[2])*s3 ) / detJ;
	  femuhy = ( femsol(ind[0])*s4 + femsol(ind[1])*s5 + femsol(ind[2])*s6 ) / detJ;

	  double ppuhx = 0.0;
	  double ppuhy = 0.0;
	  ppuhx = ( ppsol[0][i]*s1 + ppsol[1][i]*s2 + ppsol[2][i]*s3 ) / detJ;
	  ppuhy = ( ppsol[0][i]*s4 + ppsol[1][i]*s5 + ppsol[2][i]*s6 ) / detJ;
	  
	  double temx = ux - femuhx;
	  double temy = uy - femuhy;
	  errs[0] += gw[k] * (temx*temx + temy*temy) * detJ;

	  temx = ux - ppuhx;
	  temy = uy - ppuhy;
	  errs[1] += gw[k] * (temx*temx + temy*temy) * detJ;
	  
	  temx = femuhx - ppuhx;
	  temy = femuhy - ppuhy;
	  errs[2] += gw[k] * (temx*temx + temy*temy) * detJ;
 	}
    }

  errs[0] = sqrt(errs[0]);
  errs[1] = sqrt(errs[1]);
  errs[2] = sqrt(errs[2]);
}
