#include "fem_header.h"
#include "../general/general_header.h"

Postprocessing::Postprocessing(DualMesh *_dualmesh, FEM *_fem, Data *_data)
{
  dualmesh = _dualmesh;
  mesh = dualmesh->GetMesh();
  fem = _fem;
  data = _data;
  if(mesh->GetMType()==Element::TRIANGLE)
    {
      if( dualmesh->GetDualMeshOrder()==1 )
	NumLocalDOF = 3;
      else if( dualmesh->GetDualMeshOrder()==2 )
	NumLocalDOF = 6;
      else
	NumLocalDOF = 10;
    }
  else
    {
      if( dualmesh->GetDualMeshOrder()==1 )
	NumLocalDOF = 4;
      else if( dualmesh->GetDualMeshOrder()==2 )
	NumLocalDOF = 9;
      else
	NumLocalDOF = 16;
    }
}

//=========================================================

void Postprocessing::ComputeConservativeFlux(Vector &tsol, Vector &sol, Array<double*> &flux, Array<double *> &ppsol,
					     Array<double *> &ppflux, Array<double> &ppbdrflux)
{
  int n = NumLocalDOF;
  int NE = mesh->GetNE();

  for(int j=0; j<ppbdrflux.Size(); j++) ppbdrflux[j] = 0.0;
  
  for(int j=0; j<n; j++)
    for(int i=0; i<NE; i++)
      ppsol[j][i] = 0.0;
	  
  for(int j=0; j<flux.Size(); j++)
    {
      for(int i=0; i<NE; i++)
	{
	  flux[j][i] = 0.0;
	  ppflux[j][i] = 0.0;
	}
    }
  
  if(dualmesh->GetDualMeshOrder()==1)
    {
      if(mesh->GetMType()==Element::TRIANGLE)
	ComputeLinearConservativeTriangularFlux(tsol, sol, flux, ppsol, ppflux, ppbdrflux);
      else
	ComputeLinearConservativeRectangularFlux(tsol, sol, flux, ppsol, ppflux, ppbdrflux);
    }
  else if(dualmesh->GetDualMeshOrder()==2)
    {
      if(mesh->GetMType()==Element::TRIANGLE)
	ComputeQuadraticConservativeTriangularFlux(tsol, sol, flux, ppsol, ppflux, ppbdrflux);
      else
	ComputeQuadraticConservativeRectangularFlux(tsol, sol, flux, ppsol, ppflux, ppbdrflux);
    }
  else
    {
      if(mesh->GetMType()==Element::TRIANGLE)
	ComputeCubicConservativeTriangularFlux(tsol, sol, flux, ppsol, ppflux, ppbdrflux);
      else
	ComputeCubicConservativeRectangularFlux(tsol, sol, flux, ppsol, ppflux, ppbdrflux);
    }
}

//===========================================================

void Postprocessing::ComputeLinearConservativeTriangularFlux(Vector &tsol, Vector &sol, Array<double*> &flux, Array<double*> &ppsol,
							     Array<double*> &ppflux, Array<double> &ppbdrflux)
{
  mesh->GetEdgeToElementTable();

  Array<double *> normals(3);
  Array<double> elengths(3);
  Array<double *> ecoord(4);
  Array<double> areas(3);  
  for(int k=0; k<ecoord.Size(); k++) ecoord[k] = new double[2];
  for(int k=0; k<normals.Size(); k++) normals[k] = new double[2];

  Array<double *> sflux(mesh->GetNE()); // boundary flux without phi for averaging
  for(int i=0; i<sflux.Size(); i++) sflux[i] = new double[6];
  for(int j=0; j<sflux.Size(); j++)
    for(int i=0; i<6; i++)
      sflux[j][i] = 0.0;

  Array<double *> tflux(mesh->GetNE()); // boundary flux with phi for averaging
  for(int i=0; i<tflux.Size(); i++) tflux[i] = new double[6];
  for(int j=0; j<tflux.Size(); j++)
    for(int i=0; i<6; i++)
      tflux[j][i] = 0.0;

  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      mesh->GetElementVertices(i, ind);

      Array<double *> coord(mesh->GetDim());
      for (int k=0; k<coord.Size(); k++) coord[k] = new double[3];
      
      mesh->GetElementVerticesCoord(i, coord);
      double detJ = (coord[0][1] - coord[0][0]) * (coord[1][2] - coord[1][0])
	- (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);
      
      Vector locsol(3);
      for(int k=0; k<3; k++) locsol(k) = sol(ind[k]);

      Array<double> locdk(ind.Size());
      for (int j=0; j<ind.Size(); j++)
	locdk[j] = data->GetNodalEllipticCoeff(ind[j]);
      
      Array<double *> btem(6);
      for(int k=0; k<6; k++) btem[k] = new double[3];
      for(int k=0; k<6; k++)
	for(int j=0; j<3; j++)
	  btem[k][j] = 0.0;

      double tem0 = sqrt( (coord[0][1] - coord[0][0])*(coord[0][1] - coord[0][0]) +
			  (coord[1][1] - coord[1][0])*(coord[1][1] - coord[1][0]) ) / detJ;
      
      double tem1 = sqrt( (coord[0][1] - coord[0][2])*(coord[0][1] - coord[0][2]) +
			  (coord[1][1] - coord[1][2])*(coord[1][1] - coord[1][2]) ) / detJ;
      
      double tem2 = sqrt( (coord[0][2] - coord[0][0])*(coord[0][2] - coord[0][0]) +
			  (coord[1][2] - coord[1][0])*(coord[1][2] - coord[1][0]) ) / detJ;
      
      if(i%2==0)
	{
	  btem[0][0] = ( 3.0*locdk[0] + locdk[1] ) * ( coord[0][1] - coord[0][2] ) * tem0 / 8.0;
	  btem[0][1] = ( 3.0*locdk[0] + locdk[1] ) * ( coord[0][2] - coord[0][0] ) * tem0 / 8.0;
	  btem[0][2] = ( 3.0*locdk[0] + locdk[1] ) * ( coord[0][0] - coord[0][1] ) * tem0 / 8.0;
      
	  btem[1][0] = ( 3.0*locdk[1] + locdk[0] ) * ( coord[0][1] - coord[0][2] ) * tem0 / 8.0;
	  btem[1][1] = ( 3.0*locdk[1] + locdk[0] ) * ( coord[0][2] - coord[0][0] ) * tem0 / 8.0;
	  btem[1][2] = ( 3.0*locdk[1] + locdk[0] ) * ( coord[0][0] - coord[0][1] ) * tem0 / 8.0;
      
	  btem[2][0] = ( 3.0*locdk[1] + locdk[2] ) * ( coord[1][1] - coord[1][2] + coord[0][2] - coord[0][1] ) * tem1 * sqrt(2.0) / 16.0; 
	  btem[2][1] = ( 3.0*locdk[1] + locdk[2] ) * ( coord[1][2] - coord[1][0] + coord[0][0] - coord[0][2] ) * tem1 * sqrt(2.0) / 16.0;
	  btem[2][2] = ( 3.0*locdk[1] + locdk[2] ) * ( coord[1][0] - coord[1][1] + coord[0][1] - coord[0][0] ) * tem1 * sqrt(2.0) / 16.0;
     
	  btem[3][0] = ( 3.0*locdk[2] + locdk[1] ) * ( coord[1][1] - coord[1][2] + coord[0][2] - coord[0][1] ) * tem1 * sqrt(2.0) / 16.0;
	  btem[3][1] = ( 3.0*locdk[2] + locdk[1] ) * ( coord[1][2] - coord[1][0] + coord[0][0] - coord[0][2] ) * tem1 * sqrt(2.0) / 16.0;
	  btem[3][2] = ( 3.0*locdk[2] + locdk[1] ) * ( coord[1][0] - coord[1][1] + coord[0][1] - coord[0][0] ) * tem1 * sqrt(2.0) / 16.0;
            
	  btem[4][0] = ( 3.0*locdk[2] + locdk[0] ) * ( coord[1][2] - coord[1][1] ) * tem2 / 8.0;
	  btem[4][1] = ( 3.0*locdk[2] + locdk[0] ) * ( coord[1][0] - coord[1][2] ) * tem2 / 8.0;
	  btem[4][2] = ( 3.0*locdk[2] + locdk[0] ) * ( coord[1][1] - coord[1][0] ) * tem2 / 8.0;

	  btem[5][0] = ( 3.0*locdk[0] + locdk[2] ) * ( coord[1][2] - coord[1][1] ) * tem2 / 8.0;
	  btem[5][1] = ( 3.0*locdk[0] + locdk[2] ) * ( coord[1][0] - coord[1][2] ) * tem2 / 8.0;
	  btem[5][2] = ( 3.0*locdk[0] + locdk[2] ) * ( coord[1][1] - coord[1][0] ) * tem2 / 8.0;
	}
      else
	{
	  btem[0][0] = ( 3.0*locdk[0] + locdk[1] ) * ( coord[1][1] - coord[1][2] ) * tem0 / 8.0;
	  btem[0][1] = ( 3.0*locdk[0] + locdk[1] ) * ( coord[1][2] - coord[1][0] ) * tem0 / 8.0;
	  btem[0][2] = ( 3.0*locdk[0] + locdk[1] ) * ( coord[1][0] - coord[1][1] ) * tem0 / 8.0;
      
	  btem[1][0] = ( 3.0*locdk[1] + locdk[0] ) * ( coord[1][1] - coord[1][2] ) * tem0 / 8.0;
	  btem[1][1] = ( 3.0*locdk[1] + locdk[0] ) * ( coord[1][2] - coord[1][0] ) * tem0 / 8.0;
	  btem[1][2] = ( 3.0*locdk[1] + locdk[0] ) * ( coord[1][0] - coord[1][1] ) * tem0 / 8.0;
      
	  btem[4][0] = ( 3.0*locdk[2] + locdk[0] ) * ( coord[0][1] - coord[0][2] + coord[1][2] - coord[1][1] ) * tem1 * sqrt(2.0) / 16.0; 
	  btem[4][1] = ( 3.0*locdk[2] + locdk[0] ) * ( coord[0][2] - coord[0][0] + coord[1][0] - coord[1][2] ) * tem1 * sqrt(2.0) / 16.0;
	  btem[4][2] = ( 3.0*locdk[2] + locdk[0] ) * ( coord[0][0] - coord[0][1] + coord[1][1] - coord[1][0] ) * tem1 * sqrt(2.0) / 16.0;
     
	  btem[5][0] = ( 3.0*locdk[0] + locdk[2] ) * ( coord[0][1] - coord[0][2] + coord[1][2] - coord[1][1] ) * tem1 * sqrt(2.0) / 16.0;
	  btem[5][1] = ( 3.0*locdk[0] + locdk[2] ) * ( coord[0][2] - coord[0][0] + coord[1][0] - coord[1][2] ) * tem1 * sqrt(2.0) / 16.0;
	  btem[5][2] = ( 3.0*locdk[0] + locdk[2] ) * ( coord[0][0] - coord[0][1] + coord[1][1] - coord[1][0] ) * tem1 * sqrt(2.0) / 16.0;
            
	  btem[2][0] = ( 3.0*locdk[1] + locdk[2] ) * ( coord[0][2] - coord[0][1] ) * tem2 / 8.0;
	  btem[2][1] = ( 3.0*locdk[1] + locdk[2] ) * ( coord[0][0] - coord[0][2] ) * tem2 / 8.0;
	  btem[2][2] = ( 3.0*locdk[1] + locdk[2] ) * ( coord[0][1] - coord[0][0] ) * tem2 / 8.0;

	  btem[3][0] = ( 3.0*locdk[2] + locdk[1] ) * ( coord[0][2] - coord[0][1] ) * tem2 / 8.0;
	  btem[3][1] = ( 3.0*locdk[2] + locdk[1] ) * ( coord[0][0] - coord[0][2] ) * tem2 / 8.0;
	  btem[3][2] = ( 3.0*locdk[2] + locdk[1] ) * ( coord[0][1] - coord[0][0] ) * tem2 / 8.0;
	}
      
      for(int k=0; k<6; k++)
	for(int l=0; l<3; l++)
	  sflux[i][k] += btem[k][l] * locsol(l);	  
    
      // boundary integration with phi for averaging
      if(i%2==0)
	{
	  btem[0][0] = ( 2.0*locdk[0] + locdk[1] ) * ( coord[0][1] - coord[0][2] ) * tem0 / 6.0;
	  btem[0][1] = ( 2.0*locdk[0] + locdk[1] ) * ( coord[0][2] - coord[0][0] ) * tem0 / 6.0;
	  btem[0][2] = ( 2.0*locdk[0] + locdk[1] ) * ( coord[0][0] - coord[0][1] ) * tem0 / 6.0;
      
	  btem[1][0] = ( 2.0*locdk[1] + locdk[0] ) * ( coord[0][1] - coord[0][2] ) * tem0 / 6.0;
	  btem[1][1] = ( 2.0*locdk[1] + locdk[0] ) * ( coord[0][2] - coord[0][0] ) * tem0 / 6.0;
	  btem[1][2] = ( 2.0*locdk[1] + locdk[0] ) * ( coord[0][0] - coord[0][1] ) * tem0 / 6.0;

	  btem[2][0] = ( 2.0*locdk[1] + locdk[2] ) * ( coord[1][1] - coord[1][2] + coord[0][2] - coord[0][1] ) * tem1 * sqrt(2.0) / 12.0;
	  btem[2][1] = ( 2.0*locdk[1] + locdk[2] ) * ( coord[1][2] - coord[1][0] + coord[0][0] - coord[0][2] ) * tem1 * sqrt(2.0) / 12.0;
	  btem[2][2] = ( 2.0*locdk[1] + locdk[2] ) * ( coord[1][0] - coord[1][1] + coord[0][1] - coord[0][0] ) * tem1 * sqrt(2.0) / 12.0;

	  btem[3][0] = ( 2.0*locdk[2] + locdk[1] ) * ( coord[1][1] - coord[1][2] + coord[0][2] - coord[0][1] ) * tem1 * sqrt(2.0) / 12.0;
	  btem[3][1] = ( 2.0*locdk[2] + locdk[1] ) * ( coord[1][2] - coord[1][0] + coord[0][0] - coord[0][2] ) * tem1 * sqrt(2.0) / 12.0;
	  btem[3][2] = ( 2.0*locdk[2] + locdk[1] ) * ( coord[1][0] - coord[1][1] + coord[0][1] - coord[0][0] ) * tem1 * sqrt(2.0) / 12.0;
            
	  btem[4][0] = ( 2.0*locdk[2] + locdk[0] ) * ( coord[1][2] - coord[1][1] ) * tem2 / 6.0;
	  btem[4][1] = ( 2.0*locdk[2] + locdk[0] ) * ( coord[1][0] - coord[1][2] ) * tem2 / 6.0;
	  btem[4][2] = ( 2.0*locdk[2] + locdk[0] ) * ( coord[1][1] - coord[1][0] ) * tem2 / 6.0;

	  btem[5][0] = ( 2.0*locdk[0] + locdk[2] ) * ( coord[1][2] - coord[1][1] ) * tem2 / 6.0;
	  btem[5][1] = ( 2.0*locdk[0] + locdk[2] ) * ( coord[1][0] - coord[1][2] ) * tem2 / 6.0;
	  btem[5][2] = ( 2.0*locdk[0] + locdk[2] ) * ( coord[1][1] - coord[1][0] ) * tem2 / 6.0;
	}
      else
	{
	  btem[0][0] = ( 2.0*locdk[0] + locdk[1] ) * ( coord[1][1] - coord[1][2] ) * tem0 / 6.0;
	  btem[0][1] = ( 2.0*locdk[0] + locdk[1] ) * ( coord[1][2] - coord[1][0] ) * tem0 / 6.0;
	  btem[0][2] = ( 2.0*locdk[0] + locdk[1] ) * ( coord[1][0] - coord[1][1] ) * tem0 / 6.0;
      
	  btem[1][0] = ( 2.0*locdk[1] + locdk[0] ) * ( coord[1][1] - coord[1][2] ) * tem0 / 6.0;
	  btem[1][1] = ( 2.0*locdk[1] + locdk[0] ) * ( coord[1][2] - coord[1][0] ) * tem0 / 6.0;
	  btem[1][2] = ( 2.0*locdk[1] + locdk[0] ) * ( coord[1][0] - coord[1][1] ) * tem0 / 6.0;

	  btem[4][0] = ( 2.0*locdk[2] + locdk[0] ) * ( coord[0][1] - coord[0][2] + coord[1][2] - coord[1][1] ) * tem1 * sqrt(2.0) / 12.0;
	  btem[4][1] = ( 2.0*locdk[2] + locdk[0] ) * ( coord[0][2] - coord[0][0] + coord[1][0] - coord[1][2] ) * tem1 * sqrt(2.0) / 12.0;
	  btem[4][2] = ( 2.0*locdk[2] + locdk[0] ) * ( coord[0][0] - coord[0][1] + coord[1][1] - coord[1][0] ) * tem1 * sqrt(2.0) / 12.0;

	  btem[5][0] = ( 2.0*locdk[0] + locdk[2] ) * ( coord[0][1] - coord[0][2] + coord[1][2] - coord[1][1] ) * tem1 * sqrt(2.0) / 12.0;
	  btem[5][1] = ( 2.0*locdk[0] + locdk[2] ) * ( coord[0][2] - coord[0][0] + coord[1][0] - coord[1][2] ) * tem1 * sqrt(2.0) / 12.0;
	  btem[5][2] = ( 2.0*locdk[0] + locdk[2] ) * ( coord[0][0] - coord[0][1] + coord[1][1] - coord[1][0] ) * tem1 * sqrt(2.0) / 12.0;
            
	  btem[2][0] = ( 2.0*locdk[1] + locdk[2] ) * ( coord[0][2] - coord[0][1] ) * tem2 / 6.0;
	  btem[2][1] = ( 2.0*locdk[1] + locdk[2] ) * ( coord[0][0] - coord[0][2] ) * tem2 / 6.0;
	  btem[2][2] = ( 2.0*locdk[1] + locdk[2] ) * ( coord[0][1] - coord[0][0] ) * tem2 / 6.0;

	  btem[3][0] = ( 2.0*locdk[2] + locdk[1] ) * ( coord[0][2] - coord[0][1] ) * tem2 / 6.0;
	  btem[3][1] = ( 2.0*locdk[2] + locdk[1] ) * ( coord[0][0] - coord[0][2] ) * tem2 / 6.0;
	  btem[3][2] = ( 2.0*locdk[2] + locdk[1] ) * ( coord[0][1] - coord[0][0] ) * tem2 / 6.0;
	}

      for(int k=0; k<6; k++)
	for(int l=0; l<3; l++)
	  tflux[i][k] += btem[k][l] * locsol(l);	  
    
      for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
      for (int k=0; k<btem.Size(); k++) { delete []btem[k]; }
    }

  Vector cverr(mesh->GetNV()); cverr = 0.0;
  Vector cverruh(mesh->GetNV()); cverruh = 0.0;

  double bmaxuherr = 0.0;
  double bmaxpuherr = 0.0;
  double bl2uherr = 0.0;
  double bl2puherr = 0.0;

  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      mesh->GetElementVertices(i, ind);
      
      Array<double *> coord(mesh->GetDim());
      for (int k=0; k<coord.Size(); k++) coord[k] = new double[3];
      
      mesh->GetElementVerticesCoord(i, coord);
      double detJ = (coord[0][1] - coord[0][0]) * (coord[1][2] 
	    - coord[1][0]) - (coord[1][1] - coord[1][0]) * (coord[0][2] - coord[0][0]);

      dualmesh->GetElemDualInfo(i, ecoord, normals, elengths, areas);

      //cout << "amove " << ecoord[0][0] << " " << ecoord[0][1] << endl;
      //cout << "aline " << ecoord[3][0] << " " << ecoord[3][1] << endl;
      //cout << "aline " << ecoord[1][0] << " " << ecoord[1][1] << endl;
      //cout << "amove " << ecoord[3][0] << " " << ecoord[3][1] << endl;
      //cout << "aline " << ecoord[2][0] << " " << ecoord[2][1] << endl;	  
      //cout << endl;

      Array<double *> tem(3);
      for(int k=0; k<3; k++) tem[k] = new double[3];
      for(int k=0; k<3; k++)
	for(int j=0; j<3; j++)
	  tem[k][j] = 0.0;

      Array<double> locdk(ind.Size());
      for (int j=0; j<ind.Size(); j++)
	locdk[j] = data->GetNodalEllipticCoeff(ind[j]);
      
      Array<double> locdf(ind.Size());
      for (int j=0; j<ind.Size(); j++)
	locdf[j] = data->GetNodalForceCoeff(ind[j]);

      double tem0 = sqrt( (coord[0][0] + coord[0][1] - 2.0*coord[0][2])*(coord[0][0] + coord[0][1] - 2.0*coord[0][2]) + 
			  (coord[1][0] + coord[1][1] - 2.0*coord[1][2])*(coord[1][0] + coord[1][1] - 2.0*coord[1][2]) ) / detJ;
      
      double tem1 = sqrt( (coord[0][2] + coord[0][1] - 2.0*coord[0][0])*(coord[0][2] + coord[0][1] - 2.0*coord[0][0]) + 
			  (coord[1][2] + coord[1][1] - 2.0*coord[1][0])*(coord[1][2] + coord[1][1] - 2.0*coord[1][0]) ) / detJ;
      
      double tem2 = sqrt( (coord[0][0] + coord[0][2] - 2.0*coord[0][1])*(coord[0][0] + coord[0][2] - 2.0*coord[0][1]) + 
			  (coord[1][0] + coord[1][2] - 2.0*coord[1][1])*(coord[1][0] + coord[1][2] - 2.0*coord[1][1]) ) / detJ;
      
      tem[0][0] = ( 2.5*locdk[0] + 2.5*locdk[1] + locdk[2] ) * ( (coord[1][1] - coord[1][2])*normals[0][0] + (coord[0][2] - coord[0][1])*normals[0][1] ) * tem0 / 36.0;
      tem[0][1] = ( 2.5*locdk[0] + 2.5*locdk[1] + locdk[2] ) * ( (coord[1][2] - coord[1][0])*normals[0][0] + (coord[0][0] - coord[0][2])*normals[0][1] ) * tem0 / 36.0;
      tem[0][2] = ( 2.5*locdk[0] + 2.5*locdk[1] + locdk[2] ) * ( (coord[1][0] - coord[1][1])*normals[0][0] + (coord[0][1] - coord[0][0])*normals[0][1] ) * tem0 / 36.0;

      tem[1][0] = ( 2.5*locdk[2] + 2.5*locdk[1] + locdk[0] ) * ( (coord[1][1] - coord[1][2])*normals[1][0] + (coord[0][2] - coord[0][1])*normals[1][1] ) * tem1 / 36.0;
      tem[1][1] = ( 2.5*locdk[2] + 2.5*locdk[1] + locdk[0] ) * ( (coord[1][2] - coord[1][0])*normals[1][0] + (coord[0][0] - coord[0][2])*normals[1][1] ) * tem1 / 36.0;
      tem[1][2] = ( 2.5*locdk[2] + 2.5*locdk[1] + locdk[0] ) * ( (coord[1][0] - coord[1][1])*normals[1][0] + (coord[0][1] - coord[0][0])*normals[1][1] ) * tem1 / 36.0;

      tem[2][0] = ( 2.5*locdk[0] + 2.5*locdk[2] + locdk[1] ) * ( (coord[1][1] - coord[1][2])*normals[2][0] + (coord[0][2] - coord[0][1])*normals[2][1] ) * tem2 / 36.0;
      tem[2][1] = ( 2.5*locdk[0] + 2.5*locdk[2] + locdk[1] ) * ( (coord[1][2] - coord[1][0])*normals[2][0] + (coord[0][0] - coord[0][2])*normals[2][1] ) * tem2 / 36.0;
      tem[2][2] = ( 2.5*locdk[0] + 2.5*locdk[2] + locdk[1] ) * ( (coord[1][0] - coord[1][1])*normals[2][0] + (coord[0][1] - coord[0][0])*normals[2][1] ) * tem2 / 36.0;
      

      SparseMatrix *locmat = new SparseMatrix(3,3);
      locmat->Elem(0,0) = tem[2][0] - tem[0][0];
      locmat->Elem(0,1) = tem[2][1] - tem[0][1];
      locmat->Elem(0,2) = tem[2][2] - tem[0][2];
      
      locmat->Elem(1,0) = tem[0][0] - tem[1][0];
      locmat->Elem(1,1) = tem[0][1] - tem[1][1];
      locmat->Elem(1,2) = tem[0][2] - tem[1][2];
      
      locmat->Elem(2,0) = tem[1][0] - tem[2][0];
      locmat->Elem(2,1) = tem[1][1] - tem[2][1];
      locmat->Elem(2,2) = tem[1][2] - tem[2][2];
      
      locmat->Finalize();
      
      Vector Q(3);
      Q = 0.0;
      Array<double *> femat(3);
      Vector locsol(3);
      for(int k=0; k<3; k++)
	{
	  femat[k] = new double[3];
	  locsol(k) = sol(ind[k]);
	}
      
      fem->GetEllipticLocalSystem(i, ind, femat);
 
      for(int k=0; k<3; k++)
	for(int j=0; j<3; j++)
	  Q(k) += femat[k][j] * locsol(j);

      Array<double> F(3);
      fem->GetForceLocalSystem(i, ind, F);

      Vector f(3);
      f = 0.0;

      f(0) = ( 22.0*locdf[0] + 7.0*locdf[1] + 7.0*locdf[2] ) * detJ / 216.0;
      f(1) = ( 7.0*locdf[0] + 22.0*locdf[1] + 7.0*locdf[2] ) * detJ / 216.0;
      f(2) = ( 7.0*locdf[0] + 7.0*locdf[1] + 22.0*locdf[2] ) * detJ / 216.0;

      Vector RHS(3);
      for(int k=0; k<3; k++) RHS(k) = Q(k) - F[k] + f(k);

      Array<int> edges;
      mesh->GetElementEdges(i, edges);
      if(i%2==0)
	{
	  Array<int> eind;
	  int k = edges[0];
	  mesh->GetEdgeElements(k, eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[1]][3] - tflux[i][0] + tflux[eind[1]][3] );
		  RHS(1) += 0.5 * ( sflux[i][1] - sflux[eind[1]][2] - tflux[i][1] + tflux[eind[1]][2] );
		}
	      else
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[0]][3] - tflux[i][0] + tflux[eind[0]][3] );
		  RHS(1) += 0.5 * ( sflux[i][1] - sflux[eind[0]][2] - tflux[i][1] + tflux[eind[0]][2] );
		}
	    }
	  else
	    {
	      RHS(0) += -tflux[i][0] + sflux[i][0];
	      RHS(1) += -tflux[i][1] + sflux[i][1];
	    }

	  k = edges[1];
	  mesh->GetEdgeElements(k, eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(1) += 0.5 * ( sflux[i][2] - sflux[eind[1]][5] - tflux[i][2] + tflux[eind[1]][5] );
		  RHS(2) += 0.5 * ( sflux[i][3] - sflux[eind[1]][4] - tflux[i][3] + tflux[eind[1]][4] );
		}
	      else
		{
		  RHS(1) += 0.5 * ( sflux[i][2] - sflux[eind[0]][5] - tflux[i][2] + tflux[eind[0]][5] );
		  RHS(2) += 0.5 * ( sflux[i][3] - sflux[eind[0]][4] - tflux[i][3] + tflux[eind[0]][4] );
		}
	    }
	  else
	    {
	      RHS(1) += -tflux[i][2] + sflux[i][2];
	      RHS(2) += -tflux[i][3] + sflux[i][3];
	    }

	 k = edges[2];
	 mesh->GetEdgeElements(k, eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(0) += 0.5 * ( sflux[i][5] - sflux[eind[1]][0] - tflux[i][5] + tflux[eind[1]][0] );
		  RHS(2) += 0.5 * ( sflux[i][4] - sflux[eind[1]][1] - tflux[i][4] + tflux[eind[1]][1] );
		}
	      else
		{
		  RHS(0) += 0.5 * ( sflux[i][5] - sflux[eind[0]][0] - tflux[i][5] + tflux[eind[0]][0] );
		  RHS(2) += 0.5 * ( sflux[i][4] - sflux[eind[0]][1] - tflux[i][4] + tflux[eind[0]][1] );
		}
	    }
	  else
	    {
	      RHS(0) += -tflux[i][5] + sflux[i][5];
	      RHS(2) += -tflux[i][4] + sflux[i][4];
	    }
	}
      else
	{
	  Array<int> eind;
	  int k = edges[0];
	  mesh->GetEdgeElements(k, eind);
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[1]][5] - tflux[i][0] + tflux[eind[1]][5] );
		  RHS(1) += 0.5 * ( sflux[i][1] - sflux[eind[1]][4] - tflux[i][1] + tflux[eind[1]][4] );
		}
	      else
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[0]][5] - tflux[i][0] + tflux[eind[0]][5] );
		  RHS(1) += 0.5 * ( sflux[i][1] - sflux[eind[0]][4] - tflux[i][1] + tflux[eind[0]][4] );
		}
	    }
	  else
	    {
	      RHS(0) += -tflux[i][0] + sflux[i][0];
	      RHS(1) += -tflux[i][1] + sflux[i][1];
	    }

	  k = edges[1];
	  mesh->GetEdgeElements(k, eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(1) += 0.5 * ( sflux[i][2] - sflux[eind[1]][1] - tflux[i][2] + tflux[eind[1]][1] );
		  RHS(2) += 0.5 * ( sflux[i][3] - sflux[eind[1]][0] - tflux[i][3] + tflux[eind[1]][0] );
		}
	      else
		{
		  RHS(1) += 0.5 * ( sflux[i][2] - sflux[eind[0]][1] - tflux[i][2] + tflux[eind[0]][1] );
		  RHS(2) += 0.5 * ( sflux[i][3] - sflux[eind[0]][0] - tflux[i][3] + tflux[eind[0]][0] );
		}
	    }
	  else
	    {
	      RHS(1) += -tflux[i][2] + sflux[i][2];
	      RHS(2) += -tflux[i][3] + sflux[i][3];
	    }

	 k = edges[2];
	 mesh->GetEdgeElements(k, eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(0) += 0.5 * ( sflux[i][5] - sflux[eind[1]][2] - tflux[i][5] + tflux[eind[1]][2] );
		  RHS(2) += 0.5 * ( sflux[i][4] - sflux[eind[1]][3] - tflux[i][4] + tflux[eind[1]][3] );
		}
	      else
		{
		  RHS(0) += 0.5 * ( sflux[i][5] - sflux[eind[0]][2] - tflux[i][5] + tflux[eind[0]][2] );
		  RHS(2) += 0.5 * ( sflux[i][4] - sflux[eind[0]][3] - tflux[i][4] + tflux[eind[0]][3] );
		}
	    }
	  else
	    {
	      RHS(0) += -tflux[i][5] + sflux[i][5];
	      RHS(2) += -tflux[i][4] + sflux[i][4];
	    }
	}
            
      Vector s(3);
      s = 0.0;
      Vector rhs(3);
      rhs = RHS;

      RectangularMatrix localconsmat(3);
      
      for(int j=0; j<3; j++)
	for(int k=0; k<3; k++)
	  localconsmat(j,k) = locmat->Elem(j,k);

      for(int j=1; j<3; j++) localconsmat(0, j) = 0.0;
      localconsmat(0,0) = 1.0;
      RHS(0) = locsol(0); 

      DenseMatrixInverse invmat(localconsmat);
      invmat.Mult(RHS, s);

      Vector sss(3);
      sss = 0.0;
      locmat->Mult(s, sss);
      for(int k=0; k<3; k++) cverr(ind[k]) += sss(k) - f(k); 

      locmat->Mult(locsol, sss);
      for(int k=0; k<3; k++) cverruh(ind[k]) += sss(k) - f(k); 

      Vector trueflux(3);
      trueflux = 0.0;
        
      for(int k=0; k<3; k++) { ppsol[k][i] = s(k); }
      for(int k=0; k<3; k++)
	{
	  for(int j=0; j<3; j++)
	    {
	      ppflux[k][i] += -tem[k][j] * s(j);
	      flux[k][i] += -tem[k][j] * locsol(j);
	      trueflux(k) += -tem[k][j] * tsol(ind[j]);
	    }
	  
	  //cout << i << "\t" << ppflux[k][i] << endl;
	  
	  bl2uherr += fabs( trueflux(k) - flux[k][i] ) * fabs( trueflux(k) - flux[k][i] );
	  if( fabs( trueflux(k) - flux[k][i] ) > bmaxuherr )
	    bmaxuherr = fabs( trueflux(k) - flux[k][i] );

	  bl2puherr += fabs( trueflux(k) - ppflux[k][i] ) * fabs( trueflux(k) - ppflux[k][i] );
	  if( fabs( trueflux(k) - ppflux[k][i] ) > bmaxpuherr )
	    bmaxpuherr = fabs( trueflux(k) - ppflux[k][i] );
	}
      
      for(int k=0; k<femat.Size(); k++) delete []femat[k];
      for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
      for (int k=0; k<tem.Size(); k++) { delete []tem[k]; }
      delete locmat;
    }      
 
  for(int i=0; i<mesh->GetNBE(); i++)
    {
      Array<int> bdrind(2);
      mesh->GetBdrElementVertices(i, bdrind);
      if( fabs(cverr(bdrind[0])) > 1.0e-16 )
	{
	  ppbdrflux[bdrind[0]] = -cverr(bdrind[0]);
	  cverr(bdrind[0]) = 0.0;
	  cverruh(bdrind[0]) = 0.0;	  
	}
      if( fabs(cverr(bdrind[1])) > 1.0e-16 )
	{
	  ppbdrflux[bdrind[1]] = -cverr(bdrind[1]);
	  cverr(bdrind[1]) = 0.0;
	  cverruh(bdrind[1]) = 0.0;
	}
    }

  double maxlce = 0.0;
  ofstream fileout("llce.out");
  for(int i=0; i<dualmesh->GetDualMeshNumDOF(); i++)
    {
      //fileout<<i<<"\t"<<cverruh(i)<< "\t" << cverr(i) << endl;
      if(fabs(cverr(i)) > maxlce) 
	maxlce = fabs(cverr(i));
      //cout<<setprecision(6)<<i<<"\t"<<cverruh(i)<< "\t" <<cverr(i)<<endl;
    }
  fileout.close(); 
  cout << "Max LCE: " << maxlce << endl;
 

  //cout << "Edge Flux Max Error  uh: " << bmaxuherr << endl;
  //cout << "Edge Flux Max Error puh: " << bmaxpuherr << endl;
  //cout << "Edge Flux L2 Error   uh: " << sqrt(bl2uherr) << endl;
  //cout << "Edge Flux L2 Error  puh: " << sqrt(bl2puherr) << endl;
  
  for(int i=0; i<tflux.Size(); i++) delete []tflux[i];
  for(int i=0; i<sflux.Size(); i++) delete []sflux[i];  
  for(int i=0; i<normals.Size(); i++) delete []normals[i];  
  for(int i=0; i<ecoord.Size(); i++) delete []ecoord[i];  
}

//===========================================================

void Postprocessing::ComputeQuadraticConservativeTriangularFlux(Vector &tsol, Vector &sol, Array<double*> &flux, Array<double*> &ppsol,
								Array<double*> &ppflux, Array<double> &ppbdrflux)
{
  mesh->GetEdgeToElementTable();
    
  Array<double *> nmls(3);
  Array<double> lens(3);
  Array<double *> cod(13);
  Array<double> areas(1);  
  for(int k=0; k<cod.Size(); k++) cod[k] = new double[2];
  for(int k=0; k<nmls.Size(); k++) nmls[k] = new double[2];

  int n = 4;
  Array<double *> pw(n);
  for(int i=0; i<n; i++) pw[i] = new double[2];
  GetStandLineQuadPW(n, pw);

  int np = ConvertN2PN(5);
  Array<double *> tpw(np);
  for(int i=0; i<np; i++) tpw[i] = new double[3];
  GetStandTriQuadPW(5, tpw);

  Array<double *> sflux(mesh->GetNE());
  for(int i=0; i<sflux.Size(); i++) sflux[i] = new double[9];
  for(int j=0; j<sflux.Size(); j++)
    for(int i=0; i<9; i++)
      sflux[j][i] = 0.0;

  Array<double *> tflux(mesh->GetNE());
  for(int i=0; i<tflux.Size(); i++) tflux[i] = new double[9];
  for(int j=0; j<tflux.Size(); j++)
    for(int i=0; i<9; i++)
      tflux[j][i] = 0.0;

  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      fem->getElementDOF(i, ind);

      Array<double *> coord(mesh->GetDim());
      for (int k=0; k<coord.Size(); k++) coord[k] = new double[3];
      
      mesh->GetElementVerticesCoord(i, coord);

      double t1 = coord[0][1] - coord[0][0];
      double t2 = coord[0][2] - coord[0][0];
      double t3 = coord[1][1] - coord[1][0];
      double t4 = coord[1][2] - coord[1][0];
      double detJ = t1 * t4 - t2 * t3;
      
      dualmesh->GetElemDualInfo(i, cod, nmls, lens, areas); 

      Vector locsol(6);
      for(int k=0; k<6; k++) locsol(k) = sol(ind[k]);
      //locsol.Print();

      Array<double> locdk(ind.Size());
      for (int j=0; j<ind.Size(); j++)
	locdk[j] = data->GetNodalEllipticCoeff(ind[j]);
      
      Array<double *> btem(9);
      for(int k=0; k<9; k++) btem[k] = new double[6];
      for(int k=0; k<9; k++)
	for(int j=0; j<6; j++)
	  btem[k][j] = 0.0;

      //========= boundary flux without phi for averaging ==========================
      double tem0 = sqrt( (coord[0][1] - coord[0][0])*(coord[0][1] - coord[0][0]) +
			  (coord[1][1] - coord[1][0])*(coord[1][1] - coord[1][0]) ) / detJ;
      
      double tem1 = sqrt( (coord[0][1] - coord[0][2])*(coord[0][1] - coord[0][2]) +
			  (coord[1][1] - coord[1][2])*(coord[1][1] - coord[1][2]) ) / detJ;
      
      double tem2 = sqrt( (coord[0][2] - coord[0][0])*(coord[0][2] - coord[0][0]) +
			  (coord[1][2] - coord[1][0])*(coord[1][2] - coord[1][0]) ) / detJ;

      double xa, xb, ya, yb, m0, m1, m2;
      double len = 0.0;
      double nn[2];

      //---0--------------------------
      xa = coord[0][0];
      xb = cod[0][0];
      ya = coord[1][0];
      yb = cod[0][1];
      len = sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya) );
      nn[0] = (yb - ya)/len;
      nn[1] = (xa - xb)/len;
      m0 = ( coord[1][2] - coord[1][1] ) * nn[0] + ( coord[0][1] - coord[0][2] ) * nn[1];
      m1 = ( coord[1][2] - coord[1][0] ) * nn[0] + ( coord[0][0] - coord[0][2] ) * nn[1];
      m2 = ( coord[1][0] - coord[1][1] ) * nn[0] + ( coord[0][1] - coord[0][0] ) * nn[1];

      btem[0][0] = ( -55.0*locdk[0]/128.0 + 19.0*locdk[1]/384.0 - 47.0*locdk[3]/192.0 ) * m0 * tem0;
      btem[0][1] = ( -37.0*locdk[0]/384.0 + locdk[1]/128.0 - 7.0*locdk[3]/192.0 ) * m1 * tem0;
      btem[0][2] = ( -locdk[0]/6.0 + locdk[1]/48.0 - 5.0*locdk[3]/48.0 ) * m2 * tem0;
      btem[0][3] = ( 101.0*locdk[0]/192.0 - 11.0*locdk[1]/192.0 + 9.0*locdk[3]/32.0 ) * m1 * tem0
	+ ( -9.0*locdk[0]/128.0 + 5.0*locdk[1]/384.0 - 13.0*locdk[3]/192.0 ) * m2 * tem0;
      btem[0][4] = -( -9.0*locdk[0]/128.0 + 5.0*locdk[1]/384.0 - 13.0*locdk[3]/192.0 ) * m2 * tem0;
      btem[0][5] = ( 229.0*locdk[0]/384.0 - 9.0*locdk[1]/128.0 + 67.0*locdk[3]/192.0 ) * m2 * tem0;

      btem[1][0] = ( -locdk[0]/16.0 + locdk[1]/48.0 - 11.0*locdk[3]/24.0 ) * m0 * tem0;
      btem[1][1] = ( -locdk[0]/48.0 + locdk[1]/16.0 + 11.0*locdk[3]/24.0 ) * m1 * tem0;
      btem[1][2] = ( -locdk[0]/48.0 - locdk[1]/48.0 - 11.0*locdk[3]/24.0 ) * m2 * tem0;      
      btem[1][3] = ( locdk[0]/12.0 - locdk[1]/12.0 ) * m1 * tem0
	+ ( -locdk[1]/12.0 - 11.0*locdk[3]/12.0 ) * m2 * tem0;
      btem[1][4] = -( -locdk[1]/12.0 - 11.0*locdk[3]/12.0 ) * m2 * tem0;
      btem[1][5] = ( locdk[0]/12.0 + 11.0*locdk[3]/12.0 ) * m2 * tem0;

      btem[2][0] = ( -locdk[0]/128.0 + 37.0*locdk[1]/384.0 + 7.0*locdk[3]/192.0 ) * m0 * tem0;
      btem[2][1] = ( -19.0*locdk[0]/384.0 + 55.0*locdk[1]/128.0 + 47.0*locdk[3]/192.0 ) * m1 * tem0;
      btem[2][2] = ( locdk[0]/48.0 - locdk[1]/6.0 - 5.0*locdk[3]/48.0 ) * m2 * tem0;
      btem[2][3] = ( 11.0*locdk[0]/192.0 - 101.0*locdk[1]/192.0 - 9.0*locdk[3]/32.0 ) * m1 * tem0
	+ ( 9.0*locdk[0]/128.0 - 229.0*locdk[1]/384.0 - 67.0*locdk[3]/192.0 ) * m2 * tem0;
      btem[2][4] = -( 9.0*locdk[0]/128.0 - 229.0*locdk[1]/384.0 - 67.0*locdk[3]/192.0 ) * m2 * tem0;
      btem[2][5] = ( -5.0*locdk[0]/384.0 + 9.0*locdk[1]/128.0 + 13.0*locdk[3]/192.0 ) * m2 * tem0;
      
      //---3----------------------------      
      xa = coord[0][1];
      xb = cod[4][0];
      ya = coord[1][1];
      yb = cod[4][1];
      len = sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya) );
      nn[0] = (yb - ya)/len;
      nn[1] = (xa - xb)/len;
      m0 = ( coord[1][2] - coord[1][1] ) * nn[0] + ( coord[0][1] - coord[0][2] ) * nn[1];
      m1 = ( coord[1][2] - coord[1][0] ) * nn[0] + ( coord[0][0] - coord[0][2] ) * nn[1];
      m2 = ( coord[1][0] - coord[1][1] ) * nn[0] + ( coord[0][1] - coord[0][0] ) * nn[1];      

      btem[3][0] = ( locdk[1]/6.0 - locdk[2]/48.0 + 5.0*locdk[4]/48.0 ) * m0 * tem1;
      btem[3][1] = ( 55.0*locdk[1]/128.0 - 19.0*locdk[2]/384.0 + 47.0*locdk[4]/192.0 ) * m1 * tem1;
      btem[3][2] = ( -37.0*locdk[1]/384.0 + locdk[2]/128.0 - 7.0*locdk[4]/192.0 ) * m2 * tem1;
      btem[3][3] = ( -229.0*locdk[1]/384.0 + 9.0*locdk[2]/128.0 - 67.0*locdk[4]/192.0 ) * m1 * tem1
	+ ( -229.0*locdk[1]/384.0 + 9.0*locdk[2]/128.0 - 67.0*locdk[4]/192.0 ) * m2 * tem1;
      btem[3][4] = ( 9.0*locdk[1]/128.0 - 5.0*locdk[2]/384.0 + 13.0*locdk[4]/192.0 ) * m1 * tem1
	- ( -229.0*locdk[1]/384.0 + 9.0*locdk[2]/128.0 - 67.0*locdk[4]/192.0 ) * m2 * tem1;
      btem[3][5] = -( 9.0*locdk[1]/128.0 - 5.0*locdk[2]/384.0 + 13.0*locdk[4]/192.0 ) * m1 * tem1
	+ ( -9.0*locdk[1]/128.0 + 5.0*locdk[2]/384.0 - 13.0*locdk[4]/192.0 ) * m2 * tem1;

      btem[4][0] = ( locdk[1]/48.0 + locdk[2]/48.0 + 11.0*locdk[4]/24.0 ) * m0 * tem1;
      btem[4][1] = ( locdk[1]/16.0 - locdk[2]/48.0 + 11.0*locdk[4]/24.0 ) * m1 * tem1;
      btem[4][2] = ( -locdk[1]/48.0 + locdk[2]/16.0 + 11.0*locdk[4]/24.0 ) * m2 * tem1;
      btem[4][3] = -( locdk[1]/12.0 + 11.0*locdk[4]/12.0 ) * m1 * tem1 
	- ( locdk[1]/12.0 + 11.0*locdk[4]/12.0 ) * m2 * tem1;
      btem[4][4] = ( locdk[2]/12.0 + 11.0*locdk[4]/12.0 ) * m1 * tem1
	+ ( locdk[1]/12.0 + 11.0*locdk[4]/12.0 ) * m2 * tem1;
      btem[4][5] = -( locdk[2]/12.0 + 11.0*locdk[4]/12.0 ) * m1 * tem1
	- ( locdk[2]/12.0 + 11.0*locdk[4]/12.0 ) * m2 * tem1;
      
      btem[5][0] = ( -locdk[1]/48.0 + locdk[2]/6.0 + 5.0*locdk[4]/48.0 ) * m0 * tem1;
      btem[5][1] = ( locdk[1]/128.0 - 37.0*locdk[2]/384.0 - 7.0*locdk[4]/192.0 ) * m1 * tem1;
      btem[5][2] = ( -19.0*locdk[1]/384.0 + 55.0*locdk[2]/128.0 + 47.0*locdk[4]/192.0 ) * m2 * tem1;
      btem[5][3] = -( -5.0*locdk[1]/384.0 + 9.0*locdk[2]/128.0 + 13.0*locdk[4]/192.0 ) * m1 * tem1
	- ( -5.0*locdk[1]/384.0 + 9.0*locdk[2]/128.0 + 13.0*locdk[4]/192.0 ) * m2 * tem1;
      btem[5][4] = ( -9.0*locdk[1]/128.0 + 229.0*locdk[2]/384.0 + 67.0*locdk[4]/192.0 ) * m1 * tem1
	+ ( -5.0*locdk[1]/384.0 + 9.0*locdk[2]/128.0 + 13.0*locdk[4]/192.0 ) * m2 * tem1;
      btem[5][5] = -( -9.0*locdk[1]/128.0 + 229.0*locdk[2]/384.0 + 67.0*locdk[4]/192.0 ) * m1 * tem1
	- ( -9.0*locdk[1]/128.0 + 229.0*locdk[2]/384.0 + 67.0*locdk[4]/192.0 ) * m2 * tem1;
      
      //---6----------------------------
      xa = coord[0][2];
      xb = cod[8][0];
      ya = coord[1][2];
      yb = cod[8][1];      
      len = sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya) );
      nn[0] = (yb - ya)/len;
      nn[1] = (xa - xb)/len;
      m0 = ( coord[1][2] - coord[1][1] ) * nn[0] + ( coord[0][1] - coord[0][2] ) * nn[1];
      m1 = ( coord[1][2] - coord[1][0] ) * nn[0] + ( coord[0][0] - coord[0][2] ) * nn[1];
      m2 = ( coord[1][0] - coord[1][1] ) * nn[0] + ( coord[0][1] - coord[0][0] ) * nn[1];
      
      btem[6][0] = ( -locdk[0]/128.0 + 37.0*locdk[2]/384.0 + 7.0*locdk[5]/192.0 ) * m0 * tem2;
      btem[6][1] = ( locdk[0]/48.0 - locdk[2]/6.0 - 5.0*locdk[5]/48.0 ) * m1 * tem2;
      btem[6][2] = ( -19.0*locdk[0]/384.0 + 55.0*locdk[2]/128.0 + 47.0*locdk[5]/192.0 ) * m2 * tem2;
      btem[6][3] = ( -5.0*locdk[0]/384.0 + 9.0*locdk[2]/128.0 + 13.0*locdk[5]/192.0 ) * m1 * tem2;
      btem[6][4] = ( -9.0*locdk[0]/128.0 + 229.0*locdk[2]/384.0 + 67.0*locdk[5]/192.0 ) * m1 * tem2;
      btem[6][5] = -( -9.0*locdk[0]/128.0 + 229.0*locdk[2]/384.0 + 67.0*locdk[5]/192.0 ) * m1 * tem2
	+ ( 11.0*locdk[0]/192.0 - 101.0*locdk[2]/192.0 - 9.0*locdk[5]/32.0 ) * m2 * tem2;

      btem[7][0] = ( -locdk[0]/16.0 + locdk[2]/48.0 - 11.0*locdk[5]/24.0 ) * m0 * tem2;
      btem[7][1] = ( -locdk[0]/48.0 - locdk[2]/48.0 - 11.0*locdk[5]/24.0 ) * m1 * tem2;
      btem[7][2] = ( -locdk[0]/48.0 + locdk[2]/16.0 + 11.0*locdk[5]/24.0 ) * m2 * tem2;
      btem[7][3] = ( locdk[0]/12.0 + 11.0*locdk[5]/12.0 ) * m1 * tem2;
      btem[7][4] = ( locdk[2]/12.0 + 11.0*locdk[5]/12.0 ) * m1 * tem2;
      btem[7][5] = ( locdk[0]/12.0 - locdk[2]/12.0 ) * m2 * tem2
	- ( locdk[2]/12.0 + 11.0*locdk[5]/12.0 ) * m1 * tem2;;
      
      btem[8][0] = ( -55.0*locdk[0]/128.0 + 19.0*locdk[2]/384.0 - 47.0*locdk[5]/192.0 ) * m0 * tem2;
      btem[8][1] = ( -locdk[0]/6.0 + locdk[2]/48.0 - 5.0*locdk[5]/48.0 ) * m1 * tem2;
      btem[8][2] = ( -37.0*locdk[0]/384.0 + locdk[2]/128.0 - 7.0*locdk[5]/192.0 ) * m2 * tem2;
      btem[8][3] = ( 229.0*locdk[0]/384.0 - 9.0*locdk[2]/128.0 + 67.0*locdk[5]/192.0 ) * m1 * tem2;
      btem[8][4] = ( 9.0*locdk[0]/128.0 - 5.0*locdk[2]/384.0 + 13.0*locdk[5]/192.0 ) * m1 * tem2;
      btem[8][5] = -( 9.0*locdk[0]/128.0 - 5.0*locdk[2]/384.0 + 13.0*locdk[5]/192.0 ) * m1 * tem2
	+ ( 101.0*locdk[0]/192.0 - 11.0*locdk[2]/192.0 + 9.0*locdk[5]/32.0 ) * m2 * tem2;
      
      for(int k=0; k<9; k++)
	for(int l=0; l<6; l++)
	  sflux[i][k] += btem[k][l] * locsol(l); 

      //---0-----------------
      xa = coord[0][0];
      xb = coord[0][1];
      ya = coord[1][0];
      yb = coord[1][1];
      len = sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya) );
      nn[0] = (yb - ya)/len;
      nn[1] = (xa - xb)/len;
      m0 = ( coord[1][2] - coord[1][1] ) * nn[0] + ( coord[0][1] - coord[0][2] ) * nn[1];
      m1 = ( coord[1][2] - coord[1][0] ) * nn[0] + ( coord[0][0] - coord[0][2] ) * nn[1];
      m2 = ( coord[1][0] - coord[1][1] ) * nn[0] + ( coord[0][1] - coord[0][0] ) * nn[1];

      btem[0][0] = ( -locdk[0]/3.0 + locdk[1]/30.0 - locdk[3]/5.0 ) * m0 * tem0;
      btem[0][1] = ( -locdk[0]/15.0 - locdk[1]/30.0 - locdk[3]/15.0 ) * m1 * tem0;
      btem[0][2] = ( -2.0*locdk[0]/15.0 + locdk[1]/30.0 - locdk[3]/15.0 ) * m2 * tem0;
      btem[0][3] = ( 2.0*locdk[0]/5.0 + 4.0*locdk[3]/15.0 ) * m1 * tem0
	+ ( -locdk[0]/15.0 + locdk[1]/15.0 ) * m2 * tem0;
      btem[0][4] = -( -locdk[0]/15.0 + locdk[1]/15.0 ) * m2 * tem0;
      btem[0][5] = ( 7.0*locdk[0]/15.0 - locdk[1]/15.0 + 4.0*locdk[3]/15.0 ) * m2 * tem0;

      btem[1][0] = ( -locdk[0]/5.0 + locdk[1]/15.0 - 8.0*locdk[3]/15.0 ) * m0 * tem0;
      btem[1][1] = ( -locdk[0]/15.0 + locdk[1]/5.0 + 8.0*locdk[3]/15.0 ) * m1 * tem0;
      btem[1][2] = ( -locdk[0]/15.0 - locdk[1]/15.0 - 8.0*locdk[3]/15.0 ) * m2 * tem0;
      btem[1][3] = ( 4.0*locdk[0]/15.0 - 4.0*locdk[1]/15.0 ) * m1 * tem0
	+ ( -4.0*locdk[1]/15.0 - 16.0*locdk[3]/15.0 ) * m2 * tem0;
      btem[1][4] = ( 4.0*locdk[1]/15.0 + 16.0*locdk[3]/15.0 ) * m2 * tem0;
      btem[1][5] = ( 4.0*locdk[0]/15.0 + 16.0*locdk[3]/15.0 ) * m2 * tem0;

      btem[2][0] = ( locdk[0]/30.0 + locdk[1]/15.0 + locdk[3]/15.0 ) * m0 * tem0;
      btem[2][1] = ( -locdk[0]/30.0 + locdk[1]/3.0 + locdk[3]/5.0 ) * m1 * tem0;
      btem[2][2] = ( locdk[0]/30.0 - 2.0*locdk[1]/15.0 - locdk[3]/15.0 ) * m2 * tem0;
      btem[2][3] = ( -2.0*locdk[1]/5.0 - 4.0*locdk[3]/15.0 ) * m1 * tem0
	+ ( locdk[0]/15.0 - 7.0*locdk[1]/15.0 - 4.0*locdk[3]/15.0 ) * m2 * tem0;
      btem[2][4] = ( -locdk[0]/15.0 + 7.0*locdk[1]/15.0 + 4.0*locdk[3]/15.0 ) * m2 * tem0;
      btem[2][5] = ( -locdk[0]/15.0 + locdk[1]/15.0 ) * m2 * tem0;

      //---1------------------
      xa = coord[0][1];
      xb = coord[0][2];
      ya = coord[1][1];
      yb = coord[1][2];
      len = sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya) );
      nn[0] = (yb - ya)/len;
      nn[1] = (xa - xb)/len;
      m0 = ( coord[1][2] - coord[1][1] ) * nn[0] + ( coord[0][1] - coord[0][2] ) * nn[1];
      m1 = ( coord[1][2] - coord[1][0] ) * nn[0] + ( coord[0][0] - coord[0][2] ) * nn[1];
      m2 = ( coord[1][0] - coord[1][1] ) * nn[0] + ( coord[0][1] - coord[0][0] ) * nn[1];

      btem[3][0] = ( 2.0*locdk[1]/15.0 - locdk[2]/30.0 + locdk[4]/15.0 ) * m0 * tem1;
      btem[3][1] = ( locdk[1]/3.0 - locdk[2]/30.0 + locdk[4]/5.0 ) * m1 * tem1;
      btem[3][2] = ( -locdk[1]/15.0 - locdk[2]/30.0 - locdk[4]/15.0 ) * m2 * tem1;
      btem[3][3] = ( -7.0*locdk[1]/15.0 + locdk[2]/15.0 - 4.0*locdk[4]/15.0 ) * (m1 + m2) * tem1;
      btem[3][4] = ( locdk[1]/15.0 - locdk[2]/15.0 ) * m1 * tem1
	+ ( 7.0*locdk[1]/15.0 - locdk[2]/15.0 + 4.0*locdk[4]/15.0 ) * m2 * tem1;
      btem[3][5] = ( -locdk[1]/15.0 + locdk[2]/15.0 ) * (m1 + m2) * tem1;

      btem[4][0] = ( locdk[1]/15.0 + locdk[2]/15.0 + 8.0*locdk[4]/15.0 ) * m0 * tem1;
      btem[4][1] = ( locdk[1]/5.0 - locdk[2]/15.0 + 8.0*locdk[4]/15.0 ) * m1 * tem1;
      btem[4][2] = ( -locdk[1]/15.0 + locdk[2]/5.0 + 8.0*locdk[4]/15.0 ) * m2 * tem1;
      btem[4][3] = ( -4.0*locdk[1]/15.0 - 16.0*locdk[4]/15.0 ) * (m1 + m2) * tem1;
      btem[4][4] = ( 4.0*locdk[2]/15.0 + 16.0*locdk[4]/15.0 ) * m1 * tem1
	+ ( 4.0*locdk[1]/15.0 + 16.0*locdk[4]/15.0 ) * m2 * tem1;
      btem[4][5] = ( -4.0*locdk[2]/15.0 - 16.0*locdk[4]/15.0 ) * (m1 + m2) * tem1;

      btem[5][0] = ( -locdk[1]/30.0 + 2.0*locdk[2]/15.0 + locdk[4]/15.0 ) * m0 * tem1;
      btem[5][1] = ( -locdk[1]/30.0 - locdk[2]/15.0 - locdk[4]/15.0 ) * m1 * tem1;
      btem[5][2] = ( -locdk[1]/30.0 + locdk[2]/3.0 + locdk[4]/5.0 ) * m2 * tem1;
      btem[5][3] = ( locdk[1]/15.0 - locdk[2]/15.0 ) * (m1 + m2) * tem1;
      btem[5][4] = ( -locdk[1]/15.0 + 7.0*locdk[2]/15.0 + 4.0*locdk[4]/15.0 ) * m1 * tem1
	+ ( -locdk[1]/15.0 + locdk[2]/15.0 ) * m2 * tem1;
      btem[5][5] = ( locdk[1]/15.0 - 7.0*locdk[2]/15.0 - 4.0*locdk[4]/15.0 ) * (m1 + m2) * tem1;

      //---2---------------
      xa = coord[0][2];
      xb = coord[0][0];
      ya = coord[1][2];
      yb = coord[1][0];
      len = sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya) );
      nn[0] = (yb - ya)/len;
      nn[1] = (xa - xb)/len;
      m0 = ( coord[1][2] - coord[1][1] ) * nn[0] + ( coord[0][1] - coord[0][2] ) * nn[1];
      m1 = ( coord[1][2] - coord[1][0] ) * nn[0] + ( coord[0][0] - coord[0][2] ) * nn[1];
      m2 = ( coord[1][0] - coord[1][1] ) * nn[0] + ( coord[0][1] - coord[0][0] ) * nn[1];

      btem[6][0] = ( locdk[0]/30.0 + locdk[2]/15.0 + locdk[5]/15.0 ) * m0 * tem2;
      btem[6][1] = ( locdk[0]/30.0 - 2.0*locdk[2]/15.0 - locdk[5]/15.0 ) * m1 * tem2;
      btem[6][2] = ( -locdk[0]/30.0 + locdk[2]/3.0 + locdk[5]/5.0 ) * m2 * tem2;
      btem[6][3] = ( -locdk[0]/15.0 + locdk[2]/15.0 ) * m1 * tem2;
      btem[6][4] = ( -locdk[0]/15.0 + 7.0*locdk[2]/15.0 + 4.0*locdk[5]/15.0 ) * m1 * tem2;
      btem[6][5] = -( -locdk[0]/15.0 + 7.0*locdk[2]/15.0 + 4.0*locdk[5]/15.0 ) * m1 * tem2
	+ ( -2.0*locdk[2]/5.0 - 4.0*locdk[5]/15.0 ) * m2 * tem2;

      btem[7][0] = ( -locdk[0]/5.0 + locdk[2]/15.0 - 8.0*locdk[5]/15.0 ) * m0 * tem2;
      btem[7][1] = ( -locdk[0]/15.0 - locdk[2]/15.0 - 8.0*locdk[5]/15.0 ) * m1 * tem2;
      btem[7][2] = ( -locdk[0]/15.0 + locdk[2]/5.0 + 8.0*locdk[5]/15.0 ) * m2 * tem2;
      btem[7][3] = ( 4.0*locdk[0]/15.0 + 16.0*locdk[5]/15.0 ) * m1 * tem2;
      btem[7][4] = ( 4.0*locdk[2]/15.0 + 16.0*locdk[5]/15.0 ) * m1 * tem2;
      btem[7][5] = -( 4.0*locdk[2]/15.0 + 16.0*locdk[5]/15.0 ) * m1 * tem2
	+ ( 4.0*locdk[0]/15.0 - 4.0*locdk[2]/15.0 ) * m2 * tem2;

      btem[8][0] = ( -locdk[0]/3.0 + locdk[2]/30.0 - locdk[5]/5.0 ) * m0 * tem2;
      btem[8][1] = ( -2.0*locdk[0]/15.0 + locdk[2]/30.0 - locdk[5]/15.0 ) * m1 * tem2;
      btem[8][2] = ( -locdk[0]/15.0 - locdk[2]/30.0 - locdk[5]/15.0 ) * m2 * tem2;
      btem[8][3] = ( 7.0*locdk[0]/15.0 - locdk[2]/15.0 + 4.0*locdk[5]/15.0 ) * m1 * tem2;
      btem[8][4] = ( locdk[0]/15.0 - locdk[2]/15.0 ) * m1 * tem2;
      btem[8][5] = ( -locdk[0]/15.0 + locdk[2]/15.0 ) * m1 * tem2
	+ ( 2.0*locdk[0]/5.0 + 4.0*locdk[5]/15.0 ) * m2 * tem2;

      for(int k=0; k<9; k++)
	for(int l=0; l<6; l++)
	  tflux[i][k] += btem[k][l] * locsol(l);	  
      
      for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
      for (int k=0; k<btem.Size(); k++) { delete []btem[k]; }
    }

  Vector cverr(dualmesh->GetDualMeshNumDOF());
  cverr = 0.0;
  
  Vector cverruh(dualmesh->GetDualMeshNumDOF());
  cverruh = 0.0;
  
  double bmaxuherr = 0.0;
  double bmaxpuherr = 0.0;
  double bl2uherr = 0.0;
  double bl2puherr = 0.0;

  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      fem->getElementDOF(i, ind);
      Array<int> edges;
      mesh->GetElementEdges(i, edges);
      
      Array<double *> coord(mesh->GetDim());
      for (int k=0; k<coord.Size(); k++) coord[k] = new double[3];
      
      mesh->GetElementVerticesCoord(i, coord);

      double t1 = coord[0][1] - coord[0][0];
      double t2 = coord[0][2] - coord[0][0];
      double t3 = coord[1][1] - coord[1][0];
      double t4 = coord[1][2] - coord[1][0];
      double detJ = t1 * t4 - t2 * t3;
      
      dualmesh->GetElemDualInfo(i, cod, nmls, lens, areas);

      Array<double> locdk(ind.Size());
      for (int j=0; j<ind.Size(); j++)
	locdk[j] = data->GetNodalEllipticCoeff(ind[j]);
      
      Array<double *> qtem(9);
      for(int k=0; k<9; k++) qtem[k] = new double[6];
      for(int k=0; k<9; k++)
	for(int j=0; j<6; j++)
	  qtem[k][j] = 0.0;

      double tem0 = sqrt( (coord[0][1] + coord[0][0] - 2.0*coord[0][2])*(coord[0][1] + coord[0][0] - 2.0*coord[0][2]) +
			  (coord[1][1] + coord[1][0] - 2.0*coord[1][2])*(coord[1][1] + coord[1][0] - 2.0*coord[1][2]) ) / detJ;
      
      double tem1 = sqrt( (coord[0][1] + coord[0][2] - 2.0*coord[0][0])*(coord[0][1] + coord[0][2] - 2.0*coord[0][0]) +
			  (coord[1][1] + coord[1][2] - 2.0*coord[1][0])*(coord[1][1] + coord[1][2] - 2.0*coord[1][0]) ) / detJ;
      
      double tem2 = sqrt( (coord[0][2] + coord[0][0] - 2.0*coord[0][1])*(coord[0][2] + coord[0][0] - 2.0*coord[0][1]) +
			  (coord[1][2] + coord[1][0] - 2.0*coord[1][1])*(coord[1][2] + coord[1][0] - 2.0*coord[1][1]) ) / detJ;

      double m0, m1, m2;

      //---0, 5, 7--------------------------
      m0 = ( coord[1][2] - coord[1][1] ) * nmls[0][0] + ( coord[0][1] - coord[0][2] ) * nmls[0][1];
      m1 = ( coord[1][2] - coord[1][0] ) * nmls[0][0] + ( coord[0][0] - coord[0][2] ) * nmls[0][1];
      m2 = ( coord[1][0] - coord[1][1] ) * nmls[0][0] + ( coord[0][1] - coord[0][0] ) * nmls[0][1];
      
      qtem[0][0] = ( -473.0*locdk[0]/10368.0 + 191.0*locdk[1]/10368.0 + 25.0*locdk[2]/2592.0 
		     - 473.0*locdk[3]/5184.0 - 25.0*locdk[4]/2592.0 - 89.0*locdk[5]/2592.0 ) * m0 * tem0;
      qtem[0][1] = ( -13.0*locdk[0]/3456.0 + 17.0*locdk[1]/10368.0 + 1.0*locdk[2]/864.0 
		     - 13.0*locdk[3]/1728.0 - 1.0*locdk[4]/864.0 - 11.0*locdk[5]/2592.0 ) * m1 * tem0;
      qtem[0][2] = ( -89.0*locdk[0]/5184.0 + 35.0*locdk[1]/5184.0 + 1.0*locdk[2]/324.0 
		     - 89.0*locdk[3]/2592.0 - 1.0*locdk[4]/324.0 - 7.0*locdk[5]/648.0 ) * m2 * tem0;

      qtem[0][3] = ( 4.0*locdk[0]/81.0 - 13.0*locdk[1]/648.0 - 7.0*locdk[2]/648.0 
		     + 8.0*locdk[3]/81.0 + 7.0*locdk[4]/648.0 + 25.0*locdk[5]/648.0 ) * m1 * tem0
	+ ( -217.0*locdk[0]/10368.0 + 29.0*locdk[1]/3456.0 + 11.0*locdk[2]/2592.0 
	    - 217.0*locdk[3]/5184.0 - 11.0*locdk[4]/2592.0 - 13.0*locdk[5]/864.0 ) * m2 * tem0;

      qtem[0][4] = ( 13.0*locdk[0]/1728.0 - 17.0*locdk[1]/5184.0 - 1.0*locdk[2]/432.0 
		     + 13.0*locdk[3]/864.0 + 1.0*locdk[4]/432.0 + 11.0*locdk[5]/1296.0 ) * m1 * tem0
	- ( -217.0*locdk[0]/10368.0 + 29.0*locdk[1]/3456.0 + 11.0*locdk[2]/2592.0 
	    - 217.0*locdk[3]/5184.0 - 11.0*locdk[4]/2592.0 - 13.0*locdk[5]/864.0 ) * m2 * tem0;

      qtem[0][5] = -( 13.0*locdk[0]/1728.0 - 17.0*locdk[1]/5184.0 - 1.0*locdk[2]/432.0 
		     + 13.0*locdk[3]/864.0 + 1.0*locdk[4]/432.0 + 11.0*locdk[5]/1296.0 ) * m1 * tem0
	+ ( 217.0*locdk[0]/3456.0 - 29.0*locdk[1]/1152.0 - 11.0*locdk[2]/864.0 
	    + 217.0*locdk[3]/1728.0 + 11.0*locdk[4]/864.0 + 13.0*locdk[5]/288.0 ) * m2 * tem0;

      
      qtem[5][0] = ( -17.0*locdk[0]/10368.0 + 13.0*locdk[1]/3456.0 - 1.0*locdk[2]/864.0 
		     + 13.0*locdk[3]/1728.0 + 11.0*locdk[4]/2592.0 + 1.0*locdk[5]/864.0 ) * m0 * tem0;
      qtem[5][1] = ( -191.0*locdk[0]/10368.0 + 473.0*locdk[1]/10368.0 - 25.0*locdk[2]/2592.0 
		     + 473.0*locdk[3]/5184.0 + 89.0*locdk[4]/2592.0 + 25.0*locdk[5]/2592.0 ) * m1 * tem0;
      qtem[5][2] = ( 35.0*locdk[0]/5184.0 - 89.0*locdk[1]/5184.0 + 1.0*locdk[2]/324.0 
		     - 89.0*locdk[3]/2592.0 - 7.0*locdk[4]/648.0 - 1.0*locdk[5]/324.0 ) * m2 * tem0;

      qtem[5][3] = ( 13.0*locdk[0]/648.0 - 4.0*locdk[1]/81.0 + 7.0*locdk[2]/648.0 
		     - 8.0*locdk[3]/81.0 - 25.0*locdk[4]/648.0 - 7.0*locdk[5]/648.0 ) * m1 * tem0
	- ( -295.0*locdk[0]/10368.0 + 9.0*locdk[1]/128.0 - 13.0*locdk[2]/864.0 
	    + 9.0*locdk[3]/64.0 + 139.0*locdk[4]/2592.0 + 13.0*locdk[5]/864.0 ) * m2 * tem0;

      qtem[5][4] = ( -17.0*locdk[0]/5184.0 + 13.0*locdk[1]/1728.0 - 1.0*locdk[2]/432.0 
		     + 13.0*locdk[3]/864.0 + 11.0*locdk[4]/1296.0 + 1.0*locdk[5]/432.0 ) * m1 * tem0
	+ ( -295.0*locdk[0]/10368.0 + 9.0*locdk[1]/128.0 - 13.0*locdk[2]/864.0 
	    + 9.0*locdk[3]/64.0 + 139.0*locdk[4]/2592.0 + 13.0*locdk[5]/864.0 ) * m2 * tem0;

      qtem[5][5] = -( -17.0*locdk[0]/5184.0 + 13.0*locdk[1]/1728.0 - 1.0*locdk[2]/432.0 
		      + 13.0*locdk[3]/864.0 + 11.0*locdk[4]/1296.0 + 1.0*locdk[5]/432.0 ) * m1 * tem0
	+ ( -53.0*locdk[0]/10368.0 + 139.0*locdk[1]/10368.0 - 5.0*locdk[2]/2592.0 
	    + 139.0*locdk[3]/5184.0 + 17.0*locdk[4]/2592.0 + 5.0*locdk[5]/2592.0 ) * m2 * tem0;

      
      qtem[7][0] = ( 1.0*locdk[2]/324.0 - 1.0*locdk[3]/324.0 ) * m0 * tem0;
      qtem[7][1] = ( -1.0*locdk[2]/324.0 + 1.0*locdk[3]/324.0 ) * m1 * tem0;
      qtem[7][2] = ( -13.0*locdk[0]/648.0 - 13.0*locdk[1]/648.0 + 1.0*locdk[2]/108.0 
		     + 1.0*locdk[3]/27.0 + 13.0*locdk[4]/162.0 + 13.0*locdk[5]/162.0 ) * m2 * tem0;

      qtem[7][3] = -( -13.0*locdk[0]/648.0 - 13.0*locdk[1]/648.0
		      + 5.0*locdk[3]/108.0 + 13.0*locdk[4]/162.0 + 13.0*locdk[5]/162.0 ) * m2 * tem0;

      qtem[7][4] = ( -13.0*locdk[0]/324.0 - 13.0*locdk[1]/324.0 + 1.0*locdk[2]/81.0 
		     + 13.0*locdk[3]/162.0 + 13.0*locdk[4]/81.0 + 13.0*locdk[5]/81.0 ) * m1 * tem0
	+ ( -13.0*locdk[0]/648.0 - 13.0*locdk[1]/648.0
	    + 5.0*locdk[3]/108.0 + 13.0*locdk[4]/162.0 + 13.0*locdk[5]/162.0 ) * m2 * tem0;
     
      qtem[7][5] = -( -13.0*locdk[0]/324.0 - 13.0*locdk[1]/324.0 + 1.0*locdk[2]/81.0 
		      + 13.0*locdk[3]/162.0 + 13.0*locdk[4]/81.0 + 13.0*locdk[5]/81.0 ) * m1 * tem0
	- ( -13.0*locdk[0]/648.0 - 13.0*locdk[1]/648.0 + 1.0*locdk[2]/81.0
	    + 11.0*locdk[3]/324.0 + 13.0*locdk[4]/162.0 + 13.0*locdk[5]/162.0 ) * m2 * tem0;
     
      //---1, 3, 8--------------------------
      m0 = ( coord[1][2] - coord[1][1] ) * nmls[1][0] + ( coord[0][1] - coord[0][2] ) * nmls[1][1];
      m1 = ( coord[1][2] - coord[1][0] ) * nmls[1][0] + ( coord[0][0] - coord[0][2] ) * nmls[1][1];
      m2 = ( coord[1][0] - coord[1][1] ) * nmls[1][0] + ( coord[0][1] - coord[0][0] ) * nmls[1][1];
     
      qtem[1][0] = ( -1.0*locdk[0]/108.0 + 13.0*locdk[1]/648.0 + 13.0*locdk[2]/648.0 
		     - 13.0*locdk[3]/162.0 - 1.0*locdk[4]/27.0 - 13.0*locdk[5]/162.0 ) * m0 * tem1;
      qtem[1][1] = ( -1.0*locdk[0]/324.0 + 1.0*locdk[4]/324.0 ) * m1 * tem1;
      qtem[1][2] = ( -1.0*locdk[0]/324.0 + 1.0*locdk[4]/324.0 ) * m2 * tem1;
     
      qtem[1][3] = ( 1.0*locdk[0]/81.0 - 13.0*locdk[1]/648.0 - 13.0*locdk[2]/648.0
		     + 13.0*locdk[3]/162.0 + 11.0*locdk[4]/324.0 + 13.0*locdk[5]/162.0 ) * m1 * tem1
	- ( -13.0*locdk[1]/648.0 - 13.0*locdk[2]/648.0 + 13.0*locdk[3]/162.0 
	    + 5.0*locdk[4]/108.0 + 13.0*locdk[5]/162.0 ) * m2 * tem1;
      
      qtem[1][4] = ( -13.0*locdk[1]/648.0 - 13.0*locdk[2]/648.0 + 13.0*locdk[3]/162.0 
		     + 5.0*locdk[4]/108.0 + 13.0*locdk[5]/162.0 ) * m1 * tem1
	+ ( -13.0*locdk[1]/648.0 - 13.0*locdk[2]/648.0 + 13.0*locdk[3]/162.0 
	    + 5.0*locdk[4]/108.0 + 13.0*locdk[5]/162.0 ) * m2 * tem1;
      
      qtem[1][5] = -( -13.0*locdk[1]/648.0 - 13.0*locdk[2]/648.0 + 13.0*locdk[3]/162.0 
		     + 5.0*locdk[4]/108.0 + 13.0*locdk[5]/162.0 ) * m1 * tem1
	+ ( 1.0*locdk[0]/81.0 - 13.0*locdk[1]/648.0 - 13.0*locdk[2]/648.0
	    + 13.0*locdk[3]/162.0 + 11.0*locdk[4]/324.0 + 13.0*locdk[5]/162.0 ) * m2 * tem1;
      
     
      qtem[3][0] = ( -1.0*locdk[0]/324.0 + 89.0*locdk[1]/5184.0 - 35.0*locdk[2]/5184.0 
		     + 7.0*locdk[3]/648.0 + 89.0*locdk[4]/2592.0 + 1.0*locdk[5]/324.0 ) * m0 * tem1;
      qtem[3][1] = ( -25.0*locdk[0]/2592.0 + 473.0*locdk[1]/10368.0 - 191.0*locdk[2]/10368.0 
		     + 89.0*locdk[3]/2592.0 + 473.0*locdk[4]/5184.0 + 25.0*locdk[5]/2592.0 ) * m1 * tem1;
      qtem[3][2] = ( 1.0*locdk[0]/864.0 - 13.0*locdk[1]/3456.0 + 17.0*locdk[2]/10368.0 
		     - 11.0*locdk[3]/2592.0 - 13.0*locdk[4]/1728.0 - 1.0*locdk[5]/864.0 ) * m2 * tem1;

      qtem[3][3] = ( 11.0*locdk[0]/864.0 - 217.0*locdk[1]/3456.0 + 29.0*locdk[2]/1152.0 
		     - 13.0*locdk[3]/288.0 - 217.0*locdk[4]/1728.0 - 11.0*locdk[5]/864.0 ) * m1 * tem1
	+ ( 13.0*locdk[0]/864.0 - 9.0*locdk[1]/128.0 + 295.0*locdk[2]/10368.0 
	    - 139.0*locdk[3]/2592.0 - 9.0*locdk[4]/64.0 - 13.0*locdk[5]/864.0 ) * m2 * tem1;

      qtem[3][4] = -( 11.0*locdk[0]/2592.0 - 217.0*locdk[1]/10368.0 + 29.0*locdk[2]/3456.0 
		      - 13.0*locdk[3]/864.0 - 217.0*locdk[4]/5184.0 - 11.0*locdk[5]/2592.0 ) * m1 * tem1
	- ( 13.0*locdk[0]/864.0 - 9.0*locdk[1]/128.0 + 295.0*locdk[2]/10368.0 
	    - 139.0*locdk[3]/2592.0 - 9.0*locdk[4]/64.0 - 13.0*locdk[5]/864.0 ) * m2 * tem1;

      qtem[3][5] = ( 11.0*locdk[0]/2592.0 - 217.0*locdk[1]/10368.0 + 29.0*locdk[2]/3456.0 
		     - 13.0*locdk[3]/864.0 - 217.0*locdk[4]/5184.0 - 11.0*locdk[5]/2592.0 ) * m1 * tem1
	+ ( 5.0*locdk[0]/2592.0 - 139.0*locdk[1]/10368.0 + 53.0*locdk[2]/10368.0 
	    - 17.0*locdk[3]/2592.0 - 139.0*locdk[4]/5184.0 - 5.0*locdk[5]/2592.0 ) * m2 * tem1;

      
      qtem[8][0] = ( -1.0*locdk[0]/324.0 + 89.0*locdk[2]/5184.0 - 35.0*locdk[1]/5184.0 
		     + 7.0*locdk[5]/648.0 + 89.0*locdk[4]/2592.0 + 1.0*locdk[3]/324.0 ) * m0 * tem1;
      qtem[8][1] = ( 1.0*locdk[0]/864.0 - 13.0*locdk[2]/3456.0 + 17.0*locdk[1]/10368.0 
		     - 11.0*locdk[5]/2592.0 - 13.0*locdk[4]/1728.0 - 1.0*locdk[3]/864.0 ) * m1 * tem1;
      qtem[8][2] = ( -25.0*locdk[0]/2592.0 + 473.0*locdk[2]/10368.0 - 191.0*locdk[1]/10368.0 
		     + 89.0*locdk[5]/2592.0 + 473.0*locdk[4]/5184.0 + 25.0*locdk[3]/2592.0 ) * m2 * tem1;

      qtem[8][3] = ( 5.0*locdk[0]/2592.0 - 139.0*locdk[2]/10368.0 + 53.0*locdk[1]/10368.0 
	    - 17.0*locdk[5]/2592.0 - 139.0*locdk[4]/5184.0 - 5.0*locdk[3]/2592.0 ) * m1 * tem1
	+ ( 11.0*locdk[0]/2592.0 - 217.0*locdk[2]/10368.0 + 29.0*locdk[1]/3456.0 
	   - 13.0*locdk[5]/864.0 - 217.0*locdk[4]/5184.0 - 11.0*locdk[3]/2592.0 ) * m2 * tem1;

      qtem[8][4] = -( 13.0*locdk[0]/864.0 - 9.0*locdk[2]/128.0 + 295.0*locdk[1]/10368.0 
		       - 139.0*locdk[5]/2592.0 - 9.0*locdk[4]/64.0 - 13.0*locdk[3]/864.0 ) * m1 * tem1
	- ( 11.0*locdk[0]/2592.0 - 217.0*locdk[2]/10368.0 + 29.0*locdk[1]/3456.0 
	   - 13.0*locdk[5]/864.0 - 217.0*locdk[4]/5184.0 - 11.0*locdk[3]/2592.0 ) * m2 * tem1;
	
      qtem[8][5] = ( 13.0*locdk[0]/864.0 - 9.0*locdk[2]/128.0 + 295.0*locdk[1]/10368.0 
		     - 139.0*locdk[5]/2592.0 - 9.0*locdk[4]/64.0 - 13.0*locdk[3]/864.0 ) * m1 * tem1
	+ ( 11.0*locdk[0]/864.0 - 217.0*locdk[2]/3456.0 + 29.0*locdk[1]/1152.0 
	    - 13.0*locdk[5]/288.0 - 217.0*locdk[4]/1728.0 - 11.0*locdk[3]/864.0 ) * m2 * tem1;

      
      //---2, 4, 6--------------------------
      m0 = ( coord[1][2] - coord[1][1] ) * nmls[2][0] + ( coord[0][1] - coord[0][2] ) * nmls[2][1];
      m1 = ( coord[1][2] - coord[1][0] ) * nmls[2][0] + ( coord[0][0] - coord[0][2] ) * nmls[2][1];
      m2 = ( coord[1][0] - coord[1][1] ) * nmls[2][0] + ( coord[0][1] - coord[0][0] ) * nmls[2][1];
      
      qtem[2][0] = ( -473.0*locdk[0]/10368.0 + 191.0*locdk[2]/10368.0 + 25.0*locdk[1]/2592.0 
		     - 473.0*locdk[5]/5184.0 - 25.0*locdk[4]/2592.0 - 89.0*locdk[3]/2592.0 ) * m0 * tem2;
      qtem[2][1] = ( -89.0*locdk[0]/5184.0 + 35.0*locdk[2]/5184.0 + 1.0*locdk[1]/324.0 
		     - 89.0*locdk[5]/2592.0 - 1.0*locdk[4]/324.0 - 7.0*locdk[3]/648.0 ) * m1 * tem2;
      qtem[2][2] = ( -13.0*locdk[0]/3456.0 + 17.0*locdk[2]/10368.0 + 1.0*locdk[1]/864.0 
		     - 13.0*locdk[5]/1728.0 - 1.0*locdk[4]/864.0 - 11.0*locdk[3]/2592.0 ) * m2 * tem2;

      qtem[2][3] = -( 13.0*locdk[0]/1728.0 - 17.0*locdk[2]/5184.0 - 1.0*locdk[1]/432.0 
		     + 13.0*locdk[5]/864.0 + 1.0*locdk[4]/432.0 + 11.0*locdk[3]/1296.0 ) * m2 * tem2
	+ ( 217.0*locdk[0]/3456.0 - 29.0*locdk[2]/1152.0 - 11.0*locdk[1]/864.0 
	    + 217.0*locdk[5]/1728.0 + 11.0*locdk[4]/864.0 + 13.0*locdk[3]/288.0 ) * m1 * tem2;

      qtem[2][4] = ( 13.0*locdk[0]/1728.0 - 17.0*locdk[2]/5184.0 - 1.0*locdk[1]/432.0 
		     + 13.0*locdk[5]/864.0 + 1.0*locdk[4]/432.0 + 11.0*locdk[3]/1296.0 ) * m2 * tem2
	- ( -217.0*locdk[0]/10368.0 + 29.0*locdk[2]/3456.0 + 11.0*locdk[1]/2592.0 
	    - 217.0*locdk[5]/5184.0 - 11.0*locdk[4]/2592.0 - 13.0*locdk[3]/864.0 ) * m1 * tem2;
 
      qtem[2][5] = ( 4.0*locdk[0]/81.0 - 13.0*locdk[2]/648.0 - 7.0*locdk[1]/648.0 
		     + 8.0*locdk[5]/81.0 + 7.0*locdk[4]/648.0 + 25.0*locdk[3]/648.0 ) * m2 * tem2
	+ ( -217.0*locdk[0]/10368.0 + 29.0*locdk[2]/3456.0 + 11.0*locdk[1]/2592.0 
	    - 217.0*locdk[5]/5184.0 - 11.0*locdk[4]/2592.0 - 13.0*locdk[3]/864.0 ) * m1 * tem2;

     
      qtem[4][0] = ( 1.0*locdk[1]/324.0 - 1.0*locdk[5]/324.0 ) * m0 * tem2;
      qtem[4][1] = ( -13.0*locdk[0]/648.0 - 13.0*locdk[2]/648.0 + 1.0*locdk[1]/108.0 
		     + 1.0*locdk[5]/27.0 + 13.0*locdk[4]/162.0 + 13.0*locdk[3]/162.0 ) * m1 * tem2;
      qtem[4][2] = ( -1.0*locdk[1]/324.0 + 1.0*locdk[5]/324.0 ) * m2 * tem2;

      qtem[4][3] = -( -13.0*locdk[0]/324.0 - 13.0*locdk[2]/324.0 + 1.0*locdk[1]/81.0 
		      + 13.0*locdk[5]/162.0 + 13.0*locdk[4]/81.0 + 13.0*locdk[3]/81.0 ) * m2 * tem2
	- ( -13.0*locdk[0]/648.0 - 13.0*locdk[2]/648.0 + 1.0*locdk[1]/81.0
	    + 11.0*locdk[5]/324.0 + 13.0*locdk[4]/162.0 + 13.0*locdk[3]/162.0 ) * m1 * tem2;

      qtem[4][4] = ( -13.0*locdk[0]/324.0 - 13.0*locdk[2]/324.0 + 1.0*locdk[1]/81.0 
		     + 13.0*locdk[5]/162.0 + 13.0*locdk[4]/81.0 + 13.0*locdk[3]/81.0 ) * m2 * tem2
	+ ( -13.0*locdk[0]/648.0 - 13.0*locdk[2]/648.0
	    + 5.0*locdk[5]/108.0 + 13.0*locdk[4]/162.0 + 13.0*locdk[3]/162.0 ) * m1 * tem2;
     
      qtem[4][5] = -( -13.0*locdk[0]/648.0 - 13.0*locdk[2]/648.0
		      + 5.0*locdk[5]/108.0 + 13.0*locdk[4]/162.0 + 13.0*locdk[3]/162.0 ) * m1 * tem2;

      
      qtem[6][0] = ( -17.0*locdk[0]/10368.0 + 13.0*locdk[2]/3456.0 - 1.0*locdk[1]/864.0 
		     + 13.0*locdk[5]/1728.0 + 11.0*locdk[4]/2592.0 + 1.0*locdk[3]/864.0 ) * m0 * tem2;
      qtem[6][1] = ( 35.0*locdk[0]/5184.0 - 89.0*locdk[2]/5184.0 + 1.0*locdk[1]/324.0 
		     - 89.0*locdk[5]/2592.0 - 7.0*locdk[4]/648.0 - 1.0*locdk[3]/324.0 ) * m1 * tem2;
      qtem[6][2] = ( -191.0*locdk[0]/10368.0 + 473.0*locdk[2]/10368.0 - 25.0*locdk[1]/2592.0 
		     + 473.0*locdk[5]/5184.0 + 89.0*locdk[4]/2592.0 + 25.0*locdk[3]/2592.0 ) * m2 * tem2;

      qtem[6][3] = -( -17.0*locdk[0]/5184.0 + 13.0*locdk[2]/1728.0 - 1.0*locdk[1]/432.0 
		      + 13.0*locdk[5]/864.0 + 11.0*locdk[4]/1296.0 + 1.0*locdk[3]/432.0 ) * m2 * tem2
	+ ( -53.0*locdk[0]/10368.0 + 139.0*locdk[2]/10368.0 - 5.0*locdk[1]/2592.0 
	    + 139.0*locdk[5]/5184.0 + 17.0*locdk[4]/2592.0 + 5.0*locdk[3]/2592.0 ) * m1 * tem2;

      qtem[6][4] = ( -17.0*locdk[0]/5184.0 + 13.0*locdk[2]/1728.0 - 1.0*locdk[1]/432.0 
		     + 13.0*locdk[5]/864.0 + 11.0*locdk[4]/1296.0 + 1.0*locdk[3]/432.0 ) * m2 * tem2
	+ ( -295.0*locdk[0]/10368.0 + 9.0*locdk[2]/128.0 - 13.0*locdk[1]/864.0 
	    + 9.0*locdk[5]/64.0 + 139.0*locdk[4]/2592.0 + 13.0*locdk[3]/864.0 ) * m1 * tem2;
      
      qtem[6][5] = ( 13.0*locdk[0]/648.0 - 4.0*locdk[2]/81.0 + 7.0*locdk[1]/648.0 
		     - 8.0*locdk[5]/81.0 - 25.0*locdk[4]/648.0 - 7.0*locdk[3]/648.0 ) * m2 * tem2
	- ( -295.0*locdk[0]/10368.0 + 9.0*locdk[2]/128.0 - 13.0*locdk[1]/864.0 
	    + 9.0*locdk[5]/64.0 + 139.0*locdk[4]/2592.0 + 13.0*locdk[3]/864.0 ) * m1 * tem2;


      SparseMatrix *locmat = new SparseMatrix(6,6);

      for(int j=0; j<6; j++)
	locmat->Elem(0,j) = qtem[2][j] - qtem[0][j];
      for(int j=0; j<6; j++)
	locmat->Elem(1,j) = qtem[5][j] - qtem[3][j];
      for(int j=0; j<6; j++)
	locmat->Elem(2,j) = qtem[8][j] - qtem[6][j];

      for(int j=0; j<6; j++)
	locmat->Elem(3,j) = qtem[0][j] - qtem[1][j] + qtem[4][j] - qtem[5][j];
      for(int j=0; j<6; j++)
	locmat->Elem(4,j) = qtem[3][j] - qtem[4][j] + qtem[7][j] - qtem[8][j];
      for(int j=0; j<6; j++)
	locmat->Elem(5,j) = qtem[6][j] - qtem[7][j] + qtem[1][j] - qtem[2][j];
            
      locmat->Finalize();
      //locmat->Print();
	
      Vector Q(6);
      Q = 0.0;
      Array<double *> femat(6);
      Vector locsol(6);
      for(int k=0; k<6; k++)
	{
	  femat[k] = new double[6];
	  locsol(k) = sol(ind[k]);
	}

      fem->GetEllipticLocalSystem(i, ind, femat);
 
      for(int k=0; k<6; k++)
	for(int j=0; j<6; j++)
	  Q(k) += femat[k][j] * locsol(j);
	  
      Array<double> F(6);
      fem->GetForceLocalSystem(i, ind, F);

      Vector f(6);
      f = 0.0;
      Function *ff = data->GetForceFunction();
      Array<double> x(3);
      Array<double> y(3);
      double gpx, gpy, xx, yy;

      for(int l=0; l<6; l++)
	{
	  if(l==0)
	    {
	      x[0] = coord[0][0];
	      y[0] = coord[1][0];
	      x[1] = cod[0][0];
	      y[1] = cod[0][1];
	      x[2] = cod[2][0];
	      y[2] = cod[2][1];
	    }
	  else if(l==1)
	    {
	      x[0] = cod[3][0];
	      y[0] = cod[3][1];
	      x[1] = cod[2][0];
	      y[1] = cod[2][1];
	      x[2] = cod[0][0];
	      y[2] = cod[0][1];
	    }
	  else if(l==2)
	    {
	      x[0] = coord[0][1];
	      y[0] = coord[1][1];
	      x[1] = cod[4][0];
	      y[1] = cod[4][1];
	      x[2] = cod[6][0];
	      y[2] = cod[6][1];
	    }
	  else if(l==3)
	    {
	      x[0] = cod[7][0];
	      y[0] = cod[7][1];
	      x[1] = cod[6][0];
	      y[1] = cod[6][1];
	      x[2] = cod[4][0];
	      y[2] = cod[4][1];
	    }
	  else if(l==4)
	    {
	      x[0] = coord[0][2];
	      y[0] = coord[1][2];
	      x[1] = cod[8][0];
	      y[1] = cod[8][1];
	      x[2] = cod[10][0];
	      y[2] = cod[10][1];
	    }
	  else
	    {
	      x[0] = cod[11][0];
	      y[0] = cod[11][1];
	      x[1] = cod[10][0];
	      y[1] = cod[10][1];
	      x[2] = cod[8][0];
	      y[2] = cod[8][1];
	    }
	  
	  t1 = x[1] - x[0];
	  t2 = x[2] - x[0];
	  t3 = y[1] - y[0];
	  t4 = y[2] - y[0];
	  detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<np; k++)
	    {
	      gpx = tpw[k][0];
	      gpy = tpw[k][1];
	      xx = x[0] + t1*gpx + t2*gpy;
	      yy = y[0] + t3*gpx + t4*gpy;
	      
	      double bcoord[2];
	      bcoord[0] = xx;
	      bcoord[1] = yy;
	      f(l/2) += tpw[k][2] * ff->Eval(bcoord) * detJ;
	    } 
	}

      for(int l=0; l<18; l++)
	{
	  if(l==0)
	    {
	      x[0] = 0.5 * (coord[0][0] + coord[0][1]);
	      y[0] = 0.5 * (coord[1][0] + coord[1][1]);
	      x[1] = cod[1][0];
	      y[1] = cod[1][1];
	      x[2] = cod[0][0];
	      y[2] = cod[0][1];
	    }
	  else if(l==1)
	    {
	      x[0] = cod[3][0];
	      y[0] = cod[3][1];
	      x[1] = cod[0][0];
	      y[1] = cod[0][1];
	      x[2] = cod[1][0];
	      y[2] = cod[1][1];
	    }
	  else if(l==2)
	    {
	      x[0] = 0.5 * (coord[0][0] + coord[0][1]);
	      y[0] = 0.5 * (coord[1][0] + coord[1][1]);
	      x[1] = cod[6][0];
	      y[1] = cod[6][1];
	      x[2] = cod[5][0];
	      y[2] = cod[5][1];
	    }
	  else if(l==3)
	    {
	      x[0] = cod[7][0];
	      y[0] = cod[7][1];
	      x[1] = cod[5][0];
	      y[1] = cod[5][1];
	      x[2] = cod[6][0];
	      y[2] = cod[6][1];
	    }
	  else if(l==4)
	    {
	      x[0] = 0.5 * (coord[0][0] + coord[0][1]);
	      y[0] = 0.5 * (coord[1][0] + coord[1][1]);
	      x[1] = cod[5][0];
	      y[1] = cod[5][1];
	      x[2] = cod[1][0];
	      y[2] = cod[1][1];
	    }
	  else if(l==5)
	    {
	      x[0] = cod[12][0];
	      y[0] = cod[12][1];
	      x[1] = cod[1][0];
	      y[1] = cod[1][1];
	      x[2] = cod[5][0];
	      y[2] = cod[5][1];
	    }
	  else if(l==6)
	    {
	      x[0] = 0.5 * (coord[0][2] + coord[0][1]);
	      y[0] = 0.5 * (coord[1][2] + coord[1][1]);
	      x[1] = cod[5][0];
	      y[1] = cod[5][1];
	      x[2] = cod[4][0];
	      y[2] = cod[4][1];
	    }
	  else if(l==7)
	    {
	      x[0] = cod[7][0];
	      y[0] = cod[7][1];
	      x[1] = cod[4][0];
	      y[1] = cod[4][1];
	      x[2] = cod[5][0];
	      y[2] = cod[5][1];
	    }
	  else if(l==8)
	    {
	      x[0] = 0.5 * (coord[0][2] + coord[0][1]);
	      y[0] = 0.5 * (coord[1][2] + coord[1][1]);
	      x[1] = cod[10][0];
	      y[1] = cod[10][1];
	      x[2] = cod[9][0];
	      y[2] = cod[9][1];
	    }
	  else if(l==9)
	    {
	      x[0] = cod[11][0];
	      y[0] = cod[11][1];
	      x[1] = cod[9][0];
	      y[1] = cod[9][1];
	      x[2] = cod[10][0];
	      y[2] = cod[10][1];
	    }
	  else if(l==10)
	    {
	      x[0] = 0.5 * (coord[0][2] + coord[0][1]);
	      y[0] = 0.5 * (coord[1][2] + coord[1][1]);
	      x[1] = cod[9][0];
	      y[1] = cod[9][1];
	      x[2] = cod[5][0];
	      y[2] = cod[5][1];
	    }
	  else if(l==11)
	    {
	      x[0] = cod[12][0];
	      y[0] = cod[12][1];
	      x[1] = cod[5][0];
	      y[1] = cod[5][1];
	      x[2] = cod[9][0];
	      y[2] = cod[9][1];
	    }
	  else if(l==12)
	    {
	      x[0] = 0.5 * (coord[0][2] + coord[0][0]);
	      y[0] = 0.5 * (coord[1][2] + coord[1][0]);
	      x[1] = cod[2][0];
	      y[1] = cod[2][1];
	      x[2] = cod[1][0];
	      y[2] = cod[1][1];
	    }
	  else if(l==13)
	    {
	      x[0] = cod[3][0];
	      y[0] = cod[3][1];
	      x[1] = cod[1][0];
	      y[1] = cod[1][1];
	      x[2] = cod[2][0];
	      y[2] = cod[2][1];

	    }
	  else if(l==14)
	    {
	      x[0] = 0.5 * (coord[0][2] + coord[0][0]);
	      y[0] = 0.5 * (coord[1][2] + coord[1][0]);
	      x[1] = cod[9][0];
	      y[1] = cod[9][1];
	      x[2] = cod[8][0];
	      y[2] = cod[8][1];
	    }
	  else if(l==15)
	    {
	      x[0] = cod[11][0];
	      y[0] = cod[11][1];
	      x[1] = cod[8][0];
	      y[1] = cod[8][1];
	      x[2] = cod[9][0];
	      y[2] = cod[9][1];
	    }
	  else if(l==16)
	    {
	      x[0] = 0.5 * (coord[0][2] + coord[0][0]);
	      y[0] = 0.5 * (coord[1][2] + coord[1][0]);
	      x[1] = cod[1][0];
	      y[1] = cod[1][1];
	      x[2] = cod[9][0];
	      y[2] = cod[9][1];
	    }
	  else
	    {
	      x[0] = cod[12][0];
	      y[0] = cod[12][1];
	      x[1] = cod[9][0];
	      y[1] = cod[9][1];
	      x[2] = cod[1][0];
	      y[2] = cod[1][1];
	    }

	  t1 = x[1] - x[0];
	  t2 = x[2] - x[0];
	  t3 = y[1] - y[0];
	  t4 = y[2] - y[0];
	  detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<np; k++)
	    {
	      gpx = tpw[k][0];
	      gpy = tpw[k][1];
	      xx = x[0] + t1*gpx + t2*gpy;
	      yy = y[0] + t3*gpx + t4*gpy;
	      
	      double bcoord[2];
	      bcoord[0] = xx;
	      bcoord[1] = yy;
	      f(3+l/6) += tpw[k][2] * ff->Eval(bcoord) * detJ;
	    } 
	}

      Vector RHS(6);
      for(int k=0; k<6; k++) RHS(k) = Q(k) - F[k] + f(k);

      Array<int> eind;
      if(i%2==0)
	{
	  mesh->GetEdgeElements(edges[0], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[1]][5] - tflux[i][0] + tflux[eind[1]][5] );
		  RHS(3) += 0.5 * ( sflux[i][1] - sflux[eind[1]][4] - tflux[i][1] + tflux[eind[1]][4] ); 
		  RHS(1) += 0.5 * ( sflux[i][2] - sflux[eind[1]][3] - tflux[i][2] + tflux[eind[1]][3] ); 
		}
	      else
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[0]][5] - tflux[i][0] + tflux[eind[0]][5] );
		  RHS(3) += 0.5 * ( sflux[i][1] - sflux[eind[0]][4] - tflux[i][1] + tflux[eind[0]][4] ); 
		  RHS(1) += 0.5 * ( sflux[i][2] - sflux[eind[0]][3] - tflux[i][2] + tflux[eind[0]][3] ); 
		}
	    }
	  else
	    {
	      RHS(0) += -tflux[i][0] + sflux[i][0];
	      RHS(3) += -tflux[i][1] + sflux[i][1];
	      RHS(1) += -tflux[i][2] + sflux[i][2];
	    }

	  mesh->GetEdgeElements(edges[1], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(1) += 0.5 * ( sflux[i][3] - sflux[eind[1]][8] - tflux[i][3] + tflux[eind[1]][8] );
		  RHS(4) += 0.5 * ( sflux[i][4] - sflux[eind[1]][7] - tflux[i][4] + tflux[eind[1]][7] ); 
		  RHS(2) += 0.5 * ( sflux[i][5] - sflux[eind[1]][6] - tflux[i][5] + tflux[eind[1]][6] ); 
		}
	      else
		{
		  RHS(1) += 0.5 * ( sflux[i][3] - sflux[eind[0]][8] - tflux[i][3] + tflux[eind[0]][8] );
		  RHS(4) += 0.5 * ( sflux[i][4] - sflux[eind[0]][7] - tflux[i][4] + tflux[eind[0]][7] ); 
		  RHS(2) += 0.5 * ( sflux[i][5] - sflux[eind[0]][6] - tflux[i][5] + tflux[eind[0]][6] ); 
		}
	    }
	  else
	    {
	      RHS(1) += -tflux[i][3] + sflux[i][3];
	      RHS(4) += -tflux[i][4] + sflux[i][4];
	      RHS(2) += -tflux[i][4] + sflux[i][5];
	    }

	  mesh->GetEdgeElements(edges[2], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(2) += 0.5 * ( sflux[i][6] - sflux[eind[1]][2] - tflux[i][6] + tflux[eind[1]][2] );
		  RHS(5) += 0.5 * ( sflux[i][7] - sflux[eind[1]][1] - tflux[i][7] + tflux[eind[1]][1] ); 
		  RHS(0) += 0.5 * ( sflux[i][8] - sflux[eind[1]][0] - tflux[i][8] + tflux[eind[1]][0] ); 
		}
	      else
		{
		  RHS(2) += 0.5 * ( sflux[i][6] - sflux[eind[0]][2] - tflux[i][6] + tflux[eind[0]][2] );
		  RHS(5) += 0.5 * ( sflux[i][7] - sflux[eind[0]][1] - tflux[i][7] + tflux[eind[0]][1] ); 
		  RHS(0) += 0.5 * ( sflux[i][8] - sflux[eind[0]][0] - tflux[i][8] + tflux[eind[0]][0] ); 
		}
	    }
	  else
	    {
	      RHS(2) += -tflux[i][6] + sflux[i][6];
	      RHS(5) += -tflux[i][7] + sflux[i][7];
	      RHS(0) += -tflux[i][8] + sflux[i][8];
	    }
	}
      else
	{
	  mesh->GetEdgeElements(edges[0], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[1]][8] - tflux[i][0] + tflux[eind[1]][8] );
		  RHS(3) += 0.5 * ( sflux[i][1] - sflux[eind[1]][7] - tflux[i][1] + tflux[eind[1]][7] ); 
		  RHS(1) += 0.5 * ( sflux[i][2] - sflux[eind[1]][6] - tflux[i][2] + tflux[eind[1]][6] ); 
		}
	      else
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[0]][8] - tflux[i][0] + tflux[eind[0]][8] );
		  RHS(3) += 0.5 * ( sflux[i][1] - sflux[eind[0]][7] - tflux[i][1] + tflux[eind[0]][7] ); 
		  RHS(1) += 0.5 * ( sflux[i][2] - sflux[eind[0]][6] - tflux[i][2] + tflux[eind[0]][6] ); 
		}
	    }
	  else
	    {
	      RHS(0) += -tflux[i][0] + sflux[i][0];
	      RHS(3) += -tflux[i][1] + sflux[i][1];
	      RHS(1) += -tflux[i][2] + sflux[i][2];
	    }

	  mesh->GetEdgeElements(edges[1], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(1) += 0.5 * ( sflux[i][3] - sflux[eind[1]][2] - tflux[i][3] + tflux[eind[1]][2] );
		  RHS(4) += 0.5 * ( sflux[i][4] - sflux[eind[1]][1] - tflux[i][4] + tflux[eind[1]][1] ); 
		  RHS(2) += 0.5 * ( sflux[i][5] - sflux[eind[1]][0] - tflux[i][5] + tflux[eind[1]][0] ); 
		}
	      else
		{
		  RHS(1) += 0.5 * ( sflux[i][3] - sflux[eind[0]][2] - tflux[i][3] + tflux[eind[0]][2] );
		  RHS(4) += 0.5 * ( sflux[i][4] - sflux[eind[0]][1] - tflux[i][4] + tflux[eind[0]][1] ); 
		  RHS(2) += 0.5 * ( sflux[i][5] - sflux[eind[0]][0] - tflux[i][5] + tflux[eind[0]][0] ); 
		}
	    }
	  else
	    {
	      RHS(1) += -tflux[i][3] + sflux[i][3];
	      RHS(4) += -tflux[i][4] + sflux[i][4];
	      RHS(2) += -tflux[i][5] + sflux[i][5];
	    }

	  mesh->GetEdgeElements(edges[2], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(2) += 0.5 * ( sflux[i][6] - sflux[eind[1]][5] - tflux[i][6] + tflux[eind[1]][5] );
		  RHS(5) += 0.5 * ( sflux[i][7] - sflux[eind[1]][4] - tflux[i][7] + tflux[eind[1]][4] ); 
		  RHS(0) += 0.5 * ( sflux[i][8] - sflux[eind[1]][3] - tflux[i][8] + tflux[eind[1]][3] ); 
		}
	      else
		{
		  RHS(2) += 0.5 * ( sflux[i][6] - sflux[eind[0]][5] - tflux[i][6] + tflux[eind[0]][5] );
		  RHS(5) += 0.5 * ( sflux[i][7] - sflux[eind[0]][4] - tflux[i][7] + tflux[eind[0]][4] ); 
		  RHS(0) += 0.5 * ( sflux[i][8] - sflux[eind[0]][3] - tflux[i][8] + tflux[eind[0]][3] ); 
		}
	    }
	  else
	    {
	      RHS(2) += -tflux[i][6] + sflux[i][6];
	      RHS(5) += -tflux[i][7] + sflux[i][7];
	      RHS(0) += -tflux[i][8] + sflux[i][8];
	    }
	}

      double sum = 0.0;
      for(int k=0; k<6; k++) sum += RHS(k);
      //cout << sum << endl;

      Vector s(6);
      s = 0.0;
      Vector rhs(6); rhs = RHS;
      RectangularMatrix localconsmat(6);      
      for(int j=0; j<6; j++)
	{
	  for(int k=0; k<6; k++)
	    localconsmat(j,k) = locmat->Elem(j,k);
	}
      
      for(int j=0; j<6; j++) localconsmat(3, j) = 0.0;
      localconsmat(3,3) = 1.0;
      rhs(3) = locsol(3);

      DenseMatrixInverse invmat(localconsmat);
      invmat.Mult(rhs, s);
      //s.Print();
      
      Vector sss(6);
      sss = 0.0;
      locmat->Mult(s, sss);
      for(int k=0; k<6; k++)
	cverr(ind[k]) += sss(k) - f(k);
      locmat->Mult(locsol, sss);
      for(int k=0; k<6; k++)
	cverruh(ind[k]) += sss(k) - f(k);
      //cverr(ind[k]) += Q(k) - F[k];
	
      Vector trueflux(9);
      trueflux = 0.0;

      for(int k=0; k<6; k++) { ppsol[k][i] = s(k); }
      for(int k=0; k<9; k++)
	{
	  for(int j=0; j<6; j++)
	    {
	      ppflux[k][i] += -qtem[k][j] * s(j);
	      flux[k][i] += -qtem[k][j] * locsol(j);
	      trueflux(k) += -qtem[k][j] * tsol(ind[j]);
	    }

	  bl2uherr += fabs( trueflux(k) - flux[k][i] ) * fabs( trueflux(k) - flux[k][i] );
	  if( fabs( trueflux(k) - flux[k][i] ) > bmaxuherr )
	    bmaxuherr = fabs( trueflux(k) - flux[k][i] );

	  bl2puherr += fabs( trueflux(k) - ppflux[k][i] ) * fabs( trueflux(k) - ppflux[k][i] );
	  if( fabs( trueflux(k) - ppflux[k][i] ) > bmaxpuherr )
	    bmaxpuherr = fabs( trueflux(k) - ppflux[k][i] );
	}
      
      for(int k=0; k<femat.Size(); k++) delete []femat[k];
      for(int k=0; k<coord.Size(); k++) { delete []coord[k]; }
      for(int k=0; k<qtem.Size(); k++) { delete []qtem[k]; }
      delete locmat;
    }      

  for(int i=0; i<mesh->GetNBE(); i++)
    {
      Array<int> bdrind(2);
      mesh->GetBdrElementVertices(i, bdrind);
      bdrind.Append(mesh->GetNV() +  mesh->GetBdrElementEdgeIndex(i));

      if( fabs(cverr(bdrind[0])) > 1.0e-16 )
	{
	  ppbdrflux[bdrind[0]] = -cverr(bdrind[0]);
	  cverr(bdrind[0]) = 0.0;
	  cverruh(bdrind[0]) = 0.0;
	}
      if( fabs(cverr(bdrind[1])) > 1.0e-16 )
	{
	  ppbdrflux[bdrind[1]] = -cverr(bdrind[1]);
	  cverr(bdrind[1]) = 0.0;
	  cverruh(bdrind[1]) = 0.0;
	}
      
      ppbdrflux[bdrind[2]] = -cverr(bdrind[2]);
      cverr(bdrind[2]) = 0.0;
      cverruh(bdrind[2]) = 0.0;
    }

  double max = 0.0;
  for(int i=0; i<dualmesh->GetDualMeshNumDOF(); i++)
    {
      if(max<fabs(cverr(i)))
	max = fabs(cverr(i));
    }
  cout << "Max LCE: " << max << endl;

  /*
  ofstream fileout("ex11qlce.out");
  for(int i=0; i<dualmesh->GetDualMeshNumDOF(); i++)
    {
      fileout<<i<<"\t"<<cverruh(i)<< "\t" << cverr(i) << endl;
      cout<<i<<"\t"<<cverruh(i)<< "\t" << cverr(i) << endl;
    }
  fileout.close();

  cout << "Edge Flux Max Error  uh: " << bmaxuherr << endl;
  cout << "Edge Flux Max Error puh: " << bmaxpuherr << endl;
  cout << "Edge Flux L2 Error   uh: " << sqrt(bl2uherr) << endl;
  cout << "Edge Flux L2 Error  puh: " << sqrt(bl2puherr) << endl;
  */

  for(int i=0; i<tflux.Size(); i++) delete []tflux[i];
  for(int i=0; i<sflux.Size(); i++) delete []sflux[i];
  for(int i=0; i<pw.Size(); i++) delete []pw[i];
  for(int i=0; i<tpw.Size(); i++) delete []tpw[i];
  for (int k=0; k<cod.Size(); k++) { delete []cod[k]; }
  for (int k=0; k<nmls.Size(); k++) { delete []nmls[k]; }
}

//===========================================================

void Postprocessing::ComputeCubicConservativeTriangularFlux(Vector &tsol, Vector &sol, Array<double*> &flux, Array<double*> &ppsol,
							    Array<double*> &ppflux, Array<double> &ppbdrflux)
{
  mesh->GetEdgeToElementTable();
    
  int nl = NumLocalDOF;
  Array<double *> nmls(18);
  Array<double> lens(18);
  Array<double> areas(10);
  Array<double *> cod(19);
  for(int k=0; k<cod.Size(); k++) cod[k] = new double[2];
  for(int k=0; k<nmls.Size(); k++) nmls[k] = new double[2];

  int n = 5;
  Array<double *> pw(n);
  for(int i=0; i<n; i++) pw[i] = new double[2];
  GetStandLineQuadPW(n, pw);

  int np = ConvertN2PN(6);
  Array<double *> tpw(np);
  for(int i=0; i<np; i++) tpw[i] = new double[3];
  GetStandTriQuadPW(6, tpw);

  Array<double *> sflux(mesh->GetNE());
  for(int i=0; i<sflux.Size(); i++) sflux[i] = new double[12];
  for(int j=0; j<sflux.Size(); j++)
    for(int i=0; i<12; i++)
      sflux[j][i] = 0.0;

  Array<double *> tflux(mesh->GetNE());
  for(int i=0; i<tflux.Size(); i++) tflux[i] = new double[12];
  for(int j=0; j<tflux.Size(); j++)
    for(int i=0; i<12; i++)
      tflux[j][i] = 0.0;

  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      fem->getElementDOF(i, ind);
      Array<int> edges;
      mesh->GetElementEdges(i, edges);

      Array<double *> coord(mesh->GetDim());
      for (int k=0; k<coord.Size(); k++) coord[k] = new double[3];
      
      mesh->GetElementVerticesCoord(i, coord);

      double t1 = coord[0][1] - coord[0][0];
      double t2 = coord[0][2] - coord[0][0];
      double t3 = coord[1][1] - coord[1][0];
      double t4 = coord[1][2] - coord[1][0];
      double detJ = t1 * t4 - t2 * t3;
      
      // Jacobian Inverse Transpose * detJ
      Array<double *> jit(2);
      jit[0] = new double[2];
      jit[1] = new double[2];
      jit[0][0] = coord[1][2] - coord[1][0];
      jit[0][1] = coord[1][0] - coord[1][1];
      jit[1][0] = coord[0][0] - coord[0][2];
      jit[1][1] = coord[0][1] - coord[0][0];

      for(int k=0; k<2; k++)
	for(int j=0; j<2; j++)
	  jit[k][j] /= detJ;
      
      dualmesh->GetElemDualInfo(i, cod, nmls, lens, areas);

      Vector locsol(nl);
      for(int k=0; k<nl; k++) locsol(k) = sol(ind[k]);

      Array<double *> btem(12);
      for(int k=0; k<12; k++) btem[k] = new double[nl];
      for(int k=0; k<12; k++)
	for(int j=0; j<nl; j++)
	  btem[k][j] = 0.0;
      
      double xa, xb, ya, yb, tx, ty, ttx, tty;

      BasisFunctions *phi = new BasisFunctions(3); 
      Function *kk = data->GetEllipticFunction();
      for(int k=0; k<12; k++)
	{
	  double len = 0.0;
	  double nn[2];
	  if(k==0)
	    {
	      xa = coord[0][0];
	      xb = cod[0][0];
	      ya = coord[1][0];
	      yb = cod[0][1];
	    }
	  else if(k==1)
	    {
	      xa = cod[0][0];
	      xb = cod[1][0];
	      ya = cod[0][1];
	      yb = cod[1][1];
	    }
	  else if(k==2)
	    {
	      xa = cod[1][0];
	      xb = cod[2][0];
	      ya = cod[1][1];
	      yb = cod[2][1];
	    }
	  else if(k==3)
	    {
	      xa = cod[2][0];
	      xb = coord[0][1];
	      ya = cod[2][1];
	      yb = coord[1][1];
	    }
	  else if(k==4)
	    {
	      xa = coord[0][1];
	      xb = cod[3][0];
	      ya = coord[1][1];
	      yb = cod[3][1];
	    }
	  else if(k==5)
	    {
	      xa = cod[3][0];
	      xb = cod[4][0];
	      ya = cod[3][1];
	      yb = cod[4][1];
	    }
	  else if(k==6)
	    {
	      xa = cod[4][0];
	      xb = cod[5][0];
	      ya = cod[4][1];
	      yb = cod[5][1];
	    }
	  else if(k==7)
	    {
	      xa = cod[5][0];
	      xb = coord[0][2];
	      ya = cod[5][1];
	      yb = coord[1][2];
	    }
	  else if(k==8)
	    {
	      xa = coord[0][2];
	      xb = cod[6][0];
	      ya = coord[1][2];
	      yb = cod[6][1];
	    }
	  else if(k==9)
	    {
	      xa = cod[6][0];
	      xb = cod[7][0];
	      ya = cod[6][1];
	      yb = cod[7][1];
	    }
	  else if(k==10)
	    {
	      xa = cod[7][0];
	      xb = cod[8][0];
	      ya = cod[7][1];
	      yb = cod[8][1];
	    }
	  else
	    {
	      xa = cod[8][0];
	      xb = coord[0][0];
	      ya = cod[8][1];
	      yb = coord[1][0];
	    }

	  len = sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya) );
	  nn[0] = (yb - ya)/len;
	  nn[1] = (xa - xb)/len;
	  
	  for(int l=0; l<n; l++)
	    {
	      tx = 0.5 * (xb - xa) * pw[l][0] + 0.5 * (xa + xb);
	      ty = 0.5 * (yb - ya) * pw[l][0] + 0.5 * (ya + yb);
	      ttx = jit[0][0] * (tx - coord[0][0]) + jit[1][0] * (ty - coord[1][0]);
	      tty = jit[0][1] * (tx - coord[0][0]) + jit[1][1] * (ty - coord[1][0]);

	      double bcoord[2];
	      bcoord[0] = tx;
	      bcoord[1] = ty;
	      for(int j=0; j<nl; j++)
		{
		  double drx = jit[0][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[0][1] * phi->GradBF2D(j, 1, ttx, tty);
		  double dry = jit[1][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[1][1] * phi->GradBF2D(j, 1, ttx, tty);
		  btem[k][j] += pw[l][1] * kk->Eval(bcoord) * (drx*nn[0] + dry*nn[1]) * 0.5 * len;
		}
	    }
	  for(int l=0; l<nl; l++) sflux[i][k] += btem[k][l] * locsol(l);	  
	}

      for(int k=0; k<12; k++)
	for(int j=0; j<nl; j++)
	  btem[k][j] = 0.0;
      
      for(int k=0; k<12; k++)
	{
	  double len = 0.0;
	  double nn[2];
	  if(k<4)
	    {
	      xa = coord[0][0];
	      xb = coord[0][1];
	      ya = coord[1][0];
	      yb = coord[1][1];
	    }
	  else if(k<8)
	    {
	      xa = coord[0][1];
	      xb = coord[0][2];
	      ya = coord[1][1];
	      yb = coord[1][2];
	    }
	  else
	    {
	      xa = coord[0][2];
	      xb = coord[0][0];
	      ya = coord[1][2];
	      yb = coord[1][0];
	    }

	  len = sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya) );
	  nn[0] = (yb - ya)/len;
	  nn[1] = (xa - xb)/len;
	  
	  for(int l=0; l<n; l++)
	    {
	      tx = 0.5 * (xb - xa) * pw[l][0] + 0.5 * (xa + xb);
	      ty = 0.5 * (yb - ya) * pw[l][0] + 0.5 * (ya + yb);
	      ttx = jit[0][0] * (tx - coord[0][0]) + jit[1][0] * (ty - coord[1][0]);
	      tty = jit[0][1] * (tx - coord[0][0]) + jit[1][1] * (ty - coord[1][0]);

	      double bcoord[2];
	      bcoord[0] = tx;
	      bcoord[1] = ty;
	      for(int j=0; j<nl; j++)
		{
		  double drx = jit[0][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[0][1] * phi->GradBF2D(j, 1, ttx, tty);
		  double dry = jit[1][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[1][1] * phi->GradBF2D(j, 1, ttx, tty);
		  double temprod = pw[l][1] * kk->Eval(bcoord) * (drx*nn[0] + dry*nn[1]) * 0.5 * len;
		  if(k==0)
		    btem[k][j] += phi->BF2D(0, ttx, tty) * temprod;
		  else if(k==1)
		    btem[k][j] += phi->BF2D(3, ttx, tty) * temprod;
		  else if(k==2)
		    btem[k][j] += phi->BF2D(4, ttx, tty) * temprod;
		  else if(k==3)
		    btem[k][j] += phi->BF2D(1, ttx, tty) * temprod;
		  else if(k==4)
		    btem[k][j] += phi->BF2D(1, ttx, tty) * temprod;
		  else if(k==5)
		    btem[k][j] += phi->BF2D(5, ttx, tty) * temprod;
		  else if(k==6)
		    btem[k][j] += phi->BF2D(6, ttx, tty) * temprod;
		  else if(k==7)
		    btem[k][j] += phi->BF2D(2, ttx, tty) * temprod;
		  else if(k==8)
		    btem[k][j] += phi->BF2D(2, ttx, tty) * temprod;
		  else if(k==9)
		    btem[k][j] += phi->BF2D(7, ttx, tty) * temprod;
		  else if(k==10)
		    btem[k][j] += phi->BF2D(8, ttx, tty) * temprod;
		  else if(k==11)
		    btem[k][j] += phi->BF2D(0, ttx, tty) * temprod;
		}
 	    }
	  for(int l=0; l<nl; l++) tflux[i][k] += btem[k][l] * locsol(l);	  
	}

      for (int k=0; k<jit.Size(); k++) { delete []jit[k]; }
      for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
      for (int k=0; k<btem.Size(); k++) { delete []btem[k]; }
    }

  Vector cverr( dualmesh->GetDualMeshNumDOF() );
  cverr = 0.0;
  
  Vector cverruh( dualmesh->GetDualMeshNumDOF() );
  cverruh = 0.0;
  
  for(int i=0; i<mesh->GetNE(); i++)
    {
      Array<int> ind;
      fem->getElementDOF(i, ind);
      Array<int> edges;
      mesh->GetElementEdges(i, edges);

      Array<double *> coord(mesh->GetDim());
      for (int k=0; k<coord.Size(); k++) coord[k] = new double[3];
      
      mesh->GetElementVerticesCoord(i, coord);

      double t1 = coord[0][1] - coord[0][0];
      double t2 = coord[0][2] - coord[0][0];
      double t3 = coord[1][1] - coord[1][0];
      double t4 = coord[1][2] - coord[1][0];
      double detJ = t1 * t4 - t2 * t3;
      
      // Jacobian Inverse Transpose * detJ
      Array<double *> jit(2);
      jit[0] = new double[2];
      jit[1] = new double[2];
      jit[0][0] = coord[1][2] - coord[1][0];
      jit[0][1] = coord[1][0] - coord[1][1];
      jit[1][0] = coord[0][0] - coord[0][2];
      jit[1][1] = coord[0][1] - coord[0][0];

      for(int k=0; k<2; k++)
	for(int j=0; j<2; j++)
	  jit[k][j] /= detJ;
      
      dualmesh->GetElemDualInfo(i, cod, nmls, lens, areas);

      double xa, xb, ya, yb, tx, ty, ttx, tty;
      
      BasisFunctions *phi = new BasisFunctions(3);
      Function *kk = data->GetEllipticFunction();
      Function *ff = data->GetForceFunction();

      Array<double *> ctem(18);
      for(int k=0; k<18; k++) ctem[k] = new double[nl];
      for(int k=0; k<18; k++)
	for(int j=0; j<nl; j++)
	  ctem[k][j] = 0.0;

      for(int k=0; k<18; k++)
	{
	  if(k==0)
	    {
	      xa = cod[9][0];
	      xb = cod[0][0];
	      ya = cod[9][1];
	      yb = cod[0][1];
	    }
	  else if(k==1)
	    {
	      xa = cod[9][0];
	      xb = cod[8][0];
	      ya = cod[9][1];
	      yb = cod[8][1];
	    }
	  else if(k==9)
	    {
	      xa = cod[9][0];
	      xb = cod[12][0];
	      ya = cod[9][1];
	      yb = cod[12][1];
	    }
	  else if(k==2)
	    {
	      xa = cod[10][0];
	      xb = cod[3][0];
	      ya = cod[10][1];
	      yb = cod[3][1];
	    }
	  else if(k==3)
	    {
	      xa = cod[10][0];
	      xb = cod[2][0];
	      ya = cod[10][1];
	      yb = cod[2][1];
	    }
	  else if(k==10)
	    {
	      xa = cod[10][0];
	      xb = cod[14][0];
	      ya = cod[10][1];
	      yb = cod[14][1];
	    }
	  else if(k==4)
	    {
	      xa = cod[11][0];
	      xb = cod[6][0];
	      ya = cod[11][1];
	      yb = cod[6][1];
	    }
	  else if(k==5)
	    {
	      xa = cod[11][0];
	      xb = cod[5][0];
	      ya = cod[11][1];
	      yb = cod[5][1];
	    }
	  else if(k==11)
	    {
	      xa = cod[11][0];
	      xb = cod[16][0];
	      ya = cod[11][1];
	      yb = cod[16][1];
	    }
	  else if(k==6)
	    {
	      xa = cod[13][0];
	      xb = cod[1][0];
	      ya = cod[13][1];
	      yb = cod[1][1];
	    }
	  else if(k==12)
	    {
	      xa = cod[13][0];
	      xb = cod[12][0];
	      ya = cod[13][1];
	      yb = cod[12][1];
	    }
	  else if(k==13)
	    {
	      xa = cod[13][0];
	      xb = cod[14][0];
	      ya = cod[13][1];
	      yb = cod[14][1];
	    }
	  else if(k==7)
	    {
	      xa = cod[15][0];
	      xb = cod[4][0];
	      ya = cod[15][1];
	      yb = cod[4][1];
	    }
	  else if(k==14)
	    {
	      xa = cod[15][0];
	      xb = cod[14][0];
	      ya = cod[15][1];
	      yb = cod[14][1];
	    }
	  else if(k==15)
	    {
	      xa = cod[15][0];
	      xb = cod[16][0];
	      ya = cod[15][1];
	      yb = cod[16][1];
	    }
	  else if(k==8)
	    {
	      xa = cod[17][0];
	      xb = cod[7][0];
	      ya = cod[17][1];
	      yb = cod[7][1];
	    }
	  else if(k==16)
	    {
	      xa = cod[17][0];
	      xb = cod[16][0];
	      ya = cod[17][1];
	      yb = cod[16][1];
	    }
	  else
	    {
	      xa = cod[17][0];
	      xb = cod[12][0];
	      ya = cod[17][1];
	      yb = cod[12][1];
	    }
	  
	  for(int l=0; l<n; l++)
	    {
	      tx = 0.5 * (xb - xa) * pw[l][0] + 0.5 * (xa + xb);
	      ty = 0.5 * (yb - ya) * pw[l][0] + 0.5 * (ya + yb);
	      ttx = jit[0][0] * (tx - coord[0][0]) + jit[1][0] * (ty - coord[1][0]);
	      tty = jit[0][1] * (tx - coord[0][0]) + jit[1][1] * (ty - coord[1][0]);
	      
	      double bcoord[2];
	      bcoord[0] = tx;
	      bcoord[1] = ty;
	      for(int j=0; j<nl; j++)
		{
		  double drx = jit[0][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[0][1] * phi->GradBF2D(j, 1, ttx, tty);
		  double dry = jit[1][0] * phi->GradBF2D(j, 0, ttx, tty) + jit[1][1] * phi->GradBF2D(j, 1, ttx, tty);
		  ctem[k][j] += pw[l][1] * kk->Eval(bcoord) * (drx*nmls[k][0] + dry*nmls[k][1]) * 0.5 * lens[k];
		}
	    }
	}
      
      SparseMatrix *locmat = new SparseMatrix(nl,nl);

      for(int j=0; j<nl; j++)
	locmat->Elem(0,j) = ctem[1][j] - ctem[0][j];
      for(int j=0; j<nl; j++)
	locmat->Elem(1,j) = ctem[3][j] - ctem[2][j];
      for(int j=0; j<nl; j++)
	locmat->Elem(2,j) = ctem[5][j] - ctem[4][j];
      
      for(int j=0; j<nl; j++)
	locmat->Elem(3,j) = ctem[0][j] - ctem[9][j] + ctem[12][j] - ctem[6][j];
      for(int j=0; j<nl; j++)
	locmat->Elem(4,j) = ctem[6][j] - ctem[13][j] + ctem[10][j] - ctem[3][j];

      for(int j=0; j<nl; j++)
	locmat->Elem(5,j) = ctem[2][j] - ctem[10][j] + ctem[14][j] - ctem[7][j];
      for(int j=0; j<nl; j++)
	locmat->Elem(6,j) = ctem[7][j] - ctem[15][j] + ctem[11][j] - ctem[5][j];

      for(int j=0; j<nl; j++)
	locmat->Elem(7,j) = ctem[4][j] - ctem[11][j] + ctem[16][j] - ctem[8][j];
      for(int j=0; j<nl; j++)
	locmat->Elem(8,j) = ctem[8][j] - ctem[17][j] + ctem[9][j] - ctem[1][j];

      for(int j=0; j<nl; j++)
	locmat->Elem(9,j) = ctem[17][j] - ctem[12][j] + ctem[13][j] - ctem[14][j] + ctem[15][j] - ctem[16][j];

      locmat->Finalize();
      //locmat->Print();

      Vector Q(nl);
      Q = 0.0;
      Array<double *> femat(nl);
      Vector locsol(nl);
      for(int k=0; k<nl; k++)
	{
	  femat[k] = new double[nl];
	  locsol(k) = sol(ind[k]);
	}

      fem->GetEllipticLocalSystem(i, ind, femat);
 
      for(int k=0; k<nl; k++)
	for(int j=0; j<nl; j++)
	  Q(k) += femat[k][j] * locsol(j);
	  
      Array<double> F(nl);
      fem->GetForceLocalSystem(i, ind, F);

      Vector f(nl);
      f = 0.0;

      Array<double> x(3);
      Array<double> y(3);
      double gpx, gpy, xx, yy;

      for(int l=0; l<6; l++)
	{
	  if(l==0)
	    {
	      x[0] = coord[0][0];
	      y[0] = coord[1][0];
	      x[1] = cod[0][0];
	      y[1] = cod[0][1];
	      x[2] = cod[8][0];
	      y[2] = cod[8][1];
	    }
	  else if(l==1)
	    {
	      x[0] = cod[0][0];
	      y[0] = cod[0][1];
	      x[2] = cod[8][0];
	      y[2] = cod[8][1];
	      x[1] = cod[9][0];
	      y[1] = cod[9][1];
	    }
	  else if(l==2)
	    {
	      x[1] = coord[0][1];
	      y[1] = coord[1][1];
	      x[0] = cod[2][0];
	      y[0] = cod[2][1];
	      x[2] = cod[3][0];
	      y[2] = cod[3][1];
	    }
	  else if(l==3)
	    {
	      x[1] = cod[2][0];
	      y[1] = cod[2][1];
	      x[2] = cod[3][0];
	      y[2] = cod[3][1];
	      x[0] = cod[10][0];
	      y[0] = cod[10][1];
	    }
	  else if(l==4)
	    {
	      x[0] = coord[0][2];
	      y[0] = coord[1][2];
	      x[1] = cod[6][0];
	      y[1] = cod[6][1];
	      x[2] = cod[5][0];
	      y[2] = cod[5][1];
	    }
	  else
	    {
	      x[0] = cod[6][0];
	      y[0] = cod[6][1];
	      x[2] = cod[5][0];
	      y[2] = cod[5][1];
	      x[1] = cod[11][0];
	      y[1] = cod[11][1];
	    }
	  
	  t1 = x[1] - x[0];
	  t2 = x[2] - x[0];
	  t3 = y[1] - y[0];
	  t4 = y[2] - y[0];
	  detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<np; k++)
	    {
	      gpx = tpw[k][0];
	      gpy = tpw[k][1];
	      xx = x[0] + t1*gpx + t2*gpy;
	      yy = y[0] + t3*gpx + t4*gpy;
	      
	      double bcoord[2];
	      bcoord[0] = xx;
	      bcoord[1] = yy;
	      f(l/2) += tpw[k][2] * ff->Eval(bcoord) * detJ;
	    } 
	}

      for(int l=0; l<24; l++)
	{
	  if(l==0)
	    {
	      x[1] = coord[0][0]*2.0/3.0 + coord[0][1]/3.0;
	      y[1] = coord[1][0]*2.0/3.0 + coord[1][1]/3.0;
	      x[0] = cod[0][0];
	      y[0] = cod[0][1];
	      x[2] = cod[12][0];
	      y[2] = cod[12][1];
	    }
	  else if(l==1)
	    {
	      x[0] = cod[0][0];
	      y[0] = cod[0][1];
	      x[1] = cod[12][0];
	      y[1] = cod[12][1];
	      x[2] = cod[9][0];
	      y[2] = cod[9][1];
	    }
	  else if(l==2)
	    {
	      x[1] = cod[1][0];
	      y[1] = cod[1][1];
	      x[2] = cod[12][0];
	      y[2] = cod[12][1];
	      x[0] = coord[0][0]*2.0/3.0 + coord[0][1]/3.0;
	      y[0] = coord[1][0]*2.0/3.0 + coord[1][1]/3.0;
	    }
	  else if(l==3)
	    {
	      x[0] = cod[1][0];
	      y[0] = cod[1][1];
	      x[2] = cod[12][0];
	      y[2] = cod[12][1];
	      x[1] = cod[13][0];
	      y[1] = cod[13][1];
	    }
	  else if(l==4)
	    {
	      x[1] = cod[2][0];
	      y[1] = cod[2][1];
	      x[2] = cod[14][0];
	      y[2] = cod[14][1];
	      x[0] = coord[0][1]*2.0/3.0 + coord[0][0]/3.0;
	      y[0] = coord[1][1]*2.0/3.0 + coord[1][0]/3.0;
	    }
	  else if(l==5)
	    {
	      x[0] = cod[2][0];
	      y[0] = cod[2][1];
	      x[2] = cod[14][0];
	      y[2] = cod[14][1];
	      x[1] = cod[10][0];
	      y[1] = cod[10][1];
	    }
	  else if(l==6)
	    {
	      x[0] = cod[1][0];
	      y[0] = cod[1][1];
	      x[2] = cod[14][0];
	      y[2] = cod[14][1];
	      x[1] = coord[0][1]*2.0/3.0 + coord[0][0]/3.0;
	      y[1] = coord[1][1]*2.0/3.0 + coord[1][0]/3.0;
	    }
	  else if(l==7)
	    {
	      x[1] = cod[1][0];
	      y[1] = cod[1][1];
	      x[2] = cod[14][0];
	      y[2] = cod[14][1];
	      x[0] = cod[13][0];
	      y[0] = cod[13][1];
	    }
	  else if(l==8)
	    {
	      x[1] = cod[3][0];
	      y[1] = cod[3][1];
	      x[0] = cod[14][0];
	      y[0] = cod[14][1];
	      x[2] = coord[0][1]*2.0/3.0 + coord[0][2]/3.0;
	      y[2] = coord[1][1]*2.0/3.0 + coord[1][2]/3.0;
	    }
	  else if(l==9)
	    {
	      x[1] = cod[3][0];
	      y[1] = cod[3][1];
	      x[2] = cod[14][0];
	      y[2] = cod[14][1];
	      x[0] = cod[10][0];
	      y[0] = cod[10][1];
	    }
	  else if(l==10)
	    {
	      x[1] = cod[4][0];
	      y[1] = cod[4][1];
	      x[2] = cod[14][0];
	      y[2] = cod[14][1];
	      x[0] = coord[0][1]*2.0/3.0 + coord[0][2]/3.0;
	      y[0] = coord[1][1]*2.0/3.0 + coord[1][2]/3.0;
	    }
	  else if(l==11)
	    {
	      x[0] = cod[4][0];
	      y[0] = cod[4][1];
	      x[2] = cod[14][0];
	      y[2] = cod[14][1];
	      x[1] = cod[15][0];
	      y[1] = cod[15][1];
	    }
	  else if(l==12)
	    {
	      x[1] = cod[5][0];
	      y[1] = cod[5][1];
	      x[2] = cod[16][0];
	      y[2] = cod[16][1];
	      x[0] = coord[0][2]*2.0/3.0 + coord[0][1]/3.0;
	      y[0] = coord[1][2]*2.0/3.0 + coord[1][1]/3.0;
	    }
	  else if(l==13)
	    {
	      x[0] = cod[5][0];
	      y[0] = cod[5][1];
	      x[2] = cod[16][0];
	      y[2] = cod[16][1];
	      x[1] = cod[11][0];
	      y[1] = cod[11][1];
	    }
	  else if(l==14)
	    {
	      x[0] = cod[4][0];
	      y[0] = cod[4][1];
	      x[2] = cod[16][0];
	      y[2] = cod[16][1];
	      x[1] = coord[0][2]*2.0/3.0 + coord[0][1]/3.0;
	      y[1] = coord[1][2]*2.0/3.0 + coord[1][1]/3.0;
	    }
	  else if(l==15)
	    {
	      x[1] = cod[4][0];
	      y[1] = cod[4][1];
	      x[2] = cod[16][0];
	      y[2] = cod[16][1];
	      x[0] = cod[15][0];
	      y[0] = cod[15][1];
	    }
	  else if(l==16)
	    {
	      x[0] = cod[6][0];
	      y[0] = cod[6][1];
	      x[2] = cod[16][0];
	      y[2] = cod[16][1];
	      x[1] = coord[0][2]*2.0/3.0 + coord[0][0]/3.0;
	      y[1] = coord[1][2]*2.0/3.0 + coord[1][0]/3.0;
	    }
	  else if(l==17)
	    {
	      x[1] = cod[6][0];
	      y[1] = cod[6][1];
	      x[2] = cod[16][0];
	      y[2] = cod[16][1];
	      x[0] = cod[11][0];
	      y[0] = cod[11][1];
	    }
	  else if(l==18)
	    {
	      x[1] = cod[7][0];
	      y[1] = cod[7][1];
	      x[2] = cod[16][0];
	      y[2] = cod[16][1];
	      x[0] = coord[0][2]*2.0/3.0 + coord[0][0]/3.0;
	      y[0] = coord[1][2]*2.0/3.0 + coord[1][0]/3.0;
	    }
	  else if(l==19)
	    {
	      x[0] = cod[7][0];
	      y[0] = cod[7][1];
	      x[2] = cod[16][0];
	      y[2] = cod[16][1];
	      x[1] = cod[17][0];
	      y[1] = cod[17][1];
	    }
	  else if(l==20)
	    {
	      x[1] = cod[8][0];
	      y[1] = cod[8][1];
	      x[2] = cod[12][0];
	      y[2] = cod[12][1];
	      x[0] = coord[0][0]*2.0/3.0 + coord[0][2]/3.0;
	      y[0] = coord[1][0]*2.0/3.0 + coord[1][2]/3.0;
	    }
	  else if(l==21)
	    {
	      x[0] = cod[8][0];
	      y[0] = cod[8][1];
	      x[2] = cod[12][0];
	      y[2] = cod[12][1];
	      x[1] = cod[9][0];
	      y[1] = cod[9][1];
	    }
	  else if(l==23)
	    {
	      x[0] = cod[7][0];
	      y[0] = cod[7][1];
	      x[2] = cod[12][0];
	      y[2] = cod[12][1];
	      x[1] = coord[0][0]*2.0/3.0 + coord[0][2]/3.0;
	      y[1] = coord[1][0]*2.0/3.0 + coord[1][2]/3.0;
	    }
	  else
	    {
	      x[1] = cod[7][0];
	      y[1] = cod[7][1];
	      x[2] = cod[12][0];
	      y[2] = cod[12][1];
	      x[0] = cod[17][0];
	      y[0] = cod[17][1];
	    }

	  t1 = x[1] - x[0];
	  t2 = x[2] - x[0];
	  t3 = y[1] - y[0];
	  t4 = y[2] - y[0];
	  detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<np; k++)
	    {
	      gpx = tpw[k][0];
	      gpy = tpw[k][1];
	      xx = x[0] + t1*gpx + t2*gpy;
	      yy = y[0] + t3*gpx + t4*gpy;

	      double bcoord[2];
	      bcoord[0] = xx;
	      bcoord[1] = yy;
	      f(3 + l/4) += tpw[k][2] * ff->Eval(bcoord) * detJ;
	    } 
	}

      for(int l=0; l<6; l++)
	{
	  x[0] = cod[18][0];
	  y[0] = cod[18][1];
	  x[1] = cod[12+l][0];
	  y[1] = cod[12+l][1];
	  x[2] = cod[12+(l+1)%6][0];
	  y[2] = cod[12+(l+1)%6][1];
	  
	  t1 = x[1] - x[0];
	  t2 = x[2] - x[0];
	  t3 = y[1] - y[0];
	  t4 = y[2] - y[0];
	  detJ = t1 * t4 - t2 * t3;
	  
	  for(int k=0; k<np; k++)
	    {
	      gpx = tpw[k][0];
	      gpy = tpw[k][1];
	      xx = x[0] + t1*gpx + t2*gpy;
	      yy = y[0] + t3*gpx + t4*gpy;

	      double bcoord[2];
	      bcoord[0] = xx;
	      bcoord[1] = yy;
	      f(9) += tpw[k][2] * ff->Eval(bcoord) * detJ;
	    } 
	}

      Vector RHS(nl);
      for(int k=0; k<nl; k++) RHS(k) = Q(k) - F[k] + f(k);

      if(i%2==0)
	{
	  Array<int> eind;
	  mesh->GetEdgeElements(edges[0], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[1]][7] - tflux[i][0] + tflux[eind[1]][7] );
		  RHS(3) += 0.5 * ( sflux[i][1] - sflux[eind[1]][6] - tflux[i][1] + tflux[eind[1]][6] ); 
		  RHS(4) += 0.5 * ( sflux[i][2] - sflux[eind[1]][5] - tflux[i][2] + tflux[eind[1]][5] );
		  RHS(1) += 0.5 * ( sflux[i][3] - sflux[eind[1]][4] - tflux[i][3] + tflux[eind[1]][4] ); 
		}
	      else
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[0]][7] - tflux[i][0] + tflux[eind[0]][7] );
		  RHS(3) += 0.5 * ( sflux[i][1] - sflux[eind[0]][6] - tflux[i][1] + tflux[eind[0]][6] ); 
		  RHS(4) += 0.5 * ( sflux[i][2] - sflux[eind[0]][5] - tflux[i][2] + tflux[eind[0]][5] );
		  RHS(1) += 0.5 * ( sflux[i][3] - sflux[eind[0]][4] - tflux[i][3] + tflux[eind[0]][4] ); 
		}
	    }
	  else
	    {
	      RHS(0) += -tflux[i][0] + sflux[i][0];
	      RHS(3) += -tflux[i][1] + sflux[i][1];
	      RHS(4) += -tflux[i][2] + sflux[i][2];
	      RHS(1) += -tflux[i][3] + sflux[i][3];
	    }

	  mesh->GetEdgeElements(edges[1], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(1) += 0.5 * ( sflux[i][4] - sflux[eind[1]][11] - tflux[i][4] + tflux[eind[1]][11] );
		  RHS(5) += 0.5 * ( sflux[i][5] - sflux[eind[1]][10] - tflux[i][5] + tflux[eind[1]][10] ); 
		  RHS(6) += 0.5 * ( sflux[i][6] - sflux[eind[1]][9] - tflux[i][6] + tflux[eind[1]][9] );
		  RHS(2) += 0.5 * ( sflux[i][7] - sflux[eind[1]][8] - tflux[i][7] + tflux[eind[1]][8] ); 
		}
	      else
		{
		  RHS(1) += 0.5 * ( sflux[i][4] - sflux[eind[0]][11] - tflux[i][4] + tflux[eind[0]][11] );
		  RHS(5) += 0.5 * ( sflux[i][5] - sflux[eind[0]][10] - tflux[i][5] + tflux[eind[0]][10] ); 
		  RHS(6) += 0.5 * ( sflux[i][6] - sflux[eind[0]][9] - tflux[i][6] + tflux[eind[0]][9] );
		  RHS(2) += 0.5 * ( sflux[i][7] - sflux[eind[0]][8] - tflux[i][7] + tflux[eind[0]][8] ); 
		}
	    }
	  else
	    {
	      RHS(1) += -tflux[i][4] + sflux[i][4];
	      RHS(5) += -tflux[i][5] + sflux[i][5];
	      RHS(6) += -tflux[i][6] + sflux[i][6];
	      RHS(2) += -tflux[i][7] + sflux[i][7];
	    }

	  mesh->GetEdgeElements(edges[2], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(2) += 0.5 * ( sflux[i][8] - sflux[eind[1]][3] - tflux[i][8] + tflux[eind[1]][3] );
		  RHS(7) += 0.5 * ( sflux[i][9] - sflux[eind[1]][2] - tflux[i][9] + tflux[eind[1]][2] ); 
		  RHS(8) += 0.5 * ( sflux[i][10] - sflux[eind[1]][1] - tflux[i][10] + tflux[eind[1]][1] );
		  RHS(0) += 0.5 * ( sflux[i][11] - sflux[eind[1]][0] - tflux[i][11] + tflux[eind[1]][0] ); 
		}
	      else
		{
		  RHS(2) += 0.5 * ( sflux[i][8] - sflux[eind[0]][3] - tflux[i][8] + tflux[eind[0]][3] );
		  RHS(7) += 0.5 * ( sflux[i][9] - sflux[eind[0]][2] - tflux[i][9] + tflux[eind[0]][2] ); 
		  RHS(8) += 0.5 * ( sflux[i][10] - sflux[eind[0]][1] - tflux[i][10] + tflux[eind[0]][1] );
		  RHS(0) += 0.5 * ( sflux[i][11] - sflux[eind[0]][0] - tflux[i][11] + tflux[eind[0]][0] ); 
		}
	    }
	  else
	    {
	      RHS(2) += -tflux[i][8] + sflux[i][8];
	      RHS(7) += -tflux[i][9] + sflux[i][9];
	      RHS(8) += -tflux[i][10] + sflux[i][10];
	      RHS(0) += -tflux[i][11] + sflux[i][11];
	    }
	}
      else
	{
	  Array<int> eind;
	  mesh->GetEdgeElements(edges[0], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[1]][11] - tflux[i][0] + tflux[eind[1]][11] );
		  RHS(3) += 0.5 * ( sflux[i][1] - sflux[eind[1]][10] - tflux[i][1] + tflux[eind[1]][10] ); 
		  RHS(4) += 0.5 * ( sflux[i][2] - sflux[eind[1]][9] - tflux[i][2] + tflux[eind[1]][9] );
		  RHS(1) += 0.5 * ( sflux[i][3] - sflux[eind[1]][8] - tflux[i][3] + tflux[eind[1]][8] ); 
		}
	      else
		{
		  RHS(0) += 0.5 * ( sflux[i][0] - sflux[eind[0]][11] - tflux[i][0] + tflux[eind[0]][11] );
		  RHS(3) += 0.5 * ( sflux[i][1] - sflux[eind[0]][10] - tflux[i][1] + tflux[eind[0]][10] ); 
		  RHS(4) += 0.5 * ( sflux[i][2] - sflux[eind[0]][9] - tflux[i][2] + tflux[eind[0]][9] );
		  RHS(1) += 0.5 * ( sflux[i][3] - sflux[eind[0]][8] - tflux[i][3] + tflux[eind[0]][8] ); 
		}
	    }
	  else
	    {
	      RHS(0) += -tflux[i][0] + sflux[i][0];
	      RHS(3) += -tflux[i][1] + sflux[i][1];
	      RHS(4) += -tflux[i][2] + sflux[i][2];
	      RHS(1) += -tflux[i][3] + sflux[i][3];
	    }

	  mesh->GetEdgeElements(edges[1], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(1) += 0.5 * ( sflux[i][4] - sflux[eind[1]][3] - tflux[i][4] + tflux[eind[1]][3] );
		  RHS(5) += 0.5 * ( sflux[i][5] - sflux[eind[1]][2] - tflux[i][5] + tflux[eind[1]][2] ); 
		  RHS(6) += 0.5 * ( sflux[i][6] - sflux[eind[1]][1] - tflux[i][6] + tflux[eind[1]][1] );
		  RHS(2) += 0.5 * ( sflux[i][7] - sflux[eind[1]][0] - tflux[i][7] + tflux[eind[1]][0] ); 
		}
	      else
		{
		  RHS(1) += 0.5 * ( sflux[i][4] - sflux[eind[0]][3] - tflux[i][4] + tflux[eind[0]][3] );
		  RHS(5) += 0.5 * ( sflux[i][5] - sflux[eind[0]][2] - tflux[i][5] + tflux[eind[0]][2] ); 
		  RHS(6) += 0.5 * ( sflux[i][6] - sflux[eind[0]][1] - tflux[i][6] + tflux[eind[0]][1] );
		  RHS(2) += 0.5 * ( sflux[i][7] - sflux[eind[0]][0] - tflux[i][7] + tflux[eind[0]][0] ); 
		}
	    }
	  else
	    {
	      RHS(1) += -tflux[i][4] + sflux[i][4];
	      RHS(5) += -tflux[i][5] + sflux[i][5];
	      RHS(6) += -tflux[i][6] + sflux[i][6];
	      RHS(2) += -tflux[i][7] + sflux[i][7];
	    }

	  mesh->GetEdgeElements(edges[2], eind);	  
	  if( eind.Size()>1 )
	    {
	      if( eind[0]==i )
		{
		  RHS(2) += 0.5 * ( sflux[i][8] - sflux[eind[1]][7] - tflux[i][8] + tflux[eind[1]][7] );
		  RHS(7) += 0.5 * ( sflux[i][9] - sflux[eind[1]][6] - tflux[i][9] + tflux[eind[1]][6] ); 
		  RHS(8) += 0.5 * ( sflux[i][10] - sflux[eind[1]][5] - tflux[i][10] + tflux[eind[1]][5] );
		  RHS(0) += 0.5 * ( sflux[i][11] - sflux[eind[1]][4] - tflux[i][11] + tflux[eind[1]][4] ); 
		}
	      else
		{
		  RHS(2) += 0.5 * ( sflux[i][8] - sflux[eind[0]][7] - tflux[i][8] + tflux[eind[0]][7] );
		  RHS(7) += 0.5 * ( sflux[i][9] - sflux[eind[0]][6] - tflux[i][9] + tflux[eind[0]][6] ); 
		  RHS(8) += 0.5 * ( sflux[i][10] - sflux[eind[0]][5] - tflux[i][10] + tflux[eind[0]][5] );
		  RHS(0) += 0.5 * ( sflux[i][11] - sflux[eind[0]][4] - tflux[i][11] + tflux[eind[0]][4] ); 
		}
	    }
	  else
	    {
	      RHS(2) += -tflux[i][8] + sflux[i][8];
	      RHS(7) += -tflux[i][9] + sflux[i][9];
	      RHS(8) += -tflux[i][10] + sflux[i][10];
	      RHS(0) += -tflux[i][11] + sflux[i][11];
	    }
	}

      Vector s(nl);
      s = 0.0;
      
      RectangularMatrix localconsmat(nl);
      
      for(int j=0; j<nl; j++)
	for(int k=0; k<nl; k++)
	  localconsmat(j,k) = locmat->Elem(j,k);

      for(int j=0; j<nl; j++) localconsmat(9, j) = 0.0;
      localconsmat(9,9) = 1.0;
      RHS(9) = locsol(9); 

      DenseMatrixInverse invmat(localconsmat);
      invmat.Mult(RHS, s);
      //s.Print();

      Vector sss(10);
      sss = 0.0;
      locmat->Mult(s, sss);
      for(int k=0; k<10; k++) cverr(ind[k]) += sss(k) - f(k);
      locmat->Mult(locsol, sss);
      for(int k=0; k<10; k++) cverruh(ind[k]) += sss(k) - f(k);
      
      for(int k=0; k<nl; k++)
	{
	  ppsol[k][i] = s(k);
	}

      for(int k=0; k<18; k++)
	{
	  for(int j=0; j<nl; j++)
	    {
	      ppflux[k][i] += -ctem[k][j] * s(j);
	      flux[k][i] += -ctem[k][j] * locsol(j);
	    }
	}

      for(int k=0; k<femat.Size(); k++) delete []femat[k];
      for(int k=0; k<jit.Size(); k++) delete []jit[k];
      for (int k=0; k<coord.Size(); k++) { delete []coord[k]; }
      for (int k=0; k<ctem.Size(); k++) { delete []ctem[k]; }
      delete locmat;
    }

  for(int i=0; i<mesh->GetNBE(); i++)
    {
      Array<int> bdrind(2);
      mesh->GetBdrElementVertices(i, bdrind);

      if(bdrind[0] < bdrind[1])
	{
	  bdrind.Append( mesh->GetNV() + 2 * mesh->GetBdrElementEdgeIndex(i) );
	  bdrind.Append( mesh->GetNV() + 2 * mesh->GetBdrElementEdgeIndex(i) + 1);
	}
      else
	{
	  bdrind.Append( mesh->GetNV() + 2 * mesh->GetBdrElementEdgeIndex(i) + 1);
	  bdrind.Append( mesh->GetNV() + 2 * mesh->GetBdrElementEdgeIndex(i) );
	}

      if( fabs(cverr(bdrind[0])) > 1.0e-15 )
	{
	  ppbdrflux[bdrind[0]] = -cverr(bdrind[0]);
	  cverr(bdrind[0]) = 0.0;
	  cverruh(bdrind[0]) = 0.0;
	}
      if( fabs(cverr(bdrind[1])) > 1.0e-15 )
	{
	  ppbdrflux[bdrind[1]] = -cverr(bdrind[1]);
	  cverr(bdrind[1]) = 0.0;
	  cverruh(bdrind[1]) = 0.0;
	}
      
      ppbdrflux[bdrind[2]] = -cverr(bdrind[2]);
      cverr(bdrind[2]) = 0.0;
      cverruh(bdrind[2]) = 0.0;
      ppbdrflux[bdrind[3]] = -cverr(bdrind[3]);
      cverr(bdrind[3]) = 0.0;
      cverruh(bdrind[3]) = 0.0;
    }

  ofstream fileout("clce.out");
  for(int i=0; i<dualmesh->GetDualMeshNumDOF(); i++)
    {
      fileout<<i<<"\t"<<cverruh(i) << "\t" << cverr(i) <<endl;
      cout<<i<<"\t"<<cverruh(i) << "\t"<<cverr(i)<<endl;
    }
  fileout.close(); 
  
  for(int i=0; i<tflux.Size(); i++) delete []tflux[i];
  for(int i=0; i<sflux.Size(); i++) delete []sflux[i];
  for(int i=0; i<pw.Size(); i++) delete []pw[i];
  for(int i=0; i<tpw.Size(); i++) delete []tpw[i];
  for (int k=0; k<cod.Size(); k++) { delete []cod[k]; }
  for (int k=0; k<nmls.Size(); k++) { delete []nmls[k]; }
}

//===========================================================

Postprocessing::~Postprocessing()
{
  ;
}

