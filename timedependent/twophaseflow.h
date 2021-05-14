#ifndef FILE_TWOPHASEFLOW
#define FILE_TWOPHASEFLOW

/*
  Author: Q. Deng
  Data: 10/01/2015
*/

using namespace std;

//============================================================================
double MinMod(double a, double b);
double SlopeInfo_0(int i, Mesh *mesh, Vector &sol);
double SlopeInfo_1(int i, Mesh *mesh, Vector &sol);
double SlopeInfo_2(int i, Mesh *mesh, Vector &sol);
  
//============================================================================
// Produce local internal element flux using upwind;
void LocalLinearTriangularEUF(Array<double> &locflux, Array<double> &locsat, 
			      Function *fff, Array<double> &locupwindflux);

void LocalLinearTriangular2ndOrderEUF(int ind, Mesh *mesh, Array<double> &locflux, Vector &sat, 
				      Function *fff, Array<double> &locupwindflux);

//============================================================================
// Produce local internal element flux using upwind;
void LocalLinearRectangularEUF(Array<double> &locflux, Array<double> &locsat, 
			       Function *fff, Array<double> &locupwindflux);

//============================================================================
// Produce local internal element flux using upwind;
void LocalQuadraticTriangularEUF(Array<double> &locflux, Array<double> &locsat, 
				 Function *fff, Array<double> &locupwindflux);

//============================================================================
// Produce local internal element flux using upwind;
void LocalQuadraticRectangularEUF(Array<double> &locflux, Array<double> &locsat, 
				  Function *fff, Array<double> &locupwindflux);

//============================================================================
// Produce local internal element flux using upwind;
void LocalCubicTriangularEUF(Array<double> &locflux, Array<double> &locsat, 
			     Function *fff, Array<double> &locupwindflux);

//============================================================================
// Produce local internal element flux using upwind;
void LocalCubicRectangularEUF(Array<double> &locflux, Array<double> &locsat, 
			      Function *fff, Array<double> &locupwindflux);

//============================================================================
// Produce internal element flux using upwind; fff: fractional flow function
void ElementUpwindFlux(DualMesh *dualmesh, Array<double*> &flux, Vector &sat, 
		       Function *fff, Array<double> &upwindflux);

//============================================================================
// Produce global boundary flux without using upwind?
void BdrElementUpwindFlux(Array<double*> &bdrflux, Vector &sat, 
			  Function *fff, Array<double> &upwindbdrflux);

//============================================================================
// TPF: Two Phase Flow 
void TPFTimeMarchSolve(double dt, int numfinetimestep, DualMesh *dualmesh, Array<double> &areas,
		       Array<double*> &flux, Array<double> &bdrflux,
		       Function *fff, Vector &satold, Vector &satnew);

#endif
