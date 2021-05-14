#ifndef FILE_POSTPROCESSING
#define FILE_POSTPROCESSING

/*
   Data type for postprocessing

   Author: Quanling Deng
   Date: 09/01/2015
*/

/* This is a container for postprocessing techniques to 
   obtain the locally conservative fluxes */
class Postprocessing
{
protected:

  Mesh *mesh;
  
  DualMesh *dualmesh;

  FEM *fem;

  Data *data;

  int NumLocalDOF;
  
public:

  Postprocessing(DualMesh *_dualmesh, FEM *_fem, Data *_data);

  void ComputeConservativeFlux(Vector &tsol, Vector &sol, Array<double*> &flux, Array<double *> &ppsol,
			       Array<double *> &ppflux, Array<double> &ppbdrflux);

  void ComputeLinearConservativeTriangularFlux(Vector &tsol, Vector &sol, Array<double*> &flux, Array<double*> &ppsol,
					       Array<double*> &ppflux, Array<double> &ppbdrflux);

  void ComputeQuadraticConservativeTriangularFlux(Vector &tsol, Vector &sol, Array<double*> &flux, Array<double*> &ppsol,
						  Array<double*> &ppflux, Array<double> &ppbdrflux);

  void ComputeCubicConservativeTriangularFlux(Vector &tsol, Vector &sol, Array<double*> &flux, Array<double*> &ppsol,
					      Array<double*> &ppflux, Array<double> &ppbdrflux);

  void ComputeLinearConservativeRectangularFlux(Vector &tsol, Vector &sol, Array<double*> &flux, Array<double*> &ppsol,
						Array<double*> &ppflux, Array<double> &ppbdrflux) { ; }

  void ComputeQuadraticConservativeRectangularFlux(Vector &tsol, Vector &sol, Array<double*> &flux, Array<double*> &ppsol,
						   Array<double*> &ppflux, Array<double> &ppbdrflux) { ; }

  void ComputeCubicConservativeRectangularFlux(Vector &tsol, Vector &sol, Array<double*> &flux, Array<double*> &ppsol,
					       Array<double*> &ppflux, Array<double> &ppbdrflux) { ; }

  ~Postprocessing();
};

#endif
