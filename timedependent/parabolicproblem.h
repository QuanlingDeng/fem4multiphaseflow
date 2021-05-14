#ifndef FILE_PARABOLICPROBLEM
#define FILE_PARABOLICPROBLEM

/*
  Author: Q. Deng
  Data: 10/15/2015
*/

class ParabolicProblem
{
protected:

  void SetDiagonalMatrix(Array<bool *> &dirichlet, double ts, Array<double> &diag,
					   Vector &rhs);

public:
  
  ParabolicProblem();

  ~ParabolicProblem();

  void TimeMarchSolve(Array<double *> &data, Array<double *> &bdrydata,
		      Array<bool *> &dirichlet, Array<double *> &bdryreaction,
                      Vector &initsol, Array<double> &time, Vector &sol, Vector &temp);

  void TimeMarchSolveConsistentMass(Array<double *> &data, Array<double *> &bdrydata,
				    Array<bool *> &dirichlet, Array<double *> &bdryreaction,
                                    Vector &initsol, Array<double> &time, Vector &sol, Vector &temp);



};

#endif
