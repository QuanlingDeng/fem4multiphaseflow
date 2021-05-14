#ifndef FILE_HYPERBOLICPROBLEM
#define FILE_HYPERBOLICPROBLEM

/*
  Author: Q. Deng
  Data: 10/01/2015
*/

class HyperbolicProblem
{
protected:

  Array<double *> flux;
  Array<Function *> fluxfunction;
  Vector sol;
  DualMesh *dualmesh;

  virtual double FluxContribution(Vector &Solold, Array<int> &dof, FData &F, int el, int k, int b) = 0;

  virtual double BdrFluxContribution(Vector &Solold, FData &F, int vind, int el, int locind, int numv, int b) = 0;

  virtual void ConvectiveFlux(Array<double> &CF);

public:
  
  HyperbolicProblem(DualMesh *dualmesh, Array<Function*> &fluxfunction);

  void UpdateFlux(const Array<double *> &_flux);

  void TimeMarch(Vector &initsol, Vector &sol, Array<FData *> &F, const Array<double> &time, Vector &porosity);

  void PoroTimeMarch(Vector &initsol, Vector &sol, Array<FData *> &F, const Array<double> &time, Vector &porosity);

  void Accumulate(Vector &Accum, Vector &solold, FData &F, double dt, int b);

  void SSPRK2(Vector &initsol, Vector &sol, Vector &porosity, Array<FData *> &F, const Array<double> &time);

  virtual ~HyperbolicProblem();
};


class FirstOrderUpwind : public HyperbolicProblem
{
 protected:
  
  virtual double FluxContribution(Vector &sol, Array<int> &dof, FData &F, int el, int k, int b);

  virtual double BdrFluxContribution(Vector &Solold, FData &F, int vind, int el, int locind, int numv, int b);
  
 public:
  
  FirstOrderUpwind(DualMesh *_dualmesh, Array<Function*> &_fluxfunction);


};


class SecondOrderUpwind : public HyperbolicProblem
{
 protected:
  
  virtual double FluxContribution(Vector &Solold, Array<int> &dof, FData &F, int el, int k, int b);

  virtual double BdrFluxContribution(Vector &Solold, FData &F, int vind, int el, int locind, int numv, int b);

  double MinMod(double a, double b);

  double SlopeInfo_x(Vector &sol, int i, int j);

  double SlopeInfo_y(Vector &sol, int i, int j);

 public:

  SecondOrderUpwind(DualMesh *_dualmesh, Array<Function*> &_fluxfunction);

};


#endif
