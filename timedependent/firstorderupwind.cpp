#include "timedependent_header.h"

//=============================================================

FirstOrderUpwind::FirstOrderUpwind(DualMesh *_dualmesh, Array<Function*> &_fluxfunction):
  HyperbolicProblem(_dualmesh, _fluxfunction)
{
  ;
}

//=============================================================

double FirstOrderUpwind::FluxContribution(Vector &Solold, Array<int> &dof, FData &F, int el, int k, int b)
{
  /*
    Smbolically, as a typical example return $\lambda(S)* \int_{\partial C_v \cap \tau} v_D \cdot n dl $
  */
  double value = 0.0;
  int numdof = dof.Size();

  //upwinding
  double fluxvalue = F.GetFVal(k, el);
  int ind = (0.0<fluxvalue) ? k : (k+1)%numdof;
  double solval[1];
  solval[0] = Solold(dof[ind]);

  double fluxfuncval = fluxfunction[b]->Eval(solval);
  value += fluxfuncval*fluxvalue; //first boundary


  //upwinding
  int kk = (k+numdof-1)%numdof;
  fluxvalue = F.GetFVal(kk, el);
  ind = (0.0<fluxvalue) ? kk : k;
  solval[0] = Solold(dof[ind]);

  fluxfuncval = fluxfunction[b]->Eval(solval);
  value -= fluxfuncval*fluxvalue; //second boundary

  return value;
}


//=============================================================
double FirstOrderUpwind::BdrFluxContribution(Vector &Solold, FData &F, int vind, int el, int locind, int numv, int b)
{
  /*
    Smbolically, as a typical example return $\lambda(S)* \int_{\partial C_v \cap \tau} v_D \cdot n dl $
  */
  double value = 0.0;

  //upwinding
  double fluxvalue = F.GetFVal(locind, el);
  double solval[1];
  solval[0] = Solold(vind);

  double fluxfuncval = fluxfunction[b]->Eval(solval);
  value += fluxfuncval*fluxvalue; //first boundary

  //upwinding
  int kk = (locind+numv-1)%numv;
  fluxvalue = F.GetFVal(kk, el);
  solval[0] = Solold(vind);

  fluxfuncval = fluxfunction[b]->Eval(solval);
  value -= fluxfuncval*fluxvalue; //second boundary

  return value;
}

