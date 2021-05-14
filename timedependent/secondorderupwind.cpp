#include "timedependent_header.h"
/*
SecondOrderUpwind::SecondOrderUpwind(int _Nx, int _Ny, double _Hx, double _Hy,
				     Coefficient *f) :
  HyperbolicProblem(_Nx, _Ny, _Hx, _Hy, f)
{
  ;
}
*/
SecondOrderUpwind::SecondOrderUpwind(DualMesh *_dualmesh, Array<Function*> &_fluxfunction):
  HyperbolicProblem(_dualmesh, _fluxfunction)
{
  ;
}

//=============================================================

double SecondOrderUpwind::FluxContribution(Vector &sol, Array<int> &dof, FData &F, int el, int k, int b)
{
  /*
  Array<int> ind(4);
  mesh->GetElementVertices(el, ind);
  int cb = ind[0] + j*(mesh->GetNx()+1);
  /*
  F.SetSize(4);
  double val[1], sig[1];
  int c = i + j*(mesh->GetNx()+1);

  // Left edge of control volume.
  if (-flux[0][c]<0.0)
    {
      val[0] = (i!=0) ? sol(c-1) : sol(c);
      sig[0] = (i>1) ? 0.5 * SlopeInfo_x(i-1, j) : 0.0;
    }
  else
    {
      val[0] = sol(c);
      sig[0] = (i>0 && i<mesh->GetNx()) ? -0.5 * SlopeInfo_x(i,j) : 0.0;
      val[0] += sig[0];
    }
  F[0] = -flux[0][c] * FluxFunction->Eval(val);


  // Right edge of control volume.
  if (flux[1][c]<0.0)
    {
      val[0] = (i!=mesh->GetNx()) ? sol(c+1) : sol(c);
      sig[0] = (i<mesh->GetNx()-1) ? -0.5 * SlopeInfo_x(i+1, j) : 0.0;
      val[0] += sig[0];
    }
  else
    {
      val[0] = sol(c);
      sig[0] = (i>0 && i<mesh->GetNx()) ? 0.5 * SlopeInfo_x(i,j) : 0.0;
      val[0] += sig[0];
    }
  F[1] = flux[1][c] * FluxFunction->Eval(val);

  // Bottom edge of control volume.
  if (-flux[2][c]<0.0)
    {
      val[0] = (j!=0) ? sol(c-mesh->GetNx()-1) : sol(c);
      sig[0] = (j>1) ? 0.5 * SlopeInfo_y(i, j-1) : 0.0;
      val[0] += sig[0];
    }
  else
    {
      val[0] = sol(c);
      sig[0] = (j>0 && j<mesh->GetNy()) ? -0.5 * SlopeInfo_y(i,j) : 0.0;
      val[0] += sig[0];
    }
  F[2] = -flux[2][c] * FluxFunction->Eval(val);

  // Top edge of control volume.
  if (flux[3][c]<0.0)
    {
      val[0] = (j!=mesh->GetNy()) ? sol(c+mesh->GetNx()+1) : sol(c);
      sig[0] = (j<mesh->GetNy()-1) ? -0.5 * SlopeInfo_y(i, j+1) : 0.0;
      val[0] += sig[0];
    }
  else
    {
      val[0] = sol(c);
      sig[0] = (j>0 && j<mesh->GetNy()) ? 0.5 * SlopeInfo_y(i,j) : 0.0;
      val[0] += sig[0];
    }
  F[3] = flux[3][c] * FluxFunction->Eval(val);
  return 0.0;*/;
}

//=============================================================

double SecondOrderUpwind::MinMod(double a, double b)
{
  double ma = fabs(a);
  double mb = fabs(b);
  double min_a_b = (ma<mb) ? ma : mb;
  double signa = ((a<0.0) ? -1 : 1);
  double signb = ((b<0.0) ? -1 : 1);
  double fact = (double) (signa + signb);
  return (0.5 * fact * min_a_b);
} 

//=============================================================

double SecondOrderUpwind::SlopeInfo_x(Vector &sol, int i, int j)
{
  /*
  int c = i+(mesh->GetNx()+1)*j;
  int l = c-1;
  int r = c+1;
  double a = sol(r) - sol(c);
  double b = sol(c) - sol(l);
  return MinMod(a, b);*/;
}

//============================================================================

double SecondOrderUpwind::SlopeInfo_y(Vector &sol, int i, int j)
{
  /*
  int c = i+(mesh->GetNx()+1)*j;
  int b = c-mesh->GetNx()-1;
  int t = c+mesh->GetNx()+1;
  double x = sol(t) - sol(c);
  double y = sol(c) - sol(b);
  return MinMod(x, y);*/;
}
