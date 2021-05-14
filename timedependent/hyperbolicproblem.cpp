#include "timedependent_header.h"

//=============================================================

HyperbolicProblem::HyperbolicProblem(DualMesh *_dualmesh, Array<Function *> &f)
{
  dualmesh = _dualmesh;
  fluxfunction.SetSize(f.Size());
  for(int b=0; b<f.Size(); b++){fluxfunction[b] = f[b]; }
}

//=============================================================

HyperbolicProblem::~HyperbolicProblem()
{
  ;
}

//=============================================================

void HyperbolicProblem::UpdateFlux(const Array<double *> &_flux)
{
  for (int i=0; i<4; i++) { flux[i] = _flux[i]; }
}

//=============================================================

void HyperbolicProblem::ConvectiveFlux(Array<double> &CF)
{
  /*
  Array<double> F;
  double cvx, cvy;

  int nv = mesh->GetNV();
  
  CF.SetSize(nv);
  int k = 0;
  for(int k=0; k<nv; k++)
  for (int j=0; j<=mesh->GetNy(); j++)
    {
      cvy = mesh->GetCVLengthy(j);
      for (int i=0; i<=mesh->GetNx(); i++)
	{
	  cvx = mesh->GetCVLengthx(i);
	  FluxContribution(i, j, F);
	  CF[k] = (F[0]+F[1])/cvx + (F[2]+F[3])/cvy;   
	  k++;
	}
    }
*/
  
}

//=============================================================

void HyperbolicProblem::TimeMarch(Vector &initsol, Vector &sol, Array<FData *> &F, const Array<double> &time, Vector &porosity)
{
  /*
   Given an initial profile (initsol) march through the time steps (time), and numerically solve the 
   hyperbolic problem. The Flux data is assumed to already be integrated over "Sub"-control volume 
   boundaries. Symbolically, S^n = S^{n-1} - c*F. 
   */

  sol = initsol;
  bool test = false;

  for(int n=0; n<time.Size()-1; n++)
    {
      double dt = time[n+1] - time[n];
      Vector Accum(sol.Size());
      for(int b=0; b<F.Size(); b++)
	{
	  Accumulate(Accum, sol, *F[b], dt, b);
	  for(int k=0; k<sol.Size(); k++)
	    {
	      sol(k) -= Accum(k)/porosity(k);
	      //cout << k << " " << Accum(k) << endl;
	      if(sol(k) < -1.0e-7 || sol(k) > 1.0+1.0e-5){ cout << "saturation at node " << k << " is " << sol(k) << endl; test = true; }
	    }
	}  
      if(test==true){break;}
    }
}

//=============================================================

void HyperbolicProblem::PoroTimeMarch(Vector &initsol, Vector &sol, Array<FData *> &F, const Array<double> &time, Vector &porosity)
{
  /*
   Given an initial profile (initsol) march through the time steps (time), and numerically solve the 
   hyperbolic problem. The Flux data is assumed to already be integrated over "Sub"-control volume 
   boundaries. Symbolically, S^n = S^{n-1} - c*F. 
   */

  sol = initsol;
  bool test = false;

  cout << endl;
  for(int n=0; n<time.Size()-1; n++)
    {
      for(int k=0; k<sol.Size(); k++){ sol(k) = ( (fabs(sol(k)-1.0)<1e-7) ) ? 1.0 : sol(k); }
      double dt = time[n] - time[n-1];
      //cout << time[n] << " " << time[n-1] << endl;
      Vector Accum(sol.Size());

      Accumulate(Accum, sol, *F[0], dt, 0);
      for(int k=0; k<sol.Size(); k++)
	{
	  sol(k) -= Accum(k)/porosity(k);
	}

      Vector phiS(sol.Size());
      for(int i=0; i<sol.Size(); i++){ phiS(i) = porosity(i)*sol(i); }

      Accumulate(Accum, phiS, *F[1], dt, 1);
      for(int k=0; k<sol.Size(); k++)
	{
	   sol(k) -= Accum(k)/porosity(k);
	  if(sol(k) < -1.0e-2 || sol(k) > 1.0+1.0e-2 || sol(k) !=sol(k)){ cout << "saturation at node " << k << " is " << sol(k) << endl; test=true; }
	}
      if(test==true){exit(2);}
    }  
  for(int k=0; k<sol.Size(); k++){ sol(k) = ( sol(k)-1.0 > 1e-16 ) ? 1.0 : sol(k); }

}

//=============================================================

void HyperbolicProblem::Accumulate(Vector &Accum, Vector &solold, FData &F, double dt, int b)
{
  /*
    Accumulate the boundary flux effects on each control volume. This function visits each 
    element and adds the flux contribution (on each "sub"-control volume boundary in that element)
    to the accumulation vector. Symbolically, Accum = c*F.
   */
 
  Accum = 0.0;    

  /*
  for(int bdry=0; bdry<dualmesh->GetNBdrs()-1; bdry++)
    {
      for(int i=0; i<=dualmesh->GetNBdryElem(bdry); i++) 
	{
	  int vindex = dualmesh->GetBdryGindex(bdry, i);
	  double sval[1];
	  sval[0] = solold(vindex);
	  Accum(vindex) = fluxfunction[b]->Eval(sval)*F.GetBdryFVal(vindex);
	}
    }
  
  Vector Areas(Accum.Size()); Areas = 0.0;
  int numel = dualmesh->GetNE();
  int locnumv = dualmesh->GetNumLocalDOF();
  Array<int> dof;
  for(int el=0; el<numel; el++)
    {
      dualmesh->GetDOF(el, dof);
      for(int k=0; k<dof.Size(); k++)
	{
	  Accum(dof[k]) += FluxContribution(solold, dof, F, el, k, b);
	  Areas(dof[k]) += dualmesh->GetArea(el, k);
	}
    }

  
    for(int bdry=0; bdry<dualmesh->GetNBdrs(); bdry++)
    {
  
    int bdry = 3;
      for(int i=0; i<=dualmesh->GetNBdryElem(bdry); i++)
	{
	  int vindex = dualmesh->GetBdryGindex(bdry, i);
	  double sval[1];
	  sval[0] = solold(vindex);
	  Accum(vindex) = fluxfunction[b]->Eval(sval)*F.GetBdryFVal(vindex);
	  
	  for(int k=0; k<dualmesh->Getcvnumel(vindex); k++)
	    {
	      int el = dualmesh->Getcvelindex(vindex, k);
	      int locind = dualmesh->Getcvlocalindex(vindex, k);
	      Accum(vindex) += BdrFluxContribution(solold, F,  vindex, el, locind, locnumv, b);
	      Accum(vindex) = 0.0;
	    }
	}
      //  }
  
  for(int k=0; k<Accum.Size(); k++)
    {
      Accum(k) = Accum(k)*dt/((double) Areas(k));
    }
  */
}

//====================================

void HyperbolicProblem::SSPRK2(Vector &initsol, Vector &sol, Vector &porosity, Array<FData *> &F, const Array<double> &time)
{
  Vector K1(sol.Size());
  Vector K2(sol.Size());
  Vector temp(sol.Size());
  temp = initsol;
  //TimeMarch(initsol, K1, porosity, F, time);
  //TimeMarch(K1, K2, porosity, F, time);
  K2 *= 0.5;
  temp *= 0.5;
  add(K2, temp, sol);
}
