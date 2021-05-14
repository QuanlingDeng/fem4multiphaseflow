#include "fem_header.h"

Data::Data()
{
  EllipticExists = false;
  AdvectionExists = false;
  ConservativeAdvectionExists = false;
  ReactionExists = false;
  ForceExists = false;
  NeumannExists = false;
  RobinExists = false;
  StableExists = false;
  StableRHSExists = false;
  NeumannBdrExists = false;
  RobinBdrExists = false;
  DirichletBdrExists = false;
  PointSourceForceExists = false;
  nbdryconds = 1;
}

//==============================================================================

void Data::SetEllipticData(const Array<double *> &_ellipticcoeff)
{
  EllipticExists = true;
  EllipticCoeff.SetSize(_ellipticcoeff.Size());
  for(int i=0; i<_ellipticcoeff.Size(); i++) { EllipticCoeff[i] = _ellipticcoeff[i]; }
}

void Data::SetEllipticFunction(Function *_elliptic)
{
  EllipticExists = true;
  elliptic = _elliptic;
}

//==============================================================================

void Data::SetForceData(const Array<double *> &_forcecoeff)
{
  ForceExists = true;
  ForceCoeff.SetSize(_forcecoeff.Size());
  for(int i=0; i<_forcecoeff.Size(); i++) { ForceCoeff[i] = _forcecoeff[i]; }
}

void Data::SetForceFunction(Function *_force)
{
  ForceExists = true;
  force = _force;
}

//==============================================================================

void Data::SetPointSourceData(const Array<double> &_fvalue, const Array<int> &_nodeindex)
{
  PointSourceForceExists = true;
  PointSourceValues.SetSize(_fvalue.Size());
  for(int i=0; i<PointSourceValues.Size(); i++) { PointSourceValues[i] = _fvalue[i]; }

  PointSourceIndex.SetSize(_nodeindex.Size());
  for(int i=0; i<PointSourceIndex.Size(); i++) { PointSourceIndex[i] = _nodeindex[i]; }
}

//==============================================================================

void Data::SetAdvectionData(const Array<double *> &_advectioncoeff)
{
  AdvectionExists = true;
  AdvectionCoeff.SetSize(_advectioncoeff.Size());
  for(int i=0; i<_advectioncoeff.Size(); i++) { AdvectionCoeff[i] = _advectioncoeff[i]; }
}

void Data::SetAdvectionFunction(Function *_advection)
{
  AdvectionExists = true;
  advection = _advection;
}

//==============================================================================

void Data::SetConservativeAdvectionData(const Array<double *> &_conservativeadvectioncoeff)
{
  ConservativeAdvectionExists = true;
  ConservativeAdvectionCoeff.SetSize(_conservativeadvectioncoeff.Size());
  for(int i=0; i<_conservativeadvectioncoeff.Size(); i++) { ConservativeAdvectionCoeff[i] = _conservativeadvectioncoeff[i]; }
}

void Data::SetConservativeAdvectionFunction(Function *_conservativeadvection)
{
  ConservativeAdvectionExists = true;
  conservativeadvection = _conservativeadvection;
}

//==============================================================================

void Data::SetReactionData(const Array<double *> &_reactioncoeff)
{
  ReactionExists = true;
  ReactionCoeff.SetSize(_reactioncoeff.Size());
  for(int i=0; i<_reactioncoeff.Size(); i++) { ReactionCoeff[i] = _reactioncoeff[i]; }
}

void Data::SetReactionFunction(Function *_reaction)
{
  ReactionExists = true;
  reaction = _reaction;
}

//==============================================================================

void Data::SetNeumannData(const Array<double *> &_neumanncoeff)
{
  NeumannExists = true;
 
  BdrNeumann.SetSize(1);
  BdrNeumann[0] = true;

  NeumannCoeff.SetSize(_neumanncoeff.Size());
  for(int i=0; i<_neumanncoeff.Size(); i++) { NeumannCoeff[i] = _neumanncoeff[i]; }
}

void Data::SetNeumannFunction(Function *_neumann)
{
  NeumannExists = true;
  neumann = _neumann;
}

//==============================================================================

void Data::SetRobinData(const Array<double *> &_robincoeff)
{
  RobinExists = true;
  
  BdrRobin.SetSize(1);
  BdrRobin[0] = true;

  RobinCoeff.SetSize(_robincoeff.Size());
  for(int i=0; i<_robincoeff.Size(); i++) { RobinCoeff[i] = _robincoeff[i]; }
}

void Data::SetRobinFunction(Function *_robin)
{
  RobinExists = true;
  robin = _robin;
}

//==============================================================================

void Data::SetNeumannData(const Array<double *> &_NeumannBdryVal, const Array<bool> &_BdrNeumann)
{
  NeumannBdrExists = true;

  BdrNeumann.SetSize(_BdrNeumann.Size());
  for(int i=0; i<_BdrNeumann.Size(); i++) { BdrNeumann[i] = _BdrNeumann[i]; }
  
  NeumannBdryVal.SetSize(_NeumannBdryVal.Size());
  for(int i=0; i<_NeumannBdryVal.Size(); i++) { NeumannBdryVal[i] = _NeumannBdryVal[i]; }
}

//==============================================================================

void Data::SetRobinData(const Array<double *> &_RobinCoeff, const Array<double *> &_RobinBdryVal, const Array<bool> &_BdrRobin)
{
  RobinBdrExists = true;
 
  BdrRobin.SetSize(_BdrRobin.Size());
  for(int i=0; i<_BdrRobin.Size(); i++) { BdrRobin[i] = _BdrRobin[i]; }

  RobinCoeff.SetSize(_RobinCoeff.Size());
  for(int i=0; i<_RobinCoeff.Size(); i++) { RobinCoeff[i] = _RobinCoeff[i]; }

  RobinBdryVal.SetSize(_RobinBdryVal.Size());
  for(int i=0; i<_RobinBdryVal.Size(); i++) { RobinBdryVal[i] = _RobinBdryVal[i]; }
}

//==============================================================================

void Data::SetDirichletData(const Array<double *> &_DirichletBdryVal, const Array<bool> &_BdrDirichlet)
{
  DirichletBdrExists = true;

  BdrDirichlet.SetSize(_BdrDirichlet.Size());
  for(int i=0; i<_BdrDirichlet.Size(); i++) { BdrDirichlet[i] = _BdrDirichlet[i]; }
  
  DirichletBdryVal.SetSize(_DirichletBdryVal.Size());
  for(int i=0; i<_DirichletBdryVal.Size(); i++) { DirichletBdryVal[i] = _DirichletBdryVal[i]; }
}

void Data::SetDirichletFunction(Function *_dirichlet)
{
  DirichletBdrExists = true;
  dirichlet = _dirichlet;
}

//==============================================================================

void Data::SetStableData(const Array<double *> &_stablecoeff)
{
  StableExists = true;
  StableCoeff.SetSize(_stablecoeff.Size());
  for(int i=0; i<_stablecoeff.Size(); i++) { StableCoeff[i] = _stablecoeff[i]; }
}

void Data::SetStableFunction(Function *_stable)
{
  StableExists = true;
  stable = _stable;
}

//==============================================================================

void Data::SetStableRHSData()
{
  StableRHSExists = true;
}

void Data::SetStableRHSFunction(Function *_stablerhs)
{
  StableRHSExists = true;
  stablerhs = _stablerhs;
}

//==============================================================================


Data::~Data()
{
  ;
}

