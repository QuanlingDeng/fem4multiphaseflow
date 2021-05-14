#include "fem_header.h"

BdrData::BdrData()
{

}

void BdrData::SetNeumannData(const Array<double> &_RHS)
{
  RHS.SetSize(_RHS.Size());
  for(int i=0; i<_RHS.Size(); i++) { RHS[i] = _RHS[i]; }

  ReactionCoeff.SetSize(_RHS.Size());
  for(int i=0; i<_RHS.Size(); i++) { ReactionCoeff[i] = 0.0; }
}

//==============================================================================

void BdrData::SetRobinData(const Array<double> &_RHS, const Array<double> &_ReactionCoeff)
{
  RHS.SetSize(_RHS.Size());
  for(int i=0; i<_RHS.Size(); i++) { RHS[i] = _RHS[i]; }

  ReactionCoeff.SetSize(_ReactionCoeff.Size());
  for(int i=0; i<_ReactionCoeff.Size(); i++) { ReactionCoeff[i] = _ReactionCoeff[i]; }
}

//==============================================================================


BdrData::~BdrData()
{
  ;
}

