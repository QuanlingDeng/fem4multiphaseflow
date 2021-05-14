#ifndef FILE_BDRDATA
#define FILE_BDRDATA

/*
   Data type for boundary data

   Author: Quanling Deng
   Date: 09/04/2015
*/

/* This is a container for boundary values of data such as neumann, robin, and dirichlet. The nodal values are constructed (outside this class)
   to follow the values of a function at points associated with vertices of mesh boundarys
   (and nodes associated with basis for higher order polynomials) */

class BdrData
{
protected:

  Array<double> ReactionCoeff;
  Array<double> RHS;
  
public:
        
  BdrData();

  /// Set nodal Robin values
  void SetNeumannData(const Array<double> &_RHS);

  /// Set elemental stable parameter values
  void SetRobinData(const Array<double> &_RHS, const Array<double> &_ReactionCoeff);

  /// Get the i-th nodal value of robin coefficient
  inline double GetNodalRHSVal(int i) { return RHS[i]; }

  /// Get the i-th elemental value of stable coefficient
  inline double GetNodalReactionCoeff(int i) { return ReactionCoeff[i]; }

  ~BdrData();
};

#endif
