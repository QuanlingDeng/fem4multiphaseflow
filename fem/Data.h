#ifndef FILE_DATA
#define FILE_DATA

/* This is a container for nodal values of data such as elliptic coefficients,
   forcing, linear reaction. The nodal values are constructed (outside this class)
   to follow the values of a function at points associated with vertices of mesh
   (and nodes associated with basis for higher order polynomials) */
class Data
{
protected:

  Array<double *> EllipticCoeff;
  Array<double *> ForceCoeff;
  Array<double *> AdvectionCoeff;
  Array<double *> ConservativeAdvectionCoeff;
  Array<double *> ReactionCoeff;
  Array<double *> NeumannCoeff;
  Array<double *> StableCoeff;

  Array<double *> RobinCoeff;
  Array<double *> RobinBdryVal;
  Array<double *> NeumannBdryVal;
  Array<double *> DirichletBdryVal;

  Array<double> PointSourceValues;
  Array<int> PointSourceIndex;

  Array<bool> BdrNeumann;
  Array<bool> BdrRobin;
  Array<bool> BdrDirichlet;

  bool EllipticExists;
  bool AdvectionExists;
  bool ConservativeAdvectionExists;
  bool ReactionExists;
  bool ForceExists;
  bool NeumannExists;
  bool RobinExists;
  bool NeumannBdrExists;
  bool RobinBdrExists;
  bool DirichletBdrExists;
  bool StableExists;
  bool StableRHSExists;
  bool PointSourceForceExists;
  int nbdryconds;

  Function *elliptic;
  Function *force;
  Function *advection;
  Function *conservativeadvection;
  Function *reaction;
  Function *neumann;
  Function *robin;
  Function *dirichlet;
  Function *stable;
  Function *stablerhs;
  
public:
        
  Data();

  /// Set the nodal values of elliptic coefficient
  void SetEllipticData(const Array<double *> &_ellipticcoeff);
  void SetEllipticFunction(Function *_elliptic);

  /// Set the nodal values of force coefficient
  void SetForceData(const Array<double *> &_forcecoeff);
  void SetForceFunction(Function *_force);

  /// Set the nodal values of advection coefficient
  void SetAdvectionData(const Array<double *> &_advectioncoeff);
  void SetAdvectionFunction(Function *_advection);

  /// Set the nodal values of conservative advection coefficient
  void SetConservativeAdvectionData(const Array<double *> &_conservativeadvectioncoeff);
  void SetConservativeAdvectionFunction(Function *_conservativeadvection);

  /// Set the nodal values of reaction coefficient
  void SetReactionData(const Array<double *> &_reactioncoeff);
  void SetReactionFunction(Function *_reaction);

  void SetPointSourceData(const Array<double> &_fvalue, const Array<int> &_nodeindex);

  /// Set nodal Neumann values
  void SetNeumannData(const Array<double *> &_neumanncoeff);
  void SetNeumannFunction(Function *_neumann);

  /// Set nodal Robin values
  void SetRobinData(const Array<double *> &_robincoeff);
  void SetRobinFunction(Function *_robin);

  /// Set Bdry Neumann values
  void SetNeumannData(const Array<double *> &_NeumannBdryVal, const Array<bool> &_BdrNeumann);

  /// Set Bdry Robin values
  void SetRobinData(const Array<double *> &_RobinCoeff, const Array<double *> &_RobinBdryVal, const Array<bool> &_BdrRobin);

  /// Set Bdry Dirichlet values
  void SetDirichletData(const Array<double *> &_DirichletBdryVal, const Array<bool> &_BdrDirichlet);
  void SetDirichletFunction(Function *_dirichlet);

  /// Set elemental stable parameter values
  void SetStableData(const Array<double *> &_stablecoeff);
  void SetStableFunction(Function *_stable);

  /// Set StableRHSExists = true
  void SetStableRHSData( );
  void SetStableRHSFunction(Function *_stablerhs);

  /// Set Number of Boundary Conditions
  void SetNumberofBoundaryConditions(int k){ nbdryconds = k;}

  /// Get the i-th nodal value of elliptic coefficient
  inline double GetNodalEllipticCoeff(int i, int j=0) { return EllipticCoeff[j][i]; }
  inline Function *GetEllipticFunction() { return elliptic; }

  /// Get the i-th nodal value of force coefficient
  inline double GetNodalForceCoeff(int i, int j=0) { return ForceCoeff[j][i]; }
  inline Function *GetForceFunction() { return force; }

  /// Get the i-th nodal value of advection coefficient
  inline double GetNodalAdvectionCoeff(int i, int j=0) { return AdvectionCoeff[j][i]; }
  inline Function *GetAdvectionFunction() { return advection; }

  /// Get the i-th nodal value of conservative advection coefficient
  inline double GetNodalConservativeAdvectionCoeff(int i, int j=0) { return ConservativeAdvectionCoeff[j][i]; }
  inline Function *GetConservativeAdvectionFunction() { return conservativeadvection; }

  /// Get the i-th nodal value of reaction coefficient
  inline double GetNodalReactionCoeff(int i, int j=0) { return ReactionCoeff[j][i]; }
  inline Function *GetReactionFunction() { return reaction; }

  /// Get the Neumann Boundary Value
  inline double GetNodalNeumannCoeff(int i, int j=0) { return NeumannCoeff[j][i]; }
  inline Function *GetNeumannFunction() { return neumann; }

  /// Get the Robin Boundary Value
  inline double GetNodalRobinCoeff(int i, int j=0) { return RobinCoeff[j][i]; }
  inline Function *GetRobinFunction() { return robin; }

  /// Get the Neumann Boundary Value
  inline double GetNeumannBdryVal(int b, int i=0) { return NeumannBdryVal[b][i]; }

  /// Get the Robin Boundary Value
  inline double GetRobinBdryVal(int b, int i=0) { return RobinBdryVal[b][i]; }

  /// Get the Robin Coeff Value
  inline double GetRobinCoeff(int b, int i=0) { return RobinCoeff[b][i]; }

  /// Get the Dirichlet Boundary Value
  inline double GetDirichletBdryVal(int b, int i=0) { return DirichletBdryVal[b][i]; }
  inline Function *GetDirichletFunction() { return dirichlet; }

  /// Get the i-th elemental value of stable coefficient
  inline double GetElementalStableCoeff(int i, int j=0) { return StableCoeff[j][i]; }
  inline Function *GetStableFunction() { return stable; }

  /// Get the PointSource Data 
  inline double GetPointSourceValue(int i) { return PointSourceValues[i]; } 

  inline int GetPointSourceIndex(int i) { return PointSourceIndex[i]; } 

  inline int GetNumberofPointSources(){ return PointSourceIndex.Size();}

  inline int GetNumberofBoundaryConditions(){ return nbdryconds;}

  inline bool EllipticPartExists() { return EllipticExists; }

  inline bool AdvectionPartExists() { return AdvectionExists; }

  inline bool ConservativeAdvectionPartExists() { return ConservativeAdvectionExists; }

  inline bool ReactionPartExists() { return ReactionExists; }

  inline bool ForcePartExists() { return ForceExists; }

  inline bool NeumannPartExists() { return NeumannExists; }

  inline bool RobinPartExists() { return RobinExists; }

  inline bool BdrNeumannExists() { return NeumannBdrExists; }

  inline bool BdrRobinExists() { return RobinBdrExists; }

  inline bool BdrDirichletExists() { return DirichletBdrExists; }

  inline bool StablePartExists() { return StableExists; }

  inline bool StableRHSPartExists() { return StableRHSExists; }

  inline bool PointSourcePartExists() { return PointSourceForceExists; }

  /// Is the jth boundary Neumann?
  inline bool BdryNeumann(int j=0) { return BdrNeumann[j]; }

  /// Is the jth boundary Robin?
  inline bool BdryRobin(int j=0) { return BdrRobin[j]; }

  /// Is the jth boundary Dirichlet?
  inline bool BdryDirichlet(int j=0) { return BdrDirichlet[j]; }


  ~Data();
};

#endif
