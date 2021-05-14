#ifndef FILE_FUNCTION
#define FILE_FUNCTION

/// Define base class function
using namespace std;

class Function {

 public:

  virtual double Eval(double *y) = 0;

  virtual void Read(istream &in) = 0;

  virtual ~Function() {;};
};

/// Subclass constant coefficient.
class ConstantFunction : public Function {

 public:
  
  /// Constructor with constant c.
  ConstantFunction(double c=1.0) { constant = c; };

  /// Evaluate the coefficient.
  virtual double Eval(double *y);

  /// Read constant from file.
  virtual void Read(istream &in);

 private:
 
  double constant;
};

/// Subclass for functions given by string.
class SymbolFunctionFunction : public Function {

 public:

  /// Constructor.
  SymbolFunctionFunction(istream &in) : Function(in) { };

  /// Evaluate the coefficient.
  virtual double Eval(double *y);

  virtual void Read(istream &in) { };

 private:
  EvalFunction Function;
};

int ReadIdentifier (istream & input, char *buffer, int buflen);

Function * ReadFunction (istream & input);


/// Vector of Functions.
class VectorFunction 
{
 private:
  Array<Function*> Func;
  int vdim;

 public:
  /// Construct vector of dim coefficients.
  VectorFunction(int dim) { vdim = dim; Func.SetSize(dim); };

  /// Returns dimension of the vector.
  int GetVDim() { return vdim; };
  
  /// Returns i'th coefficient.
  Function & GetFunc(int i) { return *Func[i]; };

  /// Sets coefficient in the vector.
  void Set(int i, Function * c) { Func[i] = c; };

  /// Evaluate the i'th component of vector coefficient.
  double Eval(int i, double *y) { return Func[i]->Eval(y); };

  /// Reads the i'th component of vector coefficient.
  void Read(int i, istream &in) { Func[i]->Read(in); };
};

class MatrixFunction
{
 private:
  Array<Function *> Func;
  int vdim;
  
 public:
  
  MatrixFunction(int dim) { vdim = dim; Func.SetSize(dim*dim); };
  
  int GetVDim() { return vdim; };
  
  /// Return the (i,j)'th coefficient.
  Function & GetFunc( int i, int j ) { return *Func[i*vdim+j]; };
  
  /// Set the coefficient in the matrix.
  void Set( int i, int j, Function * c ) { Func[i*vdim+j] = c; };

  /// Evaluate the (i,j)'th component of the matrix coefficient.
  double Eval(int i, int j, double *y) { return Func[i*vdim+j]->Eval(y); };

  /// Read the (i,j)'th component of the matrix coefficient.
  void Read(int i, int j, istream &in)  { Func[i*vdim+j]->Read(in); };
};

#endif
