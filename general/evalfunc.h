#ifndef FILE_EVALFUNC
#define FILE_EVALFUNC

/**************************************************************************/
/* File:   evalfunc.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 95                                                    */
/**************************************************************************/

/**
  Generic Function
  function is saved in reverse Polnish notation
*/

using namespace std;

class EvalScanner
{
  istream * scanin;
};

/// Parser for math function given by a string.
class EvalFunction
{
public:
  enum EVAL_TOKEN
  {
    ADD = '+', SUB = '-', MULT = '*', DIV = '/', LP ='(', RP = ')',
    AND = 100, OR, NOT, LESS, GREATER, LESSEQ, GREATEREQ, EQ, NEQ,
    COND,
    NEG, CONSTANT, VARIABLE, FUNCTION, GLOBVAR,
    COEFF_FUNC, END, STRING
  };

  ///
  EvalFunction ();
  /// Parse from input stream.
  EvalFunction (istream & aist);
  ///
  virtual ~EvalFunction ();
  ///
  void AddConstant (double val);
  
  ///
  void AddVariable (int varnum);
  ///
  void AddGlobVariable (const double * dp);
  ///
  void AddOperation (EVAL_TOKEN op);
  ///
  void AddFunction (double (*fun) (double));
  ///
  double Eval (const double * x) const;
  ///
  double operator() (double x, double y, double z) const;
  ///
  void Eval (const double * x, double * y, int ydim) const;
  ///
  double Eval (Array<double> &x, Array<double> &u) const;
  ///
  int IsConstant () const;
  ///
  void Print (ostream & ost) const;

protected:  
  ///
  class step
  {
  public:
    ///
    EVAL_TOKEN op;
    ///
    union
    {
      double val;
      const double *globvar;
      int varnum;
      double (*fun) (double);
    } 
    operand;
  };
  ///
  Array<step> program;

  // Parsing routines:
  void ParseExpression ();
  void ParseTerm ();
  void ParsePrimary ();

  istream * ist;
  EVAL_TOKEN token;
  double num_value;
  char string_value[1000];
  double (*fun_value)(double);
  char var_num;

  EVAL_TOKEN GetToken() const
    { return token; }

  double GetNumValue() const
    { return num_value; }

  int GetVariableNumber() const
    { return var_num; }

  const char * GetStringValue() const
    { return string_value; }
  
  double (*GetFunctionValue())(double)
    { return fun_value; }
  
  void ReadNext();
};

#endif
