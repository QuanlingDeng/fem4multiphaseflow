/**************************************************************************/
/* File:   evalfunc.cc                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 95                                                    */
/**************************************************************************/

/* 
   Function parser
*/

#include <math.h>
#include <iostream>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "array.hpp"
#include "evalfunc.h"

using namespace std;

EvalFunction :: EvalFunction() {}

EvalFunction :: EvalFunction (istream & aist)
{
  ist = &aist;
  ReadNext();
  ParseExpression ();
}

EvalFunction :: ~EvalFunction () {}

void EvalFunction :: AddConstant (double val)
{
  step hstep;
  hstep.op = CONSTANT;
  hstep.operand.val = val;

  program.Append (hstep);
}

void EvalFunction :: AddVariable (int varnum)
{
  step hstep;
  hstep.op = VARIABLE;
  hstep.operand.varnum = varnum;

  program.Append (hstep);
}

void EvalFunction :: AddGlobVariable (const double * aglobvar)
{
  step hstep;
  hstep.op = GLOBVAR;
  hstep.operand.globvar = aglobvar;

  program.Append (hstep);
}

void EvalFunction :: AddFunction (double (*fun) (double))
{
  step hstep;
  hstep.op = FUNCTION;
  hstep.operand.fun = fun;

  program.Append (hstep);
}

void EvalFunction :: AddOperation (EVAL_TOKEN op)
{
  step hstep;
  hstep.op = op;
  hstep.operand.val = 0;

  program.Append (hstep);
}

double EvalFunction :: Eval (const double * x) const
{
  double y;
  Eval (x, &y, 1);
  return y;
}

double EvalFunction :: Eval (Array<double> &x, Array<double> &u) const
{
  double *X = new double[6];
  for (int i=0; i<x.Size(); i++)
    X[i] = x[i];
  for (int i=3; i<6; i++)
    X[i] = u[i-3];

  double Y;
  Eval( (const double*) X, &Y, 1);
  return Y;
}

double EvalFunction :: operator() (double x,double y, double z) const
{
  double *X = new double[3];
  X[0]=x, X[1]=y, X[2]=z;

  double Y;
  Eval( (const double*) X, &Y, 1);
  return Y;
}

void EvalFunction :: Eval (const double * x, double * y, int ydim) const
{
  int i, stacksize;
  Array<double> stack;

  stack.SetSize (program.Size());

  stacksize = 0;

  for (i = 0; i < program.Size(); i++)
    {
      switch (program[i].op)
	{
	case ADD:
	  stacksize--;
	  stack[stacksize-1] += stack[stacksize];
	  break;

	case SUB:
	  stacksize--;
	  stack[stacksize-1] -= stack[stacksize];
	  break;

	case MULT:
	  stacksize--;
	  stack[stacksize-1] *= stack[stacksize];
	  break;

	case DIV:
	  stacksize--;
	  stack[stacksize-1] /= stack[stacksize];
	  break;

	case NEG:
	  stack[stacksize-1] = -stack[stacksize-1];
	  break;

	case CONSTANT:
	  stack[stacksize] = program[i].operand.val;
	  stacksize++;
	  break;

	case VARIABLE:
	  stack[stacksize] = x[program[i].operand.varnum-1];
	  stacksize++;
	  break;

	case GLOBVAR:
	  stack[stacksize] = *program[i].operand.globvar;
	  stacksize++;
	  break;

	case FUNCTION:
	  stack[stacksize-1] = (*program[i].operand.fun) (stack[stacksize-1]);
	  break;
        
	default:
	  cerr << "undefined operation for EvalFunction" << endl;
	}
    }

  if (stacksize != ydim)
    {
      cout << stacksize << endl;
      cerr << "final stacksize not matching ydim" << endl;
      return;
    }

  for (i = 0; i < ydim; i++)
    y[i] = stack[i];
}

int EvalFunction :: IsConstant () const
{
  int i;
  for (i = 0; i < program.Size(); i++)
    {
      EVAL_TOKEN op = program[i].op;

      if (op == VARIABLE || op == GLOBVAR || op == COEFF_FUNC)
	return 0;
    }
  return 1;
}

void EvalFunction :: Print (ostream & ost) const
{
  int i;
  for (i = 0; i < program.Size(); i++)
    {
      EVAL_TOKEN op = program[i].op;
      ost << "Step " << i << ": " << (int)op << " = " << (char) op 
	  << ", val = " << program[i].operand.val << endl;
    }
}

void EvalFunction :: ParseExpression ()
{
  ParseTerm ();

  while (1)
    {
      switch (GetToken())
	{
	case ADD:
	  {
	    ReadNext();
	    ParseTerm ();
	    AddOperation (ADD);
	    break;
	  }
	case SUB:
	  {
	    ReadNext();
	    ParseTerm ();
	    AddOperation (SUB);
	    break;
	  }
	default:
	  return;
	}
    }
}

void EvalFunction :: ParseTerm ()
{
  ParsePrimary();
  
  while (1)
    {
      switch (GetToken())
	{
	case MULT:
	  {
	    ReadNext();
	    ParsePrimary();
	    AddOperation (MULT);
	    break;
	  }
	case DIV:
	  {
	    ReadNext();
	    ParsePrimary();
	    AddOperation (DIV);
	    break;
	  }
	default:
	  return;
	}
    }
}

void EvalFunction :: ParsePrimary()
{
  switch (GetToken())
    {
    case CONSTANT:
      {
	ReadNext();
	AddConstant (GetNumValue());
	break;
      }
    case SUB:
      {
	ReadNext();
	ParsePrimary();
	AddConstant (-1);
	AddOperation (MULT);
	break;
      }
    case LP:
      {
	ReadNext();
	ParseExpression();
	ReadNext();
	break;
      }
    case FUNCTION:
      {
	ReadNext();
	double (*fun)(double) = GetFunctionValue();
	ParsePrimary();
	AddFunction (fun);
	break;
      }
    case VARIABLE:
      {
	ReadNext();
	AddVariable (GetVariableNumber());
	break;
      }
    }
}

void EvalFunction :: ReadNext ()
{
  char ch;
  
  do
    { // skip whitespaces
      (*ist).get(ch);
      if ((*ist).eof() || ch == ';')
	{
	  token = END;
	  return;
	}
    }
  while (isspace(ch));
  
  
  switch (ch)
    {
    case '*': case '/': 
    case '+': case '-': 
    case '(': case ')':
      {
	token = EVAL_TOKEN (ch);
	break;
      }
      
    default:
      {
	if (isdigit (ch) || ch == '.')
	  {
	    (*ist).putback (ch);
 	    (*ist) >> num_value;
	    token = CONSTANT;
	  }
	else
	  {
	    int strcnt = 0;
	    while (isalnum(ch))
	      {
		string_value[strcnt++] = ch;
		(*ist).get(ch);
	      }
	    (*ist).putback (ch);
	    string_value[strcnt]='\0';
	    
	    if (strcmp (string_value, "x") == 0)
	      {
		var_num = 1;
		token = VARIABLE;
	      }
	    else if (strcmp (string_value, "y") == 0)
	      {
		var_num = 2;
		token = VARIABLE;
	      }
	    else if (strcmp (string_value, "z") == 0)
	      {
		var_num = 3;
		token = VARIABLE;
	      }
	    else if (strcmp (string_value, "u") == 0)
	      {
		var_num = 4;
		token = VARIABLE;
	      }
	    else if (strcmp (string_value, "v") == 0)
	      {
		var_num = 5;
		token = VARIABLE;
	      }
	    else if (strcmp (string_value, "w") == 0)
	      {
		var_num = 6;
		token = VARIABLE;
	      }
	    else if (strcmp (string_value, "t") == 0)
	      {
		var_num = 7;
		token = VARIABLE;
	      }
	    else if (strcmp (string_value, "abs") == 0)
	      {
		fun_value = fabs;
		token = FUNCTION;
	      }
	    else if (strcmp (string_value, "sin") == 0)
	      {
		fun_value = sin;
		token = FUNCTION;
	      }
	    else if (strcmp (string_value, "cos") == 0)
	      {
		fun_value = cos;
		token = FUNCTION;
	      }
	    else if (strcmp (string_value, "tan") == 0)
	      {
		fun_value = tan;
		token = FUNCTION;
	      }
	    else if (strcmp (string_value, "atan") == 0)
	      {
		fun_value = atan;
		token = FUNCTION;
	      }
	    else if (strcmp (string_value, "exp") == 0)
	      {
		fun_value = exp;
		token = FUNCTION;
	      }
	    else if (strcmp (string_value, "log") == 0)
	      {
		fun_value = log;
		token = FUNCTION;
	      }
	    else if (strcmp (string_value, "sqrt") == 0)
	      {
		fun_value = sqrt;
		token = FUNCTION;
	      }
	    else if (strcmp (string_value, "sinh") == 0)
	      {
		fun_value = sinh;
		token = FUNCTION;
	      }
	    else if (strcmp (string_value, "cosh") == 0)
	      {
		fun_value = cosh;
		token = FUNCTION;
	      }
	    else if (strcmp (string_value, "tanh") == 0)
	      {
		fun_value = tanh;
		token = FUNCTION;
	      }
	    else if (strcmp (string_value, "pi") == 0)
	      {
		num_value = M_PI;
		token = CONSTANT;
	      }
	    else if (strcmp (string_value, "e") == 0)
	      {
		num_value = M_E;
		token = CONSTANT;
	      }
	    else
	      token = STRING;
	  }
      }
    }
}
