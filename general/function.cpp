#include <math.h>
#include <iostream>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "array.hpp"
#include "evalfunc.h"
#include "function.h"


using namespace std;

double ConstantFunction::Eval(double *y)
{
  return constant;
}

void ConstantFunction::Read(istream &in)
{
  in >> constant;
}

double SymbolFunctionFunction::Eval(double *y)
{
  return Function.Eval(y);
}


int ReadIdentifier (istream & input, char *buffer, int buflen)
{
   int pos = 0;

   input >> ws;
   while ( input && !input.eof() &&
           isalpha(input.peek()) && pos < buflen )
      buffer[pos++] = input.get();

   if (pos == buflen)
      pos--;

   buffer[pos] = '\0';

   return pos;
}

Function * ReadFunction (istream & input)
{
   const int buflen = 256;
   char buffer[buflen];
   Function * Func;

   ReadIdentifier (input, buffer, buflen);
   input >> ws;

   if (!strcmp(buffer, "const"))
      Func = new ConstantFunction;
   else if (!strcmp(buffer, "function"))
      Func = new SymbolFunctionFunction(input);
   else
      Func = NULL;

   if (Func)
      Func -> Read(input);

   return Func;
}
