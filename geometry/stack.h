#ifndef FILE_STACK
#define FILE_STACK
#include "geometry_header.h"

class Stack {
 protected:
  Array<Container *> stack;
  
 public:
 
  Stack(Patch & patch, int pdim, int sdim);
  
  virtual ~Stack();
};
#endif
