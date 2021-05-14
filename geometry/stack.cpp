#include "geometry_header.h"

Stack::Stack(Patch &patch, int pdim, int sdim)
{
  if (sdim==2)
    {stack.SetSize(1);}
  else
    {stack.SetSize(patch.GetKnotSize(2)+patch.GetDegs(2)-2);}
}

//============================================================================

Stack::~Stack()
{
  for (int i=0; i<stack.Size(); i++){delete stack[i];}
}
