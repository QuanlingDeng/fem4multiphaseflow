#include "geometry_header.h"
#include "../general/array.hpp"

void Container::RD(Array<Vector *> &input)
{
  int sz = input.Size();
  Array<Vector* > temp(sz);
  for (int i=0; i<sz; i++)
    {
      int s = input[i]->Size();
      temp[i] = new Vector(s);
      for (int j=0; j<s; j++)
	{
	  double d = input[i]->Elem(j);
	  temp[i]->Elem(j) = d;
	}
      data[i]->SetSize(s);
      for (int j=0; j<s; j++)
	{
	  double d = temp[i]->Elem(j);
	  data[i]->Elem(j) = d;
	}
    }
  for (int i=0; i<sz; i++)
    {
      delete temp[i];
    } 
}

//============================================================================

Container::Container(ifstream *file, int m)
{
  data.SetSize(m);
  for (int i=0; i<m; i++)
    {
      int temp;
      *file >> temp;
      data[i] = new Vector(temp);
      double temp1;
      for (int j=0; j<temp; j++)
	{
	  *file >> temp1;
	  data[i]->Elem(j) = temp1;
	}
    }
}

//============================================================================

Container::Container(Array<Vector *> &input, int altsize)
{
  data.SetSize(input.Size());
  for (int i=0; i<input.Size(); i++)
    {
      int temp = altsize;
      data[i] = new Vector(temp);
      double temp1;
      for (int j=0; j<input[i]->Size(); j++)
	{
	  temp1 = input[i]->Elem(j);
	  data[i]->Elem(j) = temp1;
	}
    }
}

//============================================================================

Container::~Container()
{
  for (int i=0; i<data.Size(); i++)
    {
      delete data[i];
    }
}
