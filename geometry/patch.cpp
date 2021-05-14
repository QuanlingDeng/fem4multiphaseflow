#include "geometry_header.h"

Patch::Patch(ifstream *file, int pdim, int sdim)
{
  patch.SetSize(4);
  patch[0] = new Container(file, pdim);
  patch[1] = new Container(file, sdim);
  patch[2] = new Container(file,1);
  patch[3] = new Container(file, pdim);
  int temp;
  *file >> temp;
  face->SetSize(temp);
  for (int i=0; i<face->Size(); i++)
    {
      *file >> temp;
      face->Elem(i) = temp;
    }
}

//============================================================================

void Patch::Refine(Array<Vector *> &knots, Array<Vector *> &ctrlpts)
{
  patch[0]->RD(knots);
  Array<Vector *> tpts(ctrlpts[0]->Size()-1);
  Array<Vector *> twts(1);
  twts[0] = new Vector(ctrlpts.Size());
  for (int i=0; i<ctrlpts[0]->Size()-1; i++)
    {
      tpts[i] = new Vector(ctrlpts.Size());
      for (int j=0; j<ctrlpts.Size(); j++)
	{
	  double a = ctrlpts[j]->Elem(i);
	  tpts[i]->Elem(j) = a;
	}
    }
  for (int i=0; i<ctrlpts.Size(); i++)
    {
      twts[0]->Elem(i) = ctrlpts[i]->Elem(ctrlpts[0]->Size()-1);
    }
  patch[1]->RD(tpts);
  patch[2]->RD(twts);
  for (int i=0; i<ctrlpts[0]->Size()-1; i++)
    {
      delete tpts[i];
    }
  delete twts[0];
} 

//============================================================================

Patch::~Patch()
{
  for (int i=0; i<4; i++){delete patch[i];}
  delete face;
}
