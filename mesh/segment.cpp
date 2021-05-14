#include "mesh_header.h"

//============================================================================

Segment::Segment(const int *ind, int attr) : Element()
{
  attribute = attr;
  for (int i=0; i<2; i++) { indices[i] = ind[i]; }
}

//============================================================================

Segment::Segment(const int *ind, int e, int attr) : Element()
{
  attribute = attr;
  for (int i=0; i<2; i++) { indices[i] = ind[i]; }
  edge = e;
}

//============================================================================

void Segment::SetVertices(const int *ind)
{
  indices[0] = ind[0];
  indices[1] = ind[1];
}

//============================================================================

void Segment::GetVertices(Array<int> &v) const
{
  v.SetSize(2);
  for (int i=0; i<2; i++) { v[i] = indices[i]; }
}

//============================================================================

void Segment::GetEdges(Array<int> &e) const
{
  e.SetSize(1);
  e[0] = edge;
}

//============================================================================

Segment::~Segment() {;}

//============================================================================

