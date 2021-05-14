#include "mesh_header.h"

Triangle::Triangle(const int *ind, int attr) : Element()
{
  attribute = attr;
  for(int i=0; i<3; i++)
    indices[i] = ind[i];
}

Triangle::Triangle(const int *ind, int *edge, int attr) : Element()
{
  attribute = attr;
  for(int i=0; i<3; i++)
    {
      indices[i] = ind[i];
      edges[i] = edge[i];
    }
}


void Triangle::SetVertices(const int *ind)
{
  indices[0] = ind[0];
  indices[1] = ind[1];
  indices[2] = ind[2];
}

void Triangle::GetVertices(Array<int> &v) const
{
  v.SetSize(3);
  for(int i=0; i<3; i++)
    {
      v[i] = indices[i];
    }
}

void Triangle::GetEdges(Array<int> &e) const
{
  e.SetSize(3);
  for(int i=0; i<3; i++)
    e[i] = edges[i];
}

Triangle::~Triangle(){;}
