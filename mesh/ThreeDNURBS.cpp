#include "mesh_header.h"

ThreeDNURBS::ThreeDNURBS(const int *ind, int attr) : Element()
{
  attribute = attr;
  for(int i=0; i<4; i++)
    indices[i] = ind[i];
}

void ThreeDNURBS::SetVertices(const int *ind)
{
  indices[0] = ind[0];
  indices[1] = ind[1];
  indices[2] = ind[2];
  indices[3] = ind[3];
}

void ThreeDNURBS::GetVertices(Array<int> &v) const
{
  v.SetSize(4);
  for( int i=0; i<4; i++)
    {
      v[i] = indices[i];
    }
}

ThreeDNURBS:: ~ThreeDNURBS() { ;}
