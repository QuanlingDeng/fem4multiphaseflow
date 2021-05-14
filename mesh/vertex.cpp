#include "mesh_header.h"

//============================================================================

template <int Dim>
Vertex<Dim>::Vertex(double x)
{
  coord[0] = x;
}

//============================================================================

template <int Dim>
Vertex<Dim>::Vertex(double x, double y)
{
  coord[0] = x;
  coord[1] = y;
}

//============================================================================

template <int Dim>
Vertex<Dim>::Vertex(double x, double y, double z)
{
  coord[0] = x;
  coord[1] = y;
  coord[2] = z;
}

//============================================================================

template <int Dim>
Vertex<Dim>::~Vertex() {;}

//============================================================================

