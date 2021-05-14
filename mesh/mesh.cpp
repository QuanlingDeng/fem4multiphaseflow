#include "mesh_header.h"

int Mesh::GetVertToEle(int i, int j)
{
  return vert_to_ele[i][j];
}

//============================================================================

Mesh::~Mesh()
{
  for (int i=0; i<NumOfElements; i++) { delete element[i]; }
 
  for (int i=0; i<NumOfBdrElements; i++) { delete bdrelement[i]; }

  for(int i=0; i<NumOfVertices; i++)
    delete []vert_to_ele[i];

  // delete []Degs;
}

//============================================================================

void Mesh::Print(ofstream &out, int type) const
{
  if (type==Element::TRIANGLE)
    {
      out << "areamesh2\n\n" << NumOfBdrElements << endl;
      for(int j=0; j<NumOfBdrElements; j++)
	{
	  int nv = GetBdrElementNVertices(j);
	  out << GetBdrAttribute(j) << "     ";
	  Array<int> ind;
	  GetBdrElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    out << ind[i]+1 << "   ";
	  out << endl;
	}

      out << NumOfElements << endl;
      for(int j=0; j<NumOfElements; j++)
	{
	  int nv = GetElementNVertices(j);
	  out << GetAttribute(j) << "  ";
	  out << nv  << "   ";
	  Array<int> ind;
	  GetElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    out << ind[i]+1 << "  ";
	  out << endl;
	}
      
      // Create vertices.
      out << NumOfVertices << endl;
      for(int j=0; j<NumOfVertices; j++)
	{
	  double *coord = GetVertex(j);
	  for(int i=0; i<GetDim(); i++)
	    {
	      out << coord[i] << " ";
	    }
	  out << endl;
	}
    }

  if (type==Element::QUADRILATERAL)
    {

      out << "areamesh2\n\n" << NumOfBdrElements << endl;
      for(int j=0; j<NumOfBdrElements; j++)
	{
	  int nv = GetBdrElementNVertices(j);
	  out << GetBdrAttribute(j) << "     ";
	  Array<int> ind;
	  GetBdrElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    out << ind[i]+1 << "   ";
	  out << endl;
	}

      out << NumOfElements << endl;
      for(int j=0; j<NumOfElements; j++)
	{
	  int nv = GetElementNVertices(j);
	  out << GetAttribute(j) << "  ";
	  out << nv  << "   ";
	  Array<int> ind;
	  GetElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    out << ind[i]+1 << "  ";
	  out << endl;
	}
      
      // Create vertices.
      out << NumOfVertices << endl;
      for(int j=0; j<NumOfVertices; j++)
	{
	  double *coord = GetVertex(j);
	  for(int i=0; i<GetDim(); i++)
	    {
	      out << coord[i] << " ";
	    }
	  out << endl;
	}
    }
  
  else if (type==Element::TETRAHEDRAL)
    {
      out << "NETGEN_Neutral_Format\n";
      out << NumOfVertices << '\n';
      // Create vertices.
      for(int j=0; j<NumOfVertices; j++)
	{
	  for(int i=0; i<GetDim(); i++)
	    {
	      double coord = GetVertex(j, i);
	      out << coord << " ";
	    }
	  out << '\n';
	}
      out << NumOfElements << '\n';
      for(int j=0; j<NumOfElements; j++)
	{
	  int nv = GetElementNVertices(j);
	  //	  out << nv  << " " <<;
	  out << GetAttribute(j)+1 << "  ";
	  Array<int> ind;
	  GetElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    {
	      out << ind[i]+1 << " ";
	    }
	  out << '\n';
	}

  
      out << NumOfBdrElements << '\n';
      for(int j=0; j<NumOfBdrElements; j++)
	{
	  int nv = GetBdrElementNVertices(j);
	  //out << nv << " ";
	  out << GetBdrAttribute(j)+1 << "  ";
	  Array<int> ind;
	  GetBdrElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    {
	      out << ind[i]+1 << " ";
	    }
	  out << '\n';
	}
      
    }
}

//============================================================================

void Mesh::ParaviewPrint(int type) const
{
  ofstream out;
  out.open("mesh.vtk");
  out << "# vtk DataFile Version 2.0" << endl;
  out << "Unstructured Grid" << endl;
  out << "ASCII" << endl;
  out << "DATASET UNSTRUCTURED_GRID" << endl;
  out << "POINTS " << NumOfVertices << " double" << endl;


  if (type==Element::QUADRILATERAL)
    {
      for(int j=0; j<NumOfVertices; j++)
	{
	  for(int i=0; i<GetDim(); i++)
	    {
	      double coord = GetVertex(j, i);
	      out << coord << " ";
	    }
	  out << 0.0 << endl;
	}
      
      out << endl;
      out << "CELLS " << NumOfElements << " " << 5*NumOfElements << endl;

      for(int j=0; j<NumOfElements; j++)
	{
	  int nv = GetElementNVertices(j);
	  out << 4 << "     ";
	  Array<int> ind;
	  GetElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    {
	      out << ind[i] << "    ";
	    }
	  out << endl;
	}
      out << endl;
      out << "CELL_TYPES " << NumOfElements << endl;
      for(int j=0; j<NumOfElements; j++)
	{
	  out << 9 << endl;
	}
    }


  if (type==Element::TRIANGLE)
    {
      for(int j=0; j<NumOfVertices; j++)
	{
	  for(int i=0; i<GetDim(); i++)
	    {
	      double coord = GetVertex(j, i);
	      out << coord << " ";
	    }
	  out << 0.0 << endl;
	}
      
      out << endl;
      out << "CELLS " << NumOfElements << " " << 4*NumOfElements << endl;

      for(int j=0; j<NumOfElements; j++)
	{
	  int nv = GetElementNVertices(j);
	  out << 3 << "     ";
	  Array<int> ind;
	  GetElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    {
	      out << ind[i] << "    ";
	    }
	  out << endl;
	}
      out << endl;
      out << "CELL_TYPES " << NumOfElements << endl;
      for(int j=0; j<NumOfElements; j++)
	{
	  out << 5 << endl;
	}
    }
  
  else if (type==Element::TETRAHEDRAL)
    {
      for(int j=0; j<NumOfVertices; j++)
	{
	  for(int i=0; i<GetDim(); i++)
	    {
	      double coord = GetVertex(j, i);
	      out << coord << " ";
	    }
	  out << endl;
	}

      out << endl;
      out << "CELLS " << NumOfElements << " " << 5*NumOfElements << endl;

      for(int j=0; j<NumOfElements; j++)
	{
	  int nv = GetElementNVertices(j);
	  out << 4 << "     ";
	  Array<int> ind;
	  GetElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    {
	      out << ind[i] << "    ";
	    }
	  out << endl;
	}
      out << endl;
      out << "CELL_TYPES " << NumOfElements << endl;
      for(int j=0; j<NumOfElements; j++)
	{
	  out << 10 << endl;
	}
    }

  out.close();

}

//============================================================================

void Mesh::ParaviewPrint(const Vector &solution, int type) const
{
  ofstream out;
  out.open("sol.vtk");
  out << "# vtk DataFile Version 2.0" << endl;
  out << "Unstructured Grid" << endl;
  out << "ASCII" << endl;
  out << "DATASET UNSTRUCTURED_GRID" << endl;
  out << "POINTS " << NumOfVertices << " double" << endl;


  if (type==Element::QUADRILATERAL)
    {
      for(int j=0; j<NumOfVertices; j++)
	{
	  for(int i=0; i<GetDim(); i++)
	    {
	      double coord = GetVertex(j, i);
	      out << coord << " ";
	    }
	  out << 0.0 << endl;
	}
      
      out << endl;
      out << "CELLS " << NumOfElements << " " << 5*NumOfElements << endl;

      for(int j=0; j<NumOfElements; j++)
	{
	  int nv = GetElementNVertices(j);
	  out << 4 << "     ";
	  Array<int> ind;
	  GetElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    {
	      out << ind[i] << "    ";
	    }
	  out << endl;
	}
      out << endl;
      out << "CELL_TYPES " << NumOfElements << endl;
      for(int j=0; j<NumOfElements; j++)
	{
	  out << 9 << endl;
	}

      out << endl;
      out << "POINT_DATA " << NumOfVertices << endl;
      out << "SCALARS pval DOUBLE" << endl;
      out << "LOOKUP_TABLE default" << endl;

      for(int j=0; j<NumOfVertices; j++)
	{
	  out << solution(j) << endl;
	}
    }


  if (type==Element::TRIANGLE)
    {
      for(int j=0; j<NumOfVertices; j++)
	{
	  for(int i=0; i<GetDim(); i++)
	    {
	      double coord = GetVertex(j, i);
	      out << coord << " ";
	    }
	  out << 0 << endl;
	}
      
      out << endl;
      out << "CELLS " << NumOfElements << " " << 4*NumOfElements << endl;

      for(int j=0; j<NumOfElements; j++)
	{
	  int nv = GetElementNVertices(j);
	  out << 3 << "     ";
	  Array<int> ind;
	  GetElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    {
	      out << ind[i] << "    ";
	    }
	  out << endl;
	}
      out << endl;
      out << "CELL_TYPES " << NumOfElements << endl;
      for(int j=0; j<NumOfElements; j++)
	{
	  out << 5 << endl;
	}

      out << endl;
      out << "POINT_DATA " << NumOfVertices << endl;
      out << "SCALARS pval DOUBLE" << endl;
      out << "LOOKUP_TABLE default" << endl;

      for(int j=0; j<NumOfVertices; j++)
	{
	  out << solution(j) << endl;
	}
    }
  
  else if (type==Element::TETRAHEDRAL)
    {
      for(int j=0; j<NumOfVertices; j++)
	{
	  for(int i=0; i<GetDim(); i++)
	    {
	      double coord = GetVertex(j, i);
	      out << coord << " ";
	    }
	  out << endl;
	}

      out << endl;
      out << "CELLS " << NumOfElements << " " << 5*NumOfElements << endl;

      for(int j=0; j<NumOfElements; j++)
	{
	  int nv = GetElementNVertices(j);
	  out << 4 << "     ";
	  Array<int> ind;
	  GetElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    {
	      out << ind[i] << "    ";
	    }
	  out << endl;
	}
      out << endl;
      out << "CELL_TYPES " << NumOfElements << endl;
      for(int j=0; j<NumOfElements; j++)
	{
	  out << 10 << endl;
	}

      out << endl;
      out << "POINT_DATA " << NumOfVertices << endl;
      out << "SCALARS pval DOUBLE" << endl;
      out << "LOOKUP_TABLE default" << endl;

      for(int j=0; j<NumOfVertices; j++)
	{
	  out << solution(j) << endl;
	}

    }

  out.close();

}

//============================================================================
void Mesh::ParaviewPrint(ofstream &out, const Vector &solution, double scalex, double scaley) const
{
  out << "# vtk DataFile Version 2.0" << endl;
  out << "Unstructured Grid" << endl;
  out << "ASCII" << endl;
  out << "DATASET UNSTRUCTURED_GRID" << endl;
  out << "POINTS " << NumOfVertices << " double" << endl;


  if (GetMType()==Element::QUADRILATERAL)
    {
      for(int j=0; j<NumOfVertices; j++)
	{
	  out << scalex*GetVertex(j, 0) << " " << scaley*GetVertex(j, 1) << " " << 0.0 << endl;
	}
      
      out << endl;
      out << "CELLS " << NumOfElements << " " << 5*NumOfElements << endl;

      for(int j=0; j<NumOfElements; j++)
	{
	  int nv = GetElementNVertices(j);
	  out << 4 << "     ";
	  Array<int> ind;
	  GetElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    {
	      out << ind[i] << "    ";
	    }
	  out << endl;
	}
      out << endl;
      out << "CELL_TYPES " << NumOfElements << endl;
      for(int j=0; j<NumOfElements; j++)
	{
	  out << 9 << endl;
	}

      out << endl;
      out << "POINT_DATA " << NumOfVertices << endl;
      out << "SCALARS pval DOUBLE" << endl;
      out << "LOOKUP_TABLE default" << endl;

      for(int j=0; j<NumOfVertices; j++)
	{
	  out << solution(j) << endl;
	}
    }


  if (GetMType()==Element::TRIANGLE)
    {
      for(int j=0; j<NumOfVertices; j++)
	{
	  out << scalex*GetVertex(j, 0) << " " << scaley*GetVertex(j, 1) << " " << 0.0 << endl;
	}
      
      out << endl;
      out << "CELLS " << NumOfElements << " " << 4*NumOfElements << endl;

      for(int j=0; j<NumOfElements; j++)
	{
	  int nv = GetElementNVertices(j);
	  out << 3 << "     ";
	  Array<int> ind;
	  GetElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    {
	      out << ind[i] << "    ";
	    }
	  out << endl;
	}
      out << endl;
      out << "CELL_TYPES " << NumOfElements << endl;
      for(int j=0; j<NumOfElements; j++)
	{
	  out << 5 << endl;
	}

      out << endl;
      out << "POINT_DATA " << NumOfVertices << endl;
      out << "SCALARS pval DOUBLE" << endl;
      out << "LOOKUP_TABLE default" << endl;

      for(int j=0; j<NumOfVertices; j++)
	{
	  out << solution(j) << endl;
	}
    }
  
  else if (GetMType()==Element::TETRAHEDRAL)
    {
      for(int j=0; j<NumOfVertices; j++)
	{
	  for(int i=0; i<GetDim(); i++)
	    {
	      double coord = GetVertex(j, i);
	      out << coord << " ";
	    }
	  out << endl;
	}

      out << endl;
      out << "CELLS " << NumOfElements << " " << 5*NumOfElements << endl;

      for(int j=0; j<NumOfElements; j++)
	{
	  int nv = GetElementNVertices(j);
	  out << 4 << "     ";
	  Array<int> ind;
	  GetElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    {
	      out << ind[i] << "    ";
	    }
	  out << endl;
	}
      out << endl;
      out << "CELL_TYPES " << NumOfElements << endl;
      for(int j=0; j<NumOfElements; j++)
	{
	  out << 10 << endl;
	}

      out << endl;
      out << "POINT_DATA " << NumOfVertices << endl;
      out << "SCALARS pval DOUBLE" << endl;
      out << "LOOKUP_TABLE default" << endl;

      for(int j=0; j<NumOfVertices; j++)
	{
	  out << solution(j) << endl;
	}

    }

  out.close();

}

//============================================================================

void Mesh::ParaviewPrintDeformedMesh(ofstream &outfile, const Vector &displacement, double scalex, double scaley) const
{
  if (GetMType()==Element::TRIANGLE)
    {
      outfile << "# vtk DataFile Version 2.0" << endl;
      outfile << "Unstructured Grid" << endl;
      outfile << "ASCII" << endl;
      outfile << "DATASET UNSTRUCTURED_GRID" << endl;
      outfile << "POINTS " << NumOfVertices << " double" << endl;
    
      Array<double *> nodes(2);
      nodes[0] = new double[NumOfVertices];
      nodes[1] = new double[NumOfVertices];
  
      for(int i=0; i<NumOfVertices; i++) 
	{
	  nodes[0][i] = 0.0;  
	  nodes[1][i] = 0.0;
	}



      double *pt;
      for(int i=0; i<NumOfVertices; i++)
	{
	  pt = GetVertex(i); 
	  for(int k=0; k<2; k++)
	    {
	      nodes[k][i] += pt[k];
	    } 
	}

      for(int i=0; i<NumOfVertices; i++)
	{ 
	  nodes[0][i] += displacement(i);
	  nodes[1][i] += displacement(i+NumOfVertices);
	}

      for(int j=0; j<NumOfVertices; j++)
	{
	  outfile << scalex*nodes[0][j] << " " << scaley*nodes[1][j] << " 0.0 " << endl;
	}
  
      outfile << endl;
      outfile << "CELLS " << NumOfElements << " " << 4*NumOfElements << endl;
  
      for(int j=0; j<NumOfElements; j++)
	{
	  int nv = GetElementNVertices(j);
	  outfile << 3 << "     ";
	  Array<int> ind;
	  GetElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    {
	      outfile << ind[i] << "    ";
	    }
	  outfile << endl;
	}
      outfile << endl;
      outfile << "CELL_TYPES " << NumOfElements << endl;
      for(int j=0; j<NumOfElements; j++)
	{
	  outfile << 5 << endl;
	}
      outfile.close();

      for(int i=0; i<2; i++)
	delete []nodes[i];

    }


  if (GetMType()==Element::QUADRILATERAL)
    {
      outfile << "# vtk DataFile Version 2.0" << endl;
      outfile << "Unstructured Grid" << endl;
      outfile << "ASCII" << endl;
      outfile << "DATASET UNSTRUCTURED_GRID" << endl;
      outfile << "POINTS " << NumOfVertices << " double" << endl;

    
      Array<double *> nodes(2);
      nodes[0] = new double[NumOfVertices];
      nodes[1] = new double[NumOfVertices];
  
      for(int i=0; i<NumOfVertices; i++)
	{
	  nodes[0][i] = 0.0;
	  nodes[1][i] = 0.0;
	}

      double *pt;
      for(int i=0; i<NumOfVertices; i++)
	{
	  pt = GetVertex(i); 
	  for(int k=0; k<2; k++)
	    {
	      nodes[k][i] += pt[k];
	    } 
	}

     for(int i=0; i<NumOfVertices; i++)
	{ 
	  nodes[0][i] += displacement(i);
	  nodes[1][i] += displacement(i+NumOfVertices);
	}

      for(int j=0; j<NumOfVertices; j++)
	{
	  outfile << scalex*nodes[0][j] << " " << scaley*nodes[1][j] << " 0.0 " << endl;
	}

      outfile << endl;
      outfile << "CELLS " << NumOfElements << " " << 5*NumOfElements << endl;
  
      for(int j=0; j<NumOfElements; j++)
	{
	  int nv = GetElementNVertices(j);
	  outfile << 4 << "     ";
	  Array<int> ind;
	  GetElementVertices(j, ind);
	  for(int i=0; i<nv; i++)
	    {
	      outfile << ind[i] << "    ";
	    }
	  outfile << endl;
	}
      outfile << endl;
      outfile << "CELL_TYPES " << NumOfElements << endl;
      for(int j=0; j<NumOfElements; j++)
	{
	  outfile << 9 << endl;
	}

      for(int i=0; i<2; i++)
	delete []nodes[i];

    }
}

//============================================================================
void Mesh::Print(ofstream &out, const Vector &solution) const
{

  out << "solution u -size=" << NumOfVertices << " -type=nodal -components=1" 
      << endl;
  
  for (int i=0; i<NumOfVertices; i++)
    out << solution(i) << endl;
}

//============================================================================

template <int Dim>
Table *TMesh<Dim>::GetEdgeToElementTable() const
{
  if (edge_to_el != NULL)
    {
      return edge_to_el;
    }

  int i, j, nv;
  Array<int> edges, cor;

  edge_to_el = new Table;
  
  edge_to_el->MakeI(NumOfEdges);
  
  for (i=0; i<NumOfElements; i++)
    {
      nv = 3;
      //if (elements[i]->GetGeometryType()==Element::TRIANGLE) nv = 3;
      GetElementEdges(i, edges);
      for (j=0; j<nv; j++) { edge_to_el->AddAColumnInRow(edges[j]); }
    }
  
  edge_to_el->MakeJ();
  
  for (i = 0; i<NumOfElements; i++)
    {
      nv = 3;
      //if (elements[i]->GetGeometryType()==Element::TRIANGLE) nv = 3;
      GetElementEdges(i, edges);
      for (j=0; j<nv; j++)
	edge_to_el->AddConnection (edges[j], i);
    }
  
  edge_to_el->ShiftUpI();
  edge_to_el->Finalize();
  
  return edge_to_el;
}

//============================================================================

template <int Dim> int TMesh<Dim>::GetElementToEdgeTable(Table &e_to_f, Array<int> &be_to_edge) const
{
  Array<int> v;
  int NumberOfEdges = -1;
  
  if (Dim==2) {
    STable v_to_v(NumOfVertices, 15);
    GetVertexToVertexTable(v_to_v);
    
    NumberOfEdges = v_to_v.Size_of_connections();
    
    // Fill the element to edge table
    e_to_f.SetSize(NumOfElements, 4);
    for (int i=0; i<NumOfElements; i++)
      {
	element[i]->GetVertices(v);
	for (int j=1; j<v.Size(); j++)
	  e_to_f.Push( i, v_to_v(v[j-1], v[j]) );
	e_to_f.Push(i, v_to_v(v[0], v[v.Size()-1]));
      }
    e_to_f.Finalize();


    // Initialize the indices for the boundary elements.
    be_to_edge.SetSize(NumOfBdrElements);
    for (int i=0; i<NumOfBdrElements; i++)
      {
	bdrelement[i]->GetVertices(v);
	be_to_edge[i] = v_to_v(v[0], v[1]);
      }

  }
  // Return the number of edges
  return NumberOfEdges;
}

//============================================================================

template <int Dim> void TMesh<Dim>::GetVertexToVertexTable( STable & v_to_v ) const
{
  Array<int> v;

  if (Dim==2) {
    for( int i=0; i<NumOfElements; i++ )
      {
	element[i]->GetVertices(v);
	for(int j=1; j<v.Size(); j++ )
	  v_to_v.Push(v[j-1],v[j]);
	v_to_v.Push(v[0], v[v.Size()-1]);
      }
    v_to_v.Finalize();
  }
}

//============================================================================

template <int Dim> void TMesh<Dim>::GetElementEdges(int i, Array<int> &edges) const
{
  if (Dim==2) {
    if (el_to_edge)
      el_to_edge->GetRow(i, edges);
    else
      {
	static int time = 0;
	if (time == 0)
	  cerr << "\nDMesh<Dim>::GetElementEdges(...) element to edge table "
	    "is not generated." << endl;
	time++;
      }
  }
}

//============================================================================

template <int Dim> void TMesh<Dim>::GetEdgeElements(int i, Array<int> &elements) const
{
  if (Dim==2) {
    if (edge_to_el)
      edge_to_el->GetRow(i, elements);
    else
      {
	static int time = 0;
	if (time == 0)
	  cerr << "\nTMesh<Dim>::GetEdgeElements(...) edge to element table "
	    "is not generated." << endl;
	time++;
      }
  }
}

//============================================================================

template <int Dim> void TMesh<Dim>::GetElementVerticesCoord(int i, Array<double *> &coord) const
{
  int nv = element[i]->GetNVertices();
  Array<int> v;
  element[i]->GetVertices(v);
  
  for (int k=0; k<Dim; k++)
    {
      for (int l=0; l<nv; l++)
	{
	  coord[k][l] = vertices[v[l]](k);
	}
    }
}

//============================================================================

template <int Dim> void TMesh<Dim>::GetBdrElementVerticesCoord(int i, Array<double *> &coord) const
{
  int nv = bdrelement[i]->GetNVertices();
  Array<int> v;
  bdrelement[i]->GetVertices(v);
  for (int k=0; k<Dim; k++)
    for (int l=0; l<nv; l++)
      coord[k][l] = vertices[v[l]](k);
}

//============================================================================
template <int Dim> void TMesh<Dim>::UpdateMesh(const Vector &displacement)
{

  for(int i=0; i<NumOfVertices; i++)
    { 
      vertices[i].UpdateVertex(0, vertices[i](0)+displacement(i));
      vertices[i].UpdateVertex(1, vertices[i](1)+displacement(i+NumOfVertices));
    }
}

//============================================================================

template <int Dim> TMesh<Dim>::TMesh(int n, double xl, double xr)
{
  double h = (xr-xl)/ ((double) n);
  double  NumOfVertices = n + 1;
  NumOfElements = n;
  NumOfBdrElements = 2;
  vertices.SetSize(NumOfVertices);
  element.SetSize(NumOfElements);
  bdrelement.SetSize(NumOfBdrElements);

  // Create vertices.
  double cx = xl;
  for (int i=0; i<NumOfVertices; i++)
    {
      vertices[i](0) = cx;
      cx += h;
    }

  // Create element.
  int *v = new int[2];
  for (int i=0; i<NumOfElements; i++)
    {
      v[0] = i;
      v[1] = i+1;
      int attr = -1;
      element[i] = new Segment(v, attr);
    }
  delete []v;

  v = new int[1];
  v[0] = 0;
  bdrelement[0] = new Point(v, 0);
  v[0] = n;
  bdrelement[1] = new Point(v, 1);
}

//============================================================================

template <int Dim> void TMesh<Dim>::ConstructRectangularElements(Array<int> &n, Array<double> &L,
								 bool generateedges)
{     
  mtype = Element::QUADRILATERAL;

  int nx = n[0];
  int ny = n[1];
  double hx = (L[1] - L[0]) / ( (double) nx );
  double hy = (L[3] - L[2]) / ( (double) ny );

  Nx = nx;
  NumOfBdrs = 4;
  NumOfVertices = (nx + 1)*(ny + 1);
  NumOfElements = nx*ny;
  NumOfBdrElements = 2 * (nx + ny);

  vertices.SetSize(NumOfVertices);
  element.SetSize(NumOfElements);
  bdrelement.SetSize(NumOfBdrElements);    
  NumOfBdryElem.SetSize(NumOfBdrs);
  BdrytoGlobalIndex.SetSize(NumOfBdrs);

  //vertice coordinates
  int ind=0; //useful index
  for(int j=0; j<=ny; j++)
    {
      for (int i=0; i<=nx; i++)
	{
	  vertices[ind](0) = L[0] + i*hx;
	  vertices[ind](1) = L[2] + j*hy;
	  ind++;
	}
    }       
   
  
  //elements
  ind = 0;
  int v[4];//the 4 indeces of the elements
  int attr = -1;

  for(int j=0; j<ny; j++)
    {
      for (int i=0; i<nx; i++)
	{
	  v[0] = i + j*(nx+1);
	  v[1] = v[0] + 1;
	  v[2] = v[1] + (nx + 1);
	  v[3] = v[2] - 1;

	  element[ind] = new Quadrilateral(v, attr);
	  ind++;
	}
    }   
  
  //boundary elements
  int vv[2];
  for (int i=0; i<nx; i++)
    {
      //bottom (attr = 0)
      vv[0] = i;
      vv[1] = i + 1;
      attr = 0;

      bdrelement[i] = new Segment(vv, attr);

      //top (attr = 2)
      vv[0] += ny*(nx+1);
      vv[1] += ny*(nx+1);
      attr = 2;

      int bdrelind = i + nx + ny;
      bdrelement[bdrelind] = new Segment(vv, attr);
    }

  for (int j=0; j<ny; j++)
    {
      //left (attr = 3)
      vv[0] = j*(nx+1);
      vv[1] = (j+1)*(nx+1);
      attr = 3;

      int bdrelind = j + 2*nx + ny;
      bdrelement[bdrelind] = new Segment(vv, attr);

      //right (attr = 1)
      vv[0] += nx;
      vv[1] += nx;
      attr = 1;

      bdrelind = j + nx;
      bdrelement[bdrelind] = new Segment(vv, attr);
    }

  BdrytoGlobalIndex[0] = new int[nx+1];
  BdrytoGlobalIndex[2] = new int[nx+1];
  for (int i=0; i<=nx; i++)
    {
      BdrytoGlobalIndex[0][i] = i;
      BdrytoGlobalIndex[2][i] = i + ny*(nx+1);
    }

  BdrytoGlobalIndex[1] = new int[ny+1];
  BdrytoGlobalIndex[3] = new int[ny+1];
  for (int j=0; j<=ny; j++)
    {
      BdrytoGlobalIndex[1][j] = nx + j*(nx+1);
      BdrytoGlobalIndex[3][j] = j*(nx+1);
    }

  NumOfBdryElem[0] = NumOfBdryElem[2] = nx;
  NumOfBdryElem[1] = NumOfBdryElem[3] = ny;


  // set vertex to element map (The following probably needs to be revised, I don't know what its for)
  vert_to_ele.SetSize(NumOfVertices);
  
  int *count = new int[NumOfVertices];
  for(int i=0; i<NumOfVertices; i++)
    count[i] = 0;

  for(int i=0; i<NumOfElements; i++)
    {
      Array<int> ind;
      GetElementVertices(i, ind);
      for(int j=0; j<4; j++)
	count[ind[j]]++;
    }

  for(int i=0; i<NumOfVertices; i++)
    vert_to_ele[i] = new int[count[i]+1];

  for(int i=0; i<NumOfVertices; i++)
    vert_to_ele[i][0] = count[i];

  int *index = new int[NumOfVertices];
  for(int i=0; i<NumOfVertices; i++)
    index[i] = 1;

  for(int i=0; i<NumOfElements; i++)
    {
      Array<int> ind;
      GetElementVertices(i, ind);
      for(int j=0; j<4; j++)
	{
	vert_to_ele[ind[j]][index[ind[j]]] = i;
	index[ind[j]]++;
	}
    }

  if (generateedges)
    {
      el_to_edge = new Table(GetNE(), 4);
      NumOfEdges = GetElementToEdgeTable(*el_to_edge, be_to_edge);
    }
  else
    {
      el_to_edge = NULL;
      NumOfEdges = 0;
    }

  delete []index;
  delete []count;

}

//=============================================================================

template <int Dim> void TMesh<Dim>::ConstructTriangularElements(Array<int> &n, Array<double> &L,
								 bool generateedges)
{

  mtype = Element::TRIANGLE;

  // Array n = (n[0], n[1]) = (n_x, n_y)
  // Array L = (L[0], L[1], L[2], L[3]) =  (x1, x2, y1, y2)
  /*

  double dx = (L[1] - L[0]) / ( (double) n[0] );
  double dy = (L[3] - L[2]) / ( (double) n[1] );
  NumOfBdrs = 4;
  NumOfVertices = (n[0] + 1) * (n[1] + 1);
  NumOfElements = 2 * n[0] * n[1];
  NumOfBdrElements = 2 * (n[0] + n[1]);
  NumOfEdges = n[0]*(n[1]+1) + (2*n[0] + 1)*n[1];
  vertices.SetSize(NumOfVertices);
  element.SetSize(NumOfElements);
  bdrelement.SetSize(NumOfBdrElements);
  
  // Create vertices
  int jthrow = 0;
  int ithcol = 0;
  for (int i=0; i<NumOfVertices; i++)
    {
      jthrow = i / (n[0]+1);
      ithcol = i % (n[0]+1);      
      vertices[i](0) = L[0] + ( (double) ithcol ) * dx;
      vertices[i](1) = L[2] + ( (double) jthrow ) * dy;
    }
  */

  
  int nx = n[0];
  int ny = n[1];
  double hx = (L[1] - L[0]) / ( (double) nx );
  double hy = (L[3] - L[2]) / ( (double) ny );

  Nx = nx;

  NumOfBdrs = 4;
  NumOfVertices = (nx + 1)*(ny + 1);
  NumOfElements = 2*nx*ny;
  NumOfBdrElements = 2 * (nx + ny);
  NumOfEdges = nx*(ny+1) + (2*nx + 1)*ny;

  vertices.SetSize(NumOfVertices);
  element.SetSize(NumOfElements);
  bdrelement.SetSize(NumOfBdrElements);    
  NumOfBdryElem.SetSize(NumOfBdrs);
  BdrytoGlobalIndex.SetSize(NumOfBdrs);
 
  //vertice coordinates
  int ind=0; //useful index
  for(int j=0; j<=ny; j++)
    {
      for (int i=0; i<=nx; i++)
	{
	  vertices[ind](0) = L[0] + i*hx;
	  vertices[ind](1) = L[2] + j*hy;
	  ind++;
	}
    }    
   
  // Create elements
  int *v = new int[3];
  int attr = -1;
  for (int i=0; i<NumOfElements; i++)
    {
      int jthrow = i / ( 2*nx );
      if(i%2 == 0) // Lower triangle
	{
	  v[0] = i/2 + jthrow;
	  v[1] = v[0] + 1;
	  v[2] = v[0] + (nx + 1);
	}
      else // Upper triangle
	{
	  v[0] = i/2 + 1 + jthrow;
	  v[1] = v[0] + (nx + 1);
	  v[2] = v[1] - 1;
	}
     
      element[i] = new Triangle(v, attr);
    }
  delete []v;
  /*
  // Create boundary elements
  v = new int[2];
  for (int i=0; i<NumOfBdrElements; i++) // Counter-clockwise
    {
      if(i < n[0]) // Bottom boundary
	{
	  v[0] = i;
	  v[1] = i + 1;          
	  bdrelement[i] = new Segment(v, 0);
	}
      else if(i < n[0]+n[1]) // Right boundary
        {
	  v[0] = n[0] + (i - n[0]) * (n[0] + 1);
	  v[1] = v[0] + (n[0] + 1);         
	  bdrelement[i] = new Segment(v, 1);
        }
      else if(i < 2*n[0]+n[1]) // Top boundary
        {
	  v[0] = NumOfVertices - 1 - ( i - n[0] - n[1] );
	  v[1] = v[0] - 1;
	  bdrelement[i] = new Segment(v, 2);
        }
	else // Left boundary; v[1] is the top index
        {
	  v[1] = n[1] * (n[0]+1) - ( i - 2*n[0] - n[1] ) * (n[0]+1);
	  v[0] = v[1] - (n[0]+1);
	  bdrelement[i] = new Segment(v, 3);
        }
    }
  delete []v;
*/
  //boundary elements
  int vv[2];
  for (int i=0; i<nx; i++)
    {
      //bottom (attr = 0)
      vv[0] = i;
      vv[1] = i + 1;
      attr = 0;

      bdrelement[i] = new Segment(vv, attr);

      //top (attr = 2)
      vv[0] += ny*(nx+1);
      vv[1] += ny*(nx+1);
      attr = 2;

      int bdrelind = i + nx + ny;
      bdrelement[bdrelind] = new Segment(vv, attr);
    }

  for (int j=0; j<ny; j++)
    {
      //left (attr = 3)
      vv[0] = j*(nx+1);
      vv[1] = (j+1)*(nx+1);
      attr = 3;

      int bdrelind = j + 2*nx + ny;
      bdrelement[bdrelind] = new Segment(vv, attr);

      //right (attr = 1)
      vv[0] += nx;
      vv[1] += nx;
      attr = 1;

      bdrelind = j + nx;
      bdrelement[bdrelind] = new Segment(vv, attr);
    }

  BdrytoGlobalIndex[0] = new int[nx+1];
  BdrytoGlobalIndex[2] = new int[nx+1];
  for (int i=0; i<=nx; i++)
    {
      BdrytoGlobalIndex[0][i] = i;
      BdrytoGlobalIndex[2][i] = i + ny*(nx+1);
    }

  BdrytoGlobalIndex[1] = new int[ny+1];
  BdrytoGlobalIndex[3] = new int[ny+1];
  for (int j=0; j<=ny; j++)
    {
      BdrytoGlobalIndex[1][j] = nx + j*(nx+1);
      BdrytoGlobalIndex[3][j] = j*(nx+1);
    }

  NumOfBdryElem[0] = NumOfBdryElem[2] = nx;
  NumOfBdryElem[1] = NumOfBdryElem[3] = ny;


  // set vertex to element map
  vert_to_ele.SetSize(NumOfVertices);
  
  int *count = new int[NumOfVertices];
  for(int i=0; i<NumOfVertices; i++)
    count[i] = 0;

  for(int i=0; i<NumOfElements; i++)
    {
      Array<int> ind;
      GetElementVertices(i, ind);
      for(int j=0; j<3; j++)
	count[ind[j]]++;
    }

  for(int i=0; i<NumOfVertices; i++)
    vert_to_ele[i] = new int[count[i]+1];

  for(int i=0; i<NumOfVertices; i++)
    vert_to_ele[i][0] = count[i];

  int *index = new int[NumOfVertices];
  for(int i=0; i<NumOfVertices; i++)
    index[i] = 1;

  for(int i=0; i<NumOfElements; i++)
    {
      Array<int> ind;
      GetElementVertices(i, ind);
      for(int j=0; j<3; j++)
	{
	vert_to_ele[ind[j]][index[ind[j]]] = i;
	index[ind[j]]++;
	}
    }

  if (generateedges)
    {
      el_to_edge = new Table(GetNE(), 4);
      NumOfEdges = GetElementToEdgeTable(*el_to_edge, be_to_edge);
    }
  else
    {
      el_to_edge = NULL;
      NumOfEdges = 0;
    }
	 
  delete []index;
  delete []count;
}

//=============================================================================

template <int Dim> void TMesh<Dim>::ConstructTetrahedralElements(Array<int> &n, Array<double> &L,
								  bool generateedges)
{
  mtype = Element::TETRAHEDRAL;
  double nx = n[0];
  double ny = n[1];
  double nz = n[2];

  double hx = (L[1] - L[0]) / ( (double) nx );
  double hy = (L[3] - L[2]) / ( (double) ny );
  double hz = (L[5] - L[4]) / ( (double) nz );

  NumOfVertices = (nx + 1) * (ny + 1) * (nz + 1);
  NumOfElements = 5 * nx * ny * nz;
  NumOfBdrElements = 4 * (nx*ny + ny*nz + nz*nx);
  NumOfBdrs = 6;

  vertices.SetSize(NumOfVertices);
  element.SetSize(NumOfElements);
  bdrelement.SetSize(NumOfBdrElements);

  int index = 0;
  for (int k=0; k<=nz; k++)
    {
      for (int j=0; j<=ny; j++)
	{
	  for (int i=0; i<=nx; i++)
	    {   
	      double x = i*hx;
	      double y = j*hy;
	      double z = k*hz;
	      vertices[index](0) = x;
	      vertices[index](1) = y;
	      vertices[index](2) = z;
	      index++;
	    }
	}
    }

  //elements
  int v[4];//the 4 indeces of the elements
  int attr = -1;

  index = 0;
  for (int k=0; k<nz; k++)
    {
      for (int j=0; j<ny; j++)
	{
	  for (int i=0; i<nx; i++)
	    {   
	      v[0] = k*(nx+1)*(ny+1) + j*(nx+1) + i;
	      v[1] = v[0]+1;
	      v[2] = v[0] + nx+1;
	      v[3] = v[0] + (nx+1)*(ny+1);

	      element[index] = new Tetrahedral(v, attr); 
	      index++;

	      v[0] = k*(nx+1)*(ny+1) + (j+1)*(nx+1) + i+1;
	      v[1] = v[0]-1;
	      v[2] = v[0]-nx-1;
	      v[3] = v[0] + (nx+1)*(ny+1);

	      element[index] = new Tetrahedral(v, attr); 
	      index++;

	      v[3] =  k*(nx+1)*(ny+1) + j*(nx+1) + i + 1;
	      v[0] = v[3] + (nx+1)*(ny+1);
	      v[1] = v[0] - 1;
	      v[2] = v[0] + nx+1;


	      element[index] = new Tetrahedral(v, attr); 
	      index++;

	      v[3] =  k*(nx+1)*(ny+1) + (j+1)*(nx+1) + i;
	      v[0] = v[3] + (nx+1)*(ny+1);
	      v[1] = v[0] + 1;
	      v[2] = v[0] - nx-1;


	      element[index] = new Tetrahedral(v, attr); 
	      index++;

	      v[2] =  k*(nx+1)*(ny+1) + j*(nx+1) + i+1;
	      v[1] = v[2] + nx;
	      v[3] = (k+1)*(nx+1)*(ny+1) + j*(nx+1) + i;
	      v[0] = v[3] + nx + 2;


	      element[index] = new Tetrahedral(v, attr); 
	      index++;
	    }
	}
    }


  //boundary elements
  int vv[3];
  int attr1;
  index=0;

  //bottom
  attr1=0;
  for(int j=0; j<ny; j++)
    {
      for(int i=0; i<nx; i++)
	{
	  int ind = j*(nx+1) + i;
	  vv[0] = ind;
	  vv[1] = ind+1;
	  vv[2] = ind + nx+1;

	  bdrelement[index] = new Triangle(vv, attr1);
	  index++;

	  vv[0] = ind+nx+2;
	  vv[1] = ind+nx+1;
	  vv[2] = ind+1;

	  bdrelement[index] = new Triangle(vv, attr1);
	  index++;
	}
    }

  //top
  attr1 = 1;
  for(int j=0; j<ny; j++)
    {
      for(int i=0; i<nx; i++)
	{
	  int ind = j*(nx+1) + i + nz*(nx+1)*(ny+1);
	  vv[0] = ind+1;
	  vv[1] = ind+nx+2;
	  vv[2] = ind;

	  bdrelement[index] = new Triangle(vv, attr1);
	  index++;

	  vv[0] = ind+nx+1;
	  vv[1] = ind+nx+2;
	  vv[2] = ind;

	  bdrelement[index] = new Triangle(vv, attr1);
	  index++;
	}
    }


  attr1 = 2;
  for(int k=0; k<nz; k++)
    {
      for(int i=0; i<nx; i++)
	{
	  int ind = k*(nx+1)*(ny+1) + i;
	  vv[0] = ind;
	  vv[1] = ind+1;
	  vv[2] = ind + (nx+1)*(ny+1);

	  bdrelement[index] = new Triangle(vv, attr1);
	  index++;

	  vv[0] = ind + (nx+1)*(ny+1) + 1;
	  vv[1] = ind + (nx+1)*(ny+1);
	  vv[2] = ind+1;

	  bdrelement[index] = new Triangle(vv, attr1);
	  index++;
	}
    }

  attr1 = 3;
  for(int k=0; k<nz; k++)
    {
      for(int i=0; i<nx; i++)
	{
	  int ind = k*(nx+1)*(ny+1) + ny*(nx+1) + i;
	  vv[0] = ind + 1;
	  vv[1] = ind + (nx+1)*(ny+1) + 1;
	  vv[2] = ind;

	  bdrelement[index] = new Triangle(vv, attr1);
	  index++;

	  vv[0] = ind + (nx+1)*(ny+1);
	  vv[1] = ind + (nx+1)*(ny+1)+1; 
	  vv[2] = ind; 

	  bdrelement[index] = new Triangle(vv, attr1);
	  index++;
	}
    }


  attr1 = 4;
  for(int k=0; k<nz; k++)
    {
      for(int j=0; j<ny; j++)
	{
	  int ind = k*(nx+1)*(ny+1) + j*(nx+1);
	  vv[0] = ind;
	  vv[1] = ind+(nx+1);
	  vv[2] = ind + (nx+1)*(ny+1);

	  bdrelement[index] = new Triangle(vv, attr1);
	  index++;

	  vv[0] = ind + (nx+1)*(ny+1) + (nx+1);
	  vv[1] = ind + (nx+1)*(ny+1);
	  vv[2] = ind+(nx+1);

	  bdrelement[index] = new Triangle(vv, attr1);
	  index++;
	}
    }

  attr1 = 5;
  for(int k=0; k<nz; k++)
    {
      for(int j=0; j<ny; j++)
	{
	  int ind = k*(nx+1)*(ny+1) + j*(nx+1) + nx;
	  vv[0] = ind + (nx+1);
	  vv[1] = ind + (nx+1)*(ny+1) + (nx+1);
	  vv[2] = ind;

	  bdrelement[index] = new Triangle(vv, attr1);
	  index++;

	  vv[0] = ind + (nx+1)*(ny+1);
	  vv[1] = ind + (nx+1)*(ny+1) + (nx+1); 
	  vv[2] = ind;

	  bdrelement[index] = new Triangle(vv, attr1);
	  index++;
	}
    }

  }

//============================================================================

template <int Dim> void TMesh<Dim>::ConstructRectangularElements(Array<ifstream *> &files,
								 bool generateedges)
{
  //we only consider the four-point polygon to set up the attr of the boundary element
  //and to set up the attr of the boundary element, here we consider the simple case [0,1] by [0,1], four boundaries
    
  mtype = Element::QUADRILATERAL;
  int temp;
    
  //set up the vertices
  *files[0] >> NumOfVertices;
  vertices.SetSize(NumOfVertices);
    
  Array<int> vertex_attr(NumOfVertices);//store the attr of each vertex
    
  for (int i=0; i<NumOfVertices; i++)
    {
      *files[0] >> vertices[i](0) >> vertices[i](1) >> vertex_attr[i];
    }
    
  //one unkown number
  *files[0] >> temp;
    
  //set up the element
  *files[0] >> temp >> NumOfElements;
  element.SetSize(NumOfElements);
    
  int v[4];
  int attr;//here the attr means if the element contains the boundary element
  for(int i=0; i<NumOfElements; i++)
    {
      *files[0] >> v[0] >> v[1] >> v[2] >> v[3];
      v[0] -= 1;
      v[1] -= 1;
      v[2] -= 1;
      v[3] -= 1;
        
      int match = 0;
      for(int j=0; j<4; j++)
        {
	  if( vertex_attr[ v[j] ] == 2 || vertex_attr[ v[j] ] == 3 )
            {
	      match += 1;
	      break;
            }
        }
        
      attr = match;
      element[i] = new Quadrilateral(v, attr);
    }
    
  //maybe we can do this: ifstream *if_temp = *files[0]
  //set up the boundary element
  NumOfBdrElements = 0;
  int edge_data[NumOfElements][5];
  for(int i=0; i<NumOfElements; i++)
    {
      *files[0] >> edge_data[i][0] >> edge_data[i][1] >> edge_data[i][2] >> edge_data[i][3] >> edge_data[i][4];
      for(int j=0; j<5; j++)
        {
	  if(edge_data[i][j] == 2)
            {
	      NumOfBdrElements += 1;
            }
        }
    }

  bdrelement.SetSize(NumOfBdrElements);
  int vv[2];
  int bdr_ele_index = 0;
  for(int i=0; i<NumOfElements; i++)
    {
      if( GetAttribute(i) == 1 )
	{
	  for(int j=0; j<5; j++)
	    {
	      if(edge_data[i][j] == 2)
		{
		  Array<int> ind(4);
		  GetElementVertices(i, ind);
                
		  vv[0] = ind[j-1];
		  vv[1] = ind[j];
                    
		  int attr1 = -1;
		  if( GetVertex(vv[0], 0) == 0.0 && GetVertex(vv[1], 0) == 0.0  )
		    {
		      attr1 = 3;
		    }
		  else if( GetVertex(vv[0], 0) == 1.0 && GetVertex(vv[1], 0) == 1.0 )
		    {
		      attr1 = 1;
		    }
		  else if( GetVertex(vv[0], 1) == 0.0 && GetVertex(vv[1], 1) == 0.0 )
		    {
		      attr1 = 0;
		    }
		  else if( GetVertex(vv[0], 1) == 1.0 && GetVertex(vv[1], 1) == 1.0) 
		    {
		      attr1 = 2;
		    }
                
		  bdrelement[bdr_ele_index] = new Segment(vv, attr1);
		  bdr_ele_index += 1;    
		}
	    }
	}
    }
    
    
    
}

//============================================================================

template <int Dim> void TMesh<Dim>::ConstructTriangularElements(Array<ifstream *> &files,
								bool generateedges)
{

  mtype = Element::TRIANGLE;
  int t1, t2, t3, t4;
  *files[0] >> NumOfVertices >> t1 >> t2 >> t3;
  *files[1] >> NumOfElements >> t1 >> t2;  
  *files[2] >> t1 >> t2 >> t3 >> t4;
  *files[2] >> NumOfBdrElements >> t1;
  vertices.SetSize(NumOfVertices);
  element.SetSize(NumOfElements);
  bdrelement.SetSize(NumOfBdrElements);
  int NumEdges=0;
  Array<int *> edges(2);

  // Create vertices.
  for(int i=0; i<NumOfVertices; i++)
    {
      int t, bdrymarker;
      *files[0] >> t >> vertices[i](0) >> vertices[i](1) >> bdrymarker;
    }

  // Create elements  
  int *v = new int[3];
  int attr = -1;
  if(files.Size() == 3)
    {
      for(int i=0; i<NumOfElements; i++)
	{
	  int t;
	  *files[1] >> t >> v[0] >> v[1] >> v[2];
	  for (int j=0; j<3; j++) { v[j] = v[j] - 1; }
	  element[i] = new Triangle(v, attr); 
	}
    }
  else if(files.Size()==4)
    {
      int bm, num;
      *files[3] >> NumEdges >> bm;
      NumOfEdges = NumEdges;
      edges[0] = new int[NumEdges];
      edges[1] = new int[NumEdges];
      for(int i=0; i<NumEdges; i++)
	{
	  *files[3] >> num >> edges[0][i] >> edges[1][i] >> bm;
	}
 
      int *e = new int[3];
      attr = -1;
      for(int i=0; i<NumOfElements; i++)
	{
	  int t;
	  *files[1] >> t >> v[0] >> v[1] >> v[2];
	  for(int j=0; j<NumEdges; j++)
	    {
	      if((v[0]==edges[0][j] && v[1]==edges[1][j]) || 
		  (v[0]==edges[1][j] && v[1]==edges[0][j]))
		e[0] = j;
	      else if((v[1]==edges[0][j] && v[2]==edges[1][j]) || 
		  (v[1]==edges[1][j] && v[2]==edges[0][j]))
		e[1] = j;
	      else if((v[0]==edges[0][j] && v[2]==edges[1][j]) || 
		  (v[0]==edges[1][j] && v[2]==edges[0][j]))
		e[2] = j;
	      //  else
	      //cout << "Edge Problem Element" << endl;
	    }

	  for (int j=0; j<3; j++) { v[j] = v[j] - 1; }
	  element[i] = new Triangle(v, e, attr); 
	}
      delete []e;
    }
  delete []v;
  
  // Create boundary elements
  v = new int[2];
  if(files.Size() == 3)
    {
      NumOfBdrs = 0;
      for(int i=0; i<NumOfBdrElements; i++)
	{
	  int t;
	  *files[2] >> t >> v[0] >> v[1] >> attr;
	  if(attr>NumOfBdrs) NumOfBdrs = attr;
	  for (int j=0; j<2; j++) { v[j] = v[j] - 1; }
	  attr -= 1;
	  bdrelement[i] = new Segment(v, attr);
	}
    }
  else if(files.Size() == 4)
    {
      NumOfBdrs = 0;
      for(int i=0; i<NumOfBdrElements; i++)
	{
	  int t;
	  *files[2] >> t >> v[0] >> v[1] >> attr;
	  if(attr>NumOfBdrs) NumOfBdrs = attr;
	  int e;
	  for(int j=0; j<NumEdges; j++)
	    {
	      if((v[0]==edges[0][j] && v[1]==edges[1][j]) || 
		  (v[0]==edges[1][j] && v[1]==edges[0][j]))
		e = j;
	      // else
	      //	cout << "Edge Problem Boundary" << endl;
	    }

	  for (int j=0; j<2; j++) { v[j] = v[j] - 1; }
	  attr -= 1;
	  bdrelement[i] = new Segment(v, e, attr);
	}
      delete []edges[0];
      delete []edges[1];
    }

  delete []v;

  // set vertex to element map
  vert_to_ele.SetSize(NumOfVertices);
  
  int *count = new int[NumOfVertices];
  for(int i=0; i<NumOfVertices; i++)
    count[i] = 0;

  for(int i=0; i<NumOfElements; i++)
    {
      Array<int> ind;
      GetElementVertices(i, ind);
      for(int j=0; j<3; j++)
	count[ind[j]]++;
    }

  for(int i=0; i<NumOfVertices; i++)
    vert_to_ele[i] = new int[count[i]+1];

  for(int i=0; i<NumOfVertices; i++)
    vert_to_ele[i][0] = count[i];

  int *index = new int[NumOfVertices];
  for(int i=0; i<NumOfVertices; i++)
    index[i] = 1;

  for(int i=0; i<NumOfElements; i++)
    {
      Array<int> ind;
      GetElementVertices(i, ind);
      for(int j=0; j<3; j++)
	{
	vert_to_ele[ind[j]][index[ind[j]]] = i;
	index[ind[j]]++;
	}
    }
	 
  if (generateedges)
    {
      el_to_edge = new Table(GetNE(), 4);
      NumOfEdges = GetElementToEdgeTable(*el_to_edge, be_to_edge);
    }

  delete []index;
  delete []count;

}

//============================================================================
/*
template <int Dim> void TMesh<Dim>::Construct2DNURBSElements(Geometry &geom,
								  bool generateedges)
{
  
  NumOfElements = 0;
  NumOfBdrElements = 0;
  NumOfVertices = 0;
  int *size = new int[2*geom.GetNumOfPatches()];
  for (int i=0; i<geom.GetNumOfPatches(); i++)
    {
      int tempdeg = geom.GetPatchDegs(i,0);
      Degs[0] = tempdeg;
      tempdeg = geom.GetPatchDegs(i,1);
      Degs[1] = tempdeg;
      tempdeg = 1;
      for (int j=0; j<2; j++)
	{
	  int count = 1;
	  for (int k=1; k<geom.GetKnotSize(i,j); k++)
	    {
	      if (geom.GetKnot(i,j,k) > geom.GetKnot(i,j,k-1))
		{
		  count++;
		}
	    }
	  size[i*2+j] = count;
	}
      tempdeg = 1;
      for (int j=0; j<2; j++) {tempdeg *= size[i*2+j]-1;}
      NumOfElements += tempdeg;
      tempdeg = 1;
      for (int j=0; j<2; j++) {tempdeg *= size[i*2+j];}
      NumOfVertices += tempdeg;

      double temp = 0.0;
      for (int j=0; j<geom.GetNumOfFaces(); j++) 
	{ 
	  temp += 0.5 * size[i*2+j/2]*geom.GetFaceVal(i,j);
	}
      NumOfVertices -= temp;
      temp = 0.0;
      for (int j=0; j<geom.GetNumOfFaces(); j++) 
	{ 
	  temp += (size[2*i+j/2]-1)*(1.0-geom.GetFaceVal(i,j));
	}
      NumOfBdrElements += temp;
    }
  vertices.SetSize(NumOfVertices);
  element.SetSize(NumOfElements);
  bdrelement.SetSize(NumOfBdrElements);

  Array<Vector *> vert(geom.GetNumOfPatches());
  int index = 0;
  vert[0] = new Vector(size[0]*size[1]);
  for (int i=0; i<vert[0]->Size(); i++)
    {
      vert[0]->Elem(i) = index;
      index++;
    }
  for (int i=1; i<geom.GetNumOfPatches(); i++)
    {
      vert[i] = new Vector(size[2*i]*size[2*i+1]);
      for (int j=0; j<vert[i]->Size(); j++){vert[i]->Elem(j) = -1;}// set up initial vertex array to check for overwrites later below
      for (int j=0; j<geom.GetFaceSize(); j++)
	{
	  if (geom.GetConnectVal(i*geom.GetFaceSize()+j,0) > geom.GetConnectVal(i*geom.GetFaceSize()+j,1))
	    {
	      cout << "hit" << endl;
	      for (int k=0; k<geom.GetKnotSize(i,j/2); k++)
		{
		  cout << k*size[2*i]+(j-2)*(size[2*i]-1) << endl;
		  int iind = (1-j/2)*(k+j*size[2*i]*(size[2*i+1]-1)) + j/2*(k*size[2*i]+(j-2)*(size[2*i]-1)); // index in new patch to copy from old patch
		  int ref = geom.GetConnectVal(i*geom.GetFaceSize()+j,3); // local face index of old patch
		  int refpatch = geom.GetConnectVal(i*geom.GetFaceSize()+j,1); //old patch
		  int irefind = (1-ref/2)*(k+ref*size[2*refpatch]*(size[2*refpatch+1]-1)) + ref/2*(k*size[2*refpatch]+(ref-2)*(size[2*refpatch]-1)); // index from old patch
		  double c = vert[refpatch]->Elem(irefind);
		  vert[i]->Elem(iind) = c;
		}
	    }
	}
      for (int j=0; j<size[2*i] * size[2*i+1]; j++)
	{
	  if (vert[i]->Elem(j)<0) //checking for overwrites
	    {
	      vert[i]->Elem(j) = index;
	      index++;
	    }
	}
    }
  cout << endl;
  for (int i=0; i<geom.GetNumOfPatches(); i++)
    {
      for (int j=0; j<vert[i]->Size(); j++)
	{
	  cout << vert[i]->Elem(j) << "  ";
	}
      cout << endl;
    }
  for (int i=0; i<vert.Size(); i++)
    {
      delete vert[i];
    }
  delete []size;
}

//============================================================================

template <int Dim> void TMesh<Dim>::ConstructSurfaceNURBSElements(Geometry &geom,
								  bool generateedges)
{
  
  NumOfVertices = 0;
  NumOfElements = 0;
  NumOfBdrElements = 0;
  for (int i=0; i<geom.GetNumOfPatches(); i++)
    {
      int tempdeg = geom.GetPatchDegs(i,0);
      Degs[0] = tempdeg;
      tempdeg = geom.GetPatchDegs(i,1);
      Degs[1] = tempdeg;
      tempdeg = 1;
      for (int j=0; j<2; j++) {tempdeg *= geom.GetKnotSize(i,j);}
      NumOfElements = tempdeg;
      double temp = NumOfElements;
      tempdeg = 1;
      for (int j=0; j<2; j++) {tempdeg *= geom.GetKnotSize(i,j)+1;}
      NumOfVertices += tempdeg;
      temp = 0.0;
      for (int j=0; j<geom.GetNumOfFaces(); j++) 
	{ 
	  temp += 0.5 * geom.GetKnotSize(i,j)*geom.GetFaceVal(i,j);
	}
      NumOfElements -= temp;

      tempdeg = 1;
      for (int j=0; j<2; j++) {tempdeg *= geom.GetKnotSize(i,j)+1;}
      NumOfVertices += tempdeg;
      temp = 0.0;
      for (int j=0; j<geom.GetNumOfFaces(); j++) 
	{ 
	  temp += geom.GetKnotSize(i,j)*(1.0-geom.GetFaceVal(i,j));
	}
      NumOfBdrElements += temp;
    }
  vertices.SetSize(NumOfVertices);
  element.SetSize(NumOfElements);
  bdrelement.SetSize(NumOfBdrElements);
  
}

//============================================================================

template <int Dim> void TMesh<Dim>::Construct3DNURBSElements(Geometry &geom,
								  bool generateedges)
{
  
  NumOfVertices = 0;
  NumOfElements = 0;
  NumOfBdrElements = 0;
  for (int i=0; i<geom.GetNumOfPatches(); i++)
    {
      int tempdeg = geom.GetPatchDegs(i,0);
      Degs[0] = tempdeg;
      tempdeg = geom.GetPatchDegs(i,1);
      Degs[1] = tempdeg;
      tempdeg = 1;
      tempdeg = geom.GetPatchDegs(i,2);
      Degs[2] = tempdeg;
      for (int j=0; j<3; j++) {tempdeg *= geom.GetKnotSize(i,j);}
      NumOfElements = tempdeg;
      double temp = NumOfElements;
      for (int j=0; j<geom.GetNumOfFaces(); j++){temp -= 0.5*geom.GetFaceVal(i,j)*geom.GetKnotSize(i,j);}
      NumOfVertices = temp;
    }
  vertices.SetSize(NumOfVertices);
  element.SetSize(NumOfElements);
  bdrelement.SetSize(NumOfBdrElements);
  
}
*/
//============================================================================

template <int Dim> void TMesh<Dim>::ConstructTetrahedralElements(Array<ifstream *> &files,
								  bool generateedges)
{
  mtype = Element::TETRAHEDRAL;
  int temp1, temp2, temp3;
  *files[0] >> NumOfVertices >> temp1 >> temp2 >> temp3;
  *files[1] >> NumOfElements >> temp1 >> temp2;
  *files[2] >> NumOfBdrElements >> temp1;

  vertices.SetSize(NumOfVertices);
  element.SetSize(NumOfElements);
  bdrelement.SetSize(NumOfBdrElements);
      
      
  // Create vertices.
  for(int j=0; j<NumOfVertices; j++)
    {
      int temp;
      *files[0] >> temp >> vertices[j](0) >> vertices[j](1) >> vertices[j](2);
    }
  
  int v[4];//the 4 indeces of the elements
  int attr = 0;
  // Create elements
  for(int j=0; j<NumOfElements; j++)
    {
      int temp;
      *files[1] >> temp >> v[0] >> v[1] >> v[2] >> v[3];

      for(int i=0; i<4; i++)
	{
	  v[i]-=1;
	}

      element[j] = new Tetrahedral(v, attr); 
    }

  //boundary elements
  int vv[3];
  NumOfBdrs = 0;
  for(int j=0; j<NumOfBdrElements; j++)
    {
      int temp, attribute;
      *files[2] >> temp >> vv[0] >> vv[1] >> vv[2] >> attribute;
      if(attribute>NumOfBdrs) { NumOfBdrs = attribute; }      

      for(int i=0; i<3; i++)
	{
	  vv[i]--;
	}

      attribute--;
      bdrelement[j] = new Triangle(vv, attribute);
    }
}

//=============================================================================


template <int Dim> TMesh<Dim>::TMesh(Array<double> &xcoord)
{
  NumOfVertices = xcoord.Size();
  NumOfElements = NumOfVertices - 1;
  NumOfBdrElements = 2;
  vertices.SetSize(NumOfVertices);
  element.SetSize(NumOfElements);
  bdrelement.SetSize(NumOfBdrElements);

  // Create vertices.
  //double cx = xl;
  for (int i=0; i<NumOfVertices; i++) { vertices[i](0) = xcoord[i]; }

  // Create element.
  int *v = new int[2];
  for (int i=0; i<NumOfElements; i++)
    {
      v[0] = i;
      v[1] = i+1;
      int attr = -1;
      element[i] = new Segment(v, attr);
    }
  delete []v;

  int attr = -1;
  v = new int[1];
  v[0] = 0;
  bdrelement[0] = new Point(v, attr);
  v[0] = NumOfElements;
  bdrelement[1] = new Point(v, attr);
}

//============================================================================

template <int Dim> TMesh<Dim>::TMesh(Array<int> &n, Array<double> &L, Element::Type type,
				      bool generateedges)
{
  if (type==Element::QUADRILATERAL)
      ConstructRectangularElements(n, L, generateedges);

  else if (type==Element::TRIANGLE)
    ConstructTriangularElements(n, L, generateedges);

  else if (type==Element::TETRAHEDRAL)
    ConstructTetrahedralElements(n, L, generateedges);
}

//============================================================================

template <int Dim> TMesh<Dim>::TMesh(Array<ifstream *> &files, Element::Type type,
				     bool generateedges)
{
  if (type==Element::QUADRILATERAL)
    {
      ConstructRectangularElements(files, generateedges);
    }
  else if (type==Element::TRIANGLE)
    {
      ConstructTriangularElements(files, generateedges);
    }
  else if (type==Element::TETRAHEDRAL)
    {
      ConstructTetrahedralElements(files, generateedges);
    }
}


//============================================================================
/*
template <int Dim> TMesh<Dim>::TMesh(Geometry &Geom, Element::Type type,
				     bool generateedges)
{
  if (type==Element::TwoDNURBS)
    {
      Construct2DNURBSElements(Geom, generateedges);
    }
 if (type==Element::SurfaceNURBS)
    {
      ConstructSurfaceNURBSElements(Geom, generateedges);
    }
 if (type==Element::ThreeDNURBS)
    {
      Construct3DNURBSElements(Geom, generateedges);
    }
}
*/

//============================================================================

template <int Dim> TMesh<Dim>::~TMesh()
{
  for(int k=0; k<BdrytoGlobalIndex.Size(); k++)
    delete []BdrytoGlobalIndex[k];

  if (Dim==2)
    {
      if (el_to_edge != NULL) { delete el_to_edge; }
    }
  
}

template class TMesh<1>;
template class TMesh<2>;
template class TMesh<3>;
