#include "mpi_mesh_header.h"

//============================================================================

Mesh_p::~Mesh_p()
{
  for (int i=0; i<NumOfElements; i++) { delete element[i]; }

  for (int i=0; i<NumOfBdrElements; i++) { delete bdrelement[i]; }
}

//============================================================================

void Mesh_p::Print(ofstream &out, int type) const
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

void Mesh_p::ParaviewPrint(int type) const
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

void Mesh_p::ParaviewPrint(const Vector &solution, int type) const
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

void Mesh_p::Print(ofstream &out, const Vector &solution) const
{
  out << "solution u -size=" << NumOfVertices << " -type=nodal -components=1" << endl;
  
  for (int i=0; i<NumOfVertices; i++)
    out << solution(i) << endl;
}


//============================================================================

template <int Dim> void TMesh_p<Dim>::GetElementVerticesCoord(int i, Array<double *> &coord) const
{
  int nv = element[i]->GetNVertices();
  Array<int> v;
  element[i]->GetVertices(v);
  
  for (int k=0; k<Dim; k++)
    for (int l=0; l<nv; l++)
      coord[k][l] = vertices[v[l]](k);
}

//============================================================================

template <int Dim> void TMesh_p<Dim>::GetBdrElementVerticesCoord(int i, Array<double *> &coord) const
{
  int nv = bdrelement[i]->GetNVertices();
  Array<int> v;
  bdrelement[i]->GetVertices(v);
  for (int k=0; k<Dim; k++)
    for (int l=0; l<nv; l++)
      coord[k][l] = vertices[v[l]](k);
}

//============================================================================

template <int Dim> TMesh_p<Dim>::TMesh_p(int n, double xl, double xr)
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

template <int Dim> void TMesh_p<Dim>::ConstructRectangularElements(Array<int> &n, Array<double> &L)
{
    //array n and array L must be allocated before this function is used
    //array n is [n_x, n_y], array L is [x1, x2, y1, y2]
    //          [ n[0], n[1] ]          [ L[0], L[1], L[2], L[3] ]
    
    //the origin is the left bottom
    
    //the index of the vertices follows the matrix in C, and x-direction is the 
    //row, and y-direction is the column
    
    //the index of the elements also follows the matrix in C, 
    //the vertex is counter-clockwise in each element
    
    //the index of the boundary elements also follows the matrix in C 
    
    double hx = (L[1] - L[0]) / ( (double) n[0] );
    double hy = (L[3] - L[2]) / ( (double) n[1] );
    NumOfVertices = (n[0] + 1) * (n[1] + 1);
    NumOfElements = n[0] * n[1];
    NumOfBdrElements = 2 * (n[0] + n[1]);
    vertices.SetSize(NumOfVertices);
    element.SetSize(NumOfElements);
    bdrelement.SetSize(NumOfBdrElements);
    
    //vertices
    for (int i=0; i<NumOfVertices; i++)
    {
        int row = i / (n[0]+1);
        int col = i - row * (n[0]+1);
        
        vertices[i](0) = col * hx;
        vertices[i](1) = row * hy;
    }
    
    //elements
    int v[4];//the 4 indeces of the elements
    int attr = -1;
    for (int i=0; i<NumOfElements; i++)
    {
        int row = i / n[0];
        v[0] = i + row;
        v[1] = v[0] + 1;
        v[2] = v[1] + (n[0] + 1);
        v[3] = v[2] - 1;
        
        element[i] = new Quadrilateral(v, attr); //Question here: element[i] points to the
        //the base class, but the new Quadrilateral points to the dirived class
        //So, is it correct here?
    }
    
    //boundary elements
    int vv[2];
    int attr1 = -1;
    for (int i=0; i<NumOfBdrElements; i++)
    {
        if (i < n[0] )
        {
            vv[0] = i;
            vv[1] = i + 1;
            
            bdrelement[i] = new Segment(vv, attr1);
        }
        else if (i >= n[0] && i < n[0]+n[1])
        {
            vv[0] = 0 + (i - n[0]) * (n[0] + 1);
            vv[1] = vv[0] + (n[0] + 1);
            
            bdrelement[i] = new Segment(vv, attr1);
        }
        else if (i >= n[0]+n[1] && i < n[0]+2*n[1])
        {
            vv[0] = n[0] + (i - (n[0]+n[1])) * (n[0] + 1);
            vv[1] = vv[0] + (n[0] + 1);
            bdrelement[i] = new Segment(vv, attr1);
        }
        else
        {
            vv[0] = n[1]*(n[0]+1) + ( i - (n[0]+2*n[1]) );
            vv[1] = vv[0] + 1;
            
            bdrelement[i] = new Segment(vv, attr1);
        }
    }
    
    
}

//=============================================================================

template <int Dim> void TMesh_p<Dim>::ConstructTriangularElements(Array<int> &n, Array<double> &L)
{
  // Array n = (n[0], n[1]) = (n_x, n_y)
  // Array L = (L[0], L[1], L[2], L[3]) =  (x1, x2, y1, y2)

  double dx = (L[1] - L[0]) / ( (double) n[0] );
  double dy = (L[3] - L[2]) / ( (double) n[1] );
  NumOfBdrs = 4;
  NumOfVertices = (n[0] + 1) * (n[1] + 1);
  NumOfElements = 2 * n[0] * n[1];
  NumOfBdrElements = 2 * (n[0] + n[1]);
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
  
  // Create elements
  int *v = new int[3];
  int attr = -1;
  for (int i=0; i<NumOfElements; i++)
    {
      jthrow = i / ( 2*n[0] );
      if(i%2 == 0) // Lower triangle
	{
	  v[0] = i/2 + jthrow;
	  v[1] = v[0] + 1;
	  v[2] = v[0] + (n[0] + 1);
	}
      else // Upper triangle
	{
	  v[0] = i/2 + 1 + jthrow;
	  v[1] = v[0] + (n[0] + 1);
	  v[2] = v[1] - 1;
	}
     
      element[i] = new Triangle(v, attr);
    }
  delete []v;
  
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
}

//=============================================================================

template <int Dim> void TMesh_p<Dim>::ConstructTetrahedralElements(Array<int> &n, Array<double> &L)
  {

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

template <int Dim> void TMesh_p<Dim>::ConstructRectangularElements(Array<ifstream *> &files)
{
    //we only consider the four-point polygon to set up the attr of the boundary element
    //and to set up the attr of the boundary element, here we consider the simple case [0,1] by [0,1], four boundaries
    
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

template <int Dim> void TMesh_p<Dim>::ConstructTriangularElements(Array<ifstream *> &files)
{
  int t1, t2, t3, t4;
  *files[0] >> NumOfVertices >> t1 >> t2 >> t3;
  *files[1] >> NumOfElements >> t1 >> t2;  
  *files[2] >> t1 >> t2 >> t3 >> t4;
  *files[2] >> NumOfBdrElements >> t1;
  vertices.SetSize(NumOfVertices);
  element.SetSize(NumOfElements);
  bdrelement.SetSize(NumOfBdrElements);
          
  // Create vertices.
  for(int i=0; i<NumOfVertices; i++)
    {
      int t, bdrymarker;
      *files[0] >> t >> vertices[i](0) >> vertices[i](1) >> bdrymarker;
    }

  // Create elements  
  int *v = new int[3];
  int attr = -1;
  for(int i=0; i<NumOfElements; i++)
    {
      int t;
      *files[1] >> t >> v[0] >> v[1] >> v[2];
      for (int j=0; j<3; j++) { v[j] = v[j] - 1; }
      element[i] = new Triangle(v, attr); 
    }
  delete []v;
  
  // Create boundary elements
  v = new int[2];
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
  delete []v;
}

//============================================================================

template <int Dim> void TMesh_p<Dim>::ConstructTetrahedralElements(Array<ifstream *> &files)
{
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


template <int Dim> TMesh_p<Dim>::TMesh_p(Array<double> &xcoord)
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

template <int Dim> TMesh_p<Dim>::TMesh_p(Array<int> &n, Array<double> &L, Element::Type type)
{
  if (type==Element::QUADRILATERAL)
    ConstructRectangularElements(n, L);

  else if (type==Element::TRIANGLE)
    ConstructTriangularElements(n, L);

  else if (type==Element::TETRAHEDRAL)
    ConstructTetrahedralElements(n, L);
}

//============================================================================

template <int Dim> TMesh_p<Dim>::TMesh_p(Array<ifstream *> &files, Element::Type type)
{
  if (type==Element::QUADRILATERAL)
    ConstructRectangularElements(files);

  else if (type==Element::TRIANGLE)
    ConstructTriangularElements(files);

  else if (type==Element::TETRAHEDRAL)
    ConstructTetrahedralElements(files);
}

//============================================================================

template <int Dim> TMesh_p<Dim>::~TMesh_p()
{
  ;
}

template class TMesh_p<1>;
template class TMesh_p<2>;
template class TMesh_p<3>;
