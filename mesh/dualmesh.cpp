#include "mesh_header.h"

DualMesh::DualMesh(Mesh *_mesh, int _dualmeshorder, int _dualmeshtype)
{
  mesh = _mesh;
  dualmeshorder = _dualmeshorder;
  dualmeshtype = _dualmeshtype; 
  if(dualmeshorder==1)
    dof = mesh->GetNV();
  else if(dualmeshorder==2)
    dof = mesh->GetNV() + mesh->GetNEdges();
  else if(dualmeshorder==3)
    dof = mesh->GetNV() + 2*mesh->GetNEdges() + mesh->GetNE();
  else
    {
      cout << "DualMesh order error: " << dualmeshorder << endl;
      exit(1);
    }
}

//=============================================================================

void DualMesh::GetElemDualInfo(int i, Array<double *> &coord, Array<double *> &normals,
			       Array<double> &elengths, Array<double> &areas)
{
  int n = (mesh->GetMType()==Element::TRIANGLE) ? 3 : 4;  

  Array<int> ind(n);
  Array<double *> v(n);
  mesh->GetElementVertices (i, ind);
  for (int k=0; k<n; k++) { v[k] = mesh->GetVertex(ind[k]); }

  if(n==3)
    {
      double *v0 = v[0];
      double *v1 = v[1];
      double *v2 = v[2];
      if(dualmeshorder==1)
	{
	  for(int k=0; k<2; k++)
	    {
	      coord[0][k] = (v1[k]+v0[k]) / 2.0;
	      coord[1][k] = (v2[k]+v1[k]) / 2.0;
	      coord[2][k] = (v2[k]+v0[k]) / 2.0;
	      coord[3][k] = (v0[k]+v1[k]+v2[k]) / 3.0;
	    }

	  for(int k=0; k<3; k++)
	    {
	      elengths[k] = sqrt((coord[k][0]-coord[3][0])*(coord[k][0]-coord[3][0]) + 
				 (coord[k][1]-coord[3][1])*(coord[k][1]-coord[3][1]));
	    }

	  for(int k=0; k<3; k++)
	    {
	      double x, y;
	      x = coord[3][0]-coord[k][0];
	      y = coord[3][1]-coord[k][1];
	      if(k==0)
		{
		  normals[k][0] = -y/elengths[k];
		  normals[k][1] = x/elengths[k];
		  if (normals[k][0]*(v1[0]-v0[0])+normals[k][1]*(v1[1]-v0[1]) < 0)
		    {
		      normals[k][0] = -normals[k][0];
		      normals[k][1] = -normals[k][1];
		    }		  
		}
	      else if(k==1)
		{
		  normals[k][0] = y/elengths[k];
		  normals[k][1] = -x/elengths[k];
		  if (normals[k][0]*(v2[0]-v1[0])+normals[k][1]*(v2[1]-v1[1]) < 0)
		    {
		      normals[k][0] = -normals[k][0];
		      normals[k][1] = -normals[k][1];
		    }		  
		}
	      else if(k==2)
		{
		  normals[k][0] = y/elengths[k];
		  normals[k][1] = -x/elengths[k];
		  if (normals[k][0]*(v0[0]-v2[0])+normals[k][1]*(v0[1]-v2[1]) < 0)
		    {
		      normals[k][0] = -normals[k][0];
		      normals[k][1] = -normals[k][1];
		    }
		}
	    }

	  double *p, *q;
	  p = new double[2];
	  q = new double[2];
	  p[0] = coord[2][0]-coord[0][0];
	  p[1] = coord[2][1]-coord[0][1];
	  q[0] = coord[3][0]-v0[0];
	  q[1] = coord[3][1]-v0[1];
	  areas[0] = 0.5*fabs(p[0]*q[1]-p[1]*q[0]);
	  areas[1] = areas[0]; 
	  areas[2] = areas[1];
	  delete []p;
	  delete []q;
	}
      else if(dualmeshorder==2)
	{
	  if(dualmeshtype==0)
	    {
	      for(int k=0; k<2; k++)
		{
		  coord[0][k] = v0[k]*3.0/4.0 + v1[k]/4.0; 
		  coord[1][k] = v0[k]/2.0 + v1[k]/4.0 + v2[k]/4.0; 
		  coord[2][k] = v0[k]*3.0/4.0 + v2[k]/4.0;
		  coord[3][k] = v0[k]*2.0/3.0 + v1[k]/6.0 + v2[k]/6.0;
		  
		  coord[4][k] = v1[k]*3.0/4.0 + v2[k]/4.0; 
		  coord[5][k] = v1[k]/2.0 + v0[k]/4.0 + v2[k]/4.0; 
		  coord[6][k] = v1[k]*3.0/4.0 + v0[k]/4.0;
		  coord[7][k] = v1[k]*2.0/3.0 + v0[k]/6.0 + v2[k]/6.0;
		  
		  coord[8][k] = v2[k]*3.0/4.0 + v0[k]/4.0; 
		  coord[9][k] = v2[k]/2.0 + v1[k]/4.0 + v0[k]/4.0; 
		  coord[10][k] = v2[k]*3.0/4.0 + v1[k]/4.0;
		  coord[11][k] = v2[k]*2.0/3.0 + v1[k]/6.0 + v0[k]/6.0;
	  	  
		  coord[12][k] = (v0[k] + v1[k] + v2[k])/3.0;
		}
	      
	      for(int k=0; k<3; k++)
		{
		  elengths[k] = sqrt((coord[k][0]-coord[3][0])*(coord[k][0]-coord[3][0]) + 
				     (coord[k][1]-coord[3][1])*(coord[k][1]-coord[3][1]));
		}

	      double x, y;
	      x = coord[3][0]-coord[0][0];
	      y = coord[3][1]-coord[0][1];
	      normals[0][0] = y/elengths[0];
	      normals[0][1] = -x/elengths[0];
	      if (normals[0][0]*(v1[0]-v0[0])+normals[0][1]*(v1[1]-v0[1]) < 0)
		{
		  normals[0][0] = -normals[0][0];
		  normals[0][1] = -normals[0][1];
		}
      
	      x = coord[3][0]-coord[1][0];
	      y = coord[3][1]-coord[1][1];
	      normals[1][0] = y/elengths[1];
	      normals[1][1] = -x/elengths[1];
	      if (normals[1][0]*(v2[0]-v1[0])+normals[1][1]*(v2[1]-v1[1]) < 0)
		{
		  normals[1][0] = -normals[1][0];
		  normals[1][1] = -normals[1][1];
		}
      
	      x = coord[3][0]-coord[2][0];
	      y = coord[3][1]-coord[2][1];
	      normals[2][0] = y/elengths[2];
	      normals[2][1] = -x/elengths[2];
	      if (normals[2][0]*(v0[0]-v2[0])+normals[2][1]*(v0[1]-v2[1]) < 0)
		{
		  normals[2][0] = -normals[2][0];
		  normals[2][1] = -normals[2][1];
		}
      
	      double *p, *q;
	      p = new double[2];
	      q = new double[2];
	      p[0] = coord[2][0]-coord[0][0];
	      p[1] = coord[2][1]-coord[0][1];
	      q[0] = coord[3][0]-v0[0];
	      q[1] = coord[3][1]-v0[1];
	      areas[0] = 0.5*fabs(p[0]*q[1]-p[1]*q[0]);
	      delete []p;
	      delete []q;
	    }
	  
	  else if(dualmeshtype==1)
	    {
	      for(int k=0; k<2; k++)
		{
		  coord[0][k] = v0[k]*5.0/6.0 + v1[k]/6.0; 
		  coord[1][k] = v0[k]/2.0 + v1[k]/4.0 + v2[k]/4.0; 
		  coord[2][k] = v0[k]*5.0/6.0 + v2[k]/6.0;
		  coord[3][k] = v0[k]*2.0/3.0 + v1[k]/6.0 + v2[k]/6.0;
		  
		  coord[4][k] = v1[k]*5.0/6.0 + v2[k]/6.0; 
		  coord[5][k] = v1[k]/2.0 + v0[k]/4.0 + v2[k]/4.0; 
		  coord[6][k] = v1[k]*5.0/6.0 + v0[k]/6.0;
		  coord[7][k] = v1[k]*2.0/3.0 + v0[k]/6.0 + v2[k]/6.0;
		  
		  coord[8][k] = v2[k]*5.0/6.0 + v0[k]/6.0; 
		  coord[9][k] = v2[k]/2.0 + v1[k]/4.0 + v0[k]/4.0; 
		  coord[10][k] = v2[k]*5.0/6.0 + v1[k]/6.0;
		  coord[11][k] = v2[k]*2.0/3.0 + v1[k]/6.0 + v0[k]/6.0;
	  	  
		  coord[12][k] = (v0[k] + v1[k] + v2[k])/3.0;
		}
	      
	      for(int j=0; j<3; j++)
		{
		  elengths[j] = sqrt((coord[j][0]-coord[3][0])*(coord[j][0]-coord[3][0]) + 
				     (coord[j][1]-coord[3][1])*(coord[j][1]-coord[3][1]));
	  
		  int k=j+4;
		  elengths[k-1] = sqrt((coord[k][0]-coord[7][0])*(coord[k][0]-coord[7][0]) + 
				       (coord[k][1]-coord[7][1])*(coord[k][1]-coord[7][1]));
	  
		  k=j+8;
		  elengths[k-2] = sqrt((coord[k][0]-coord[11][0])*(coord[k][0]-coord[11][0]) + 
				       (coord[k][1]-coord[11][1])*(coord[k][1]-coord[11][1]));
		}
	      elengths[9] = sqrt((coord[9][0]-coord[12][0])*(coord[9][0]-coord[12][0]) + 
				 (coord[9][1]-coord[12][1])*(coord[9][1]-coord[12][1]));
      
	      elengths[10] = sqrt((coord[1][0]-coord[12][0])*(coord[1][0]-coord[12][0]) + 
				  (coord[1][1]-coord[12][1])*(coord[1][1]-coord[12][1]));
      
	      elengths[11] = sqrt((coord[5][0]-coord[12][0])*(coord[5][0]-coord[12][0]) + 
				  (coord[5][1]-coord[12][1])*(coord[5][1]-coord[12][1]));
      
	      double x, y;
	      x = coord[3][0]-coord[0][0];
	      y = coord[3][1]-coord[0][1];
	      normals[0][0] = y/elengths[0];
	      normals[0][1] = -x/elengths[0];
	      if (normals[0][0]*(v1[0]-v0[0])+normals[0][1]*(v1[1]-v0[1]) < 0)
		{
		  normals[0][0] = -normals[0][0];
		  normals[0][1] = -normals[0][1];
		}
      
	      x = coord[3][0]-coord[1][0];
	      y = coord[3][1]-coord[1][1];
	      normals[1][0] = y/elengths[1];
	      normals[1][1] = -x/elengths[1];
	      if (normals[1][0]*(v2[0]-v1[0])+normals[1][1]*(v2[1]-v1[1]) < 0)
		{
		  normals[1][0] = -normals[1][0];
		  normals[1][1] = -normals[1][1];
		}
      
	      x = coord[3][0]-coord[2][0];
	      y = coord[3][1]-coord[2][1];
	      normals[2][0] = y/elengths[2];
	      normals[2][1] = -x/elengths[2];
	      if (normals[2][0]*(v0[0]-v2[0])+normals[2][1]*(v0[1]-v2[1]) < 0)
		{
		  normals[2][0] = -normals[2][0];
		  normals[2][1] = -normals[2][1];
		}

	      x = coord[7][0]-coord[4][0];
	      y = coord[7][1]-coord[4][1];
	      normals[3][0] = y/elengths[3];
	      normals[3][1] = -x/elengths[3];
	      if (normals[3][0]*(v2[0]-v1[0])+normals[3][1]*(v2[1]-v1[1]) < 0)
		{
		  normals[3][0] = -normals[3][0];
		  normals[3][1] = -normals[3][1];
		}
      
	      x = coord[7][0]-coord[5][0];
	      y = coord[7][1]-coord[5][1];
	      normals[4][0] = y/elengths[4];
	      normals[4][1] = -x/elengths[4];
	      if (normals[4][0]*(v0[0]-v2[0])+normals[4][1]*(v0[1]-v2[1]) < 0)
		{
		  normals[4][0] = -normals[4][0];
		  normals[4][1] = -normals[4][1];
		}
      
	      x = coord[7][0]-coord[6][0];
	      y = coord[7][1]-coord[6][1];
	      normals[5][0] = y/elengths[5];
	      normals[5][1] = -x/elengths[5];
	      if (normals[5][0]*(v1[0]-v0[0])+normals[5][1]*(v1[1]-v0[1]) < 0)
		{
		  normals[5][0] = -normals[5][0];
		  normals[5][1] = -normals[5][1];
		}

	      x = coord[11][0]-coord[8][0];
	      y = coord[11][1]-coord[8][1];
	      normals[6][0] = y/elengths[6];
	      normals[6][1] = -x/elengths[6];
	      if (normals[6][0]*(v0[0]-v2[0])+normals[6][1]*(v0[1]-v2[1]) < 0)
		{
		  normals[6][0] = -normals[6][0];
		  normals[6][1] = -normals[6][1];
		}
      
	      x = coord[11][0]-coord[9][0];
	      y = coord[11][1]-coord[9][1];
	      normals[7][0] = y/elengths[7];
	      normals[7][1] = -x/elengths[7];
	      if (normals[7][0]*(v1[0]-v0[0])+normals[7][1]*(v1[1]-v0[1]) < 0)
		{
		  normals[7][0] = -normals[7][0];
		  normals[7][1] = -normals[7][1];
		}
      
	      x = coord[11][0]-coord[10][0];
	      y = coord[11][1]-coord[10][1];
	      normals[8][0] = y/elengths[8];
	      normals[8][1] = -x/elengths[8];
	      if (normals[8][0]*(v2[0]-v1[0])+normals[8][1]*(v2[1]-v1[1]) < 0)
		{
		  normals[8][0] = -normals[8][0];
		  normals[8][1] = -normals[8][1];
		}

	      x = coord[12][0]-coord[9][0];
	      y = coord[12][1]-coord[9][1];
	      normals[9][0] = y/elengths[9];
	      normals[9][1] = -x/elengths[9];
	      if (normals[9][0]*(v1[0]-v0[0])+normals[9][1]*(v1[1]-v0[1]) < 0)
		{
		  normals[9][0] = -normals[9][0];
		  normals[9][1] = -normals[9][1];
		}
      
	      x = coord[12][0]-coord[1][0];
	      y = coord[12][1]-coord[1][1];
	      normals[10][0] = y/elengths[10];
	      normals[10][1] = -x/elengths[10];
	      if (normals[10][0]*(v2[0]-v1[0])+normals[10][1]*(v2[1]-v1[1]) < 0)
		{
		  normals[10][0] = -normals[10][0];
		  normals[10][1] = -normals[10][1];
		}
      
	      x = coord[12][0]-coord[5][0];
	      y = coord[12][1]-coord[5][1];
	      normals[11][0] = y/elengths[11];
	      normals[11][1] = -x/elengths[11];
	      if (normals[11][0]*(v0[0]-v2[0])+normals[11][1]*(v0[1]-v2[1]) < 0)
		{
		  normals[11][0] = -normals[11][0];
		  normals[11][1] = -normals[11][1];
		}
      
	      double *p, *q;
	      p = new double[2];
	      q = new double[2];
	      p[0] = coord[2][0]-coord[0][0];
	      p[1] = coord[2][1]-coord[0][1];
	      q[0] = coord[3][0]-v0[0];
	      q[1] = coord[3][1]-v0[1];
	      areas[0] = 0.5*fabs(p[0]*q[1]-p[1]*q[0]);

	      delete []p;
	      delete []q;
	    }
	  else
	    cout << "Quadratic DualMesh is not constructed for this type: " << dualmeshtype << endl;	  
	}
      else if(dualmeshorder==3)
	{
	  for(int k=0; k<2; k++)
	    {
	      coord[0][k] = v0[k]*5.0/6.0 + v1[k]/6.0;   
	      coord[1][k] = v0[k]/2.0 + v1[k]/2.0;
	      coord[2][k] = v1[k]*5.0/6.0 + v0[k]/6.0;
	      
	      coord[3][k] = v1[k]*5.0/6.0 + v2[k]/6.0;   
	      coord[4][k] = v1[k]/2.0 + v2[k]/2.0;
	      coord[5][k] = v2[k]*5.0/6.0 + v1[k]/6.0;
	      
	      coord[6][k] = v2[k]*5.0/6.0 + v0[k]/6.0;   
	      coord[7][k] = v0[k]/2.0 + v2[k]/2.0;
	      coord[8][k] = v0[k]*5.0/6.0 + v2[k]/6.0;
	      
	      coord[9][k] = v0[k]*7.0/9.0 + v1[k]/9.0 + v2[k]/9.0;
	      coord[10][k] = v1[k]*7.0/9.0 + v0[k]/9.0 + v2[k]/9.0;
	      coord[11][k] = v2[k]*7.0/9.0 + v1[k]/9.0 + v0[k]/9.0;
	      
	      coord[12][k] = v0[k]*5.0/9.0 + v1[k]*2.0/9.0 + v2[k]*2.0/9.0;
	      coord[14][k] = v1[k]*5.0/9.0 + v0[k]*2.0/9.0 + v2[k]*2.0/9.0;
	      coord[16][k] = v2[k]*5.0/9.0 + v1[k]*2.0/9.0 + v0[k]*2.0/9.0;
	      
	      coord[13][k] = v0[k]*4.0/9.0 + v1[k]*4.0/9.0 + v2[k]/9.0;
	      coord[15][k] = v2[k]*4.0/9.0 + v1[k]*4.0/9.0 + v0[k]/9.0;
	      coord[17][k] = v0[k]*4.0/9.0 + v2[k]*4.0/9.0 + v1[k]/9.0;
	      
	      coord[18][k] = v0[k]/3.0 + v1[k]/3.0 + v2[k]/3.0;
	    }

	  elengths[0] = sqrt((coord[9][0]-coord[0][0])*(coord[9][0]-coord[0][0]) + 
			     (coord[9][1]-coord[0][1])*(coord[9][1]-coord[0][1]));
	  elengths[1] = sqrt((coord[9][0]-coord[8][0])*(coord[9][0]-coord[8][0]) + 
			     (coord[9][1]-coord[8][1])*(coord[9][1]-coord[8][1]));
	  
	  elengths[2] = sqrt((coord[10][0]-coord[3][0])*(coord[10][0]-coord[3][0]) + 
			     (coord[10][1]-coord[3][1])*(coord[10][1]-coord[3][1]));
	  elengths[3] = sqrt((coord[10][0]-coord[2][0])*(coord[10][0]-coord[2][0]) + 
			     (coord[10][1]-coord[2][1])*(coord[10][1]-coord[2][1]));
	  
	  elengths[4] = sqrt((coord[11][0]-coord[6][0])*(coord[11][0]-coord[6][0]) + 
			     (coord[11][1]-coord[6][1])*(coord[11][1]-coord[6][1]));
	  elengths[5] = sqrt((coord[11][0]-coord[5][0])*(coord[11][0]-coord[5][0]) + 
			     (coord[11][1]-coord[5][1])*(coord[11][1]-coord[5][1]));
	  
	  elengths[6] = sqrt((coord[13][0]-coord[1][0])*(coord[13][0]-coord[1][0]) + 
			     (coord[13][1]-coord[1][1])*(coord[13][1]-coord[1][1]));
	  elengths[7] = sqrt((coord[15][0]-coord[4][0])*(coord[15][0]-coord[4][0]) + 
			     (coord[15][1]-coord[4][1])*(coord[15][1]-coord[4][1]));
	  elengths[8] = sqrt((coord[17][0]-coord[7][0])*(coord[17][0]-coord[7][0]) + 
			     (coord[17][1]-coord[7][1])*(coord[17][1]-coord[7][1]));
	  
	  elengths[9] = sqrt((coord[12][0]-coord[9][0])*(coord[12][0]-coord[9][0]) + 
			     (coord[12][1]-coord[9][1])*(coord[12][1]-coord[9][1]));
	  elengths[12] = sqrt((coord[12][0]-coord[13][0])*(coord[12][0]-coord[13][0]) + 
			      (coord[12][1]-coord[13][1])*(coord[12][1]-coord[13][1]));
	  elengths[17] = sqrt((coord[12][0]-coord[17][0])*(coord[12][0]-coord[17][0]) + 
			      (coord[12][1]-coord[17][1])*(coord[12][1]-coord[17][1]));
	  
	  elengths[10] = sqrt((coord[14][0]-coord[10][0])*(coord[14][0]-coord[10][0]) + 
			      (coord[14][1]-coord[10][1])*(coord[14][1]-coord[10][1]));
	  elengths[13] = sqrt((coord[14][0]-coord[13][0])*(coord[14][0]-coord[13][0]) + 
			      (coord[14][1]-coord[13][1])*(coord[14][1]-coord[13][1]));
	  elengths[14] = sqrt((coord[14][0]-coord[15][0])*(coord[14][0]-coord[15][0]) + 
			      (coord[14][1]-coord[15][1])*(coord[14][1]-coord[15][1]));
	  
	  elengths[11] = sqrt((coord[16][0]-coord[11][0])*(coord[16][0]-coord[11][0]) + 
			      (coord[16][1]-coord[11][1])*(coord[16][1]-coord[11][1]));
	  elengths[15] = sqrt((coord[16][0]-coord[15][0])*(coord[16][0]-coord[15][0]) + 
			      (coord[16][1]-coord[15][1])*(coord[16][1]-coord[15][1]));
	  elengths[16] = sqrt((coord[16][0]-coord[17][0])*(coord[16][0]-coord[17][0]) + 
			      (coord[16][1]-coord[17][1])*(coord[16][1]-coord[17][1]));
	  
	  double x, y;
	  x = coord[9][0]-coord[0][0];
	  y = coord[9][1]-coord[0][1];
	  normals[0][0] = y/elengths[0];
	  normals[0][1] = -x/elengths[0];
	  if (normals[0][0]*(v1[0]-v0[0])+normals[0][1]*(v1[1]-v0[1]) < 0)
	    {
	      normals[0][0] = -normals[0][0];
	      normals[0][1] = -normals[0][1];
	    }
	  normals[3][0] = normals[0][0];
	  normals[3][1] = normals[0][1];
	  normals[6][0] = normals[0][0];
	  normals[6][1] = normals[0][1];
	  
	  x = coord[9][0]-coord[8][0];
	  y = coord[9][1]-coord[8][1];
	  normals[1][0] = y/elengths[1];
	  normals[1][1] = -x/elengths[1];
	  if (normals[1][0]*(v0[0]-v2[0])+normals[1][1]*(v0[1]-v2[1]) < 0)
	    {
	      normals[1][0] = -normals[1][0];
	      normals[1][1] = -normals[1][1];
	    }
	  normals[4][0] = normals[1][0];
	  normals[4][1] = normals[1][1];
	  normals[8][0] = normals[1][0];
	  normals[8][1] = normals[1][1];
	  
	  x = coord[10][0]-coord[3][0];
	  y = coord[10][1]-coord[3][1];
	  normals[2][0] = y/elengths[2];
	  normals[2][1] = -x/elengths[2];
	  if (normals[2][0]*(v2[0]-v1[0])+normals[2][1]*(v2[1]-v1[1]) < 0)
	    {
	      normals[2][0] = -normals[2][0];
	      normals[2][1] = -normals[2][1];
	    }
	  normals[5][0] = normals[2][0];
	  normals[5][1] = normals[2][1];
	  normals[7][0] = normals[2][0];
	  normals[7][1] = normals[2][1];
	  
	  x = coord[12][0]-coord[9][0];
	  y = coord[12][1]-coord[9][1];
	  normals[9][0] = y/elengths[9];
	  normals[9][1] = -x/elengths[9];
	  if (normals[9][0]*(v2[0]-v1[0])+normals[9][1]*(v2[1]-v1[1]) < 0)
	    {
	      normals[9][0] = -normals[9][0];
	      normals[9][1] = -normals[9][1];
	    }
	  normals[13][0] = normals[9][0];
	  normals[13][1] = normals[9][1];
	  normals[16][0] = normals[9][0];
	  normals[16][1] = normals[9][1];

	  x = coord[14][0]-coord[10][0];
	  y = coord[14][1]-coord[10][1];
	  normals[10][0] = y/elengths[10];
	  normals[10][1] = -x/elengths[10];
	  if (normals[10][0]*(v0[0]-v2[0])+normals[10][1]*(v0[1]-v2[1]) < 0)
	    {
	      normals[10][0] = -normals[10][0];
	      normals[10][1] = -normals[10][1];
	    }
	  normals[12][0] = normals[10][0];
	  normals[12][1] = normals[10][1];
	  normals[15][0] = normals[10][0];
	  normals[15][1] = normals[10][1];

	  x = coord[16][0]-coord[11][0];
	  y = coord[16][1]-coord[11][1];
	  normals[11][0] = y/elengths[11];
	  normals[11][1] = -x/elengths[11];
	  if (normals[11][0]*(v1[0]-v0[0])+normals[11][1]*(v1[1]-v0[1]) < 0)
	    {
	      normals[11][0] = -normals[11][0];
	      normals[11][1] = -normals[11][1];
	    }
	  normals[14][0] = normals[11][0];
	  normals[14][1] = normals[11][1];
	  normals[17][0] = normals[11][0];
	  normals[17][1] = normals[11][1];
      
	  double *p, *q;
	  p = new double[2];
	  q = new double[2];
	  p[0] = coord[8][0]-coord[0][0];
	  p[1] = coord[8][1]-coord[0][1];
	  q[0] = coord[9][0]-v0[0];
	  q[1] = coord[9][1]-v0[1];
	  areas[0] = 0.5*fabs(p[0]*q[1]-p[1]*q[0]);
	  areas[1] = areas[0];
	  areas[2] = areas[0];
	  areas[3] = 3 * areas[0];
	  areas[4] = 3 * areas[0];
	  areas[5] = 3 * areas[0];
	  areas[6] = 3 * areas[0];
	  areas[7] = 3 * areas[0];
	  areas[8] = 3 * areas[0];
	  areas[9] = 6 * areas[0];	  
	  delete []p;
	  delete []q;	  
	}
      else
	{
	  cout << "DualMesh is under construction for order: " << dualmeshorder << endl;
	}
    }
  else
    {
      double *v0 = v[0];
      double *v1 = v[1];
      double *v2 = v[2];
      double *v3 = v[3];
      if(dualmeshorder==1)
	{
	  for(int k=0; k<2; k++)
	    {
	      coord[0][k] = (v1[k]+v0[k]) / 2.0;
	      coord[1][k] = (v2[k]+v1[k]) / 2.0;
	      coord[2][k] = (v2[k]+v3[k]) / 2.0;
	      coord[3][k] = (v3[k]+v0[k]) / 2.0;
	      coord[4][k] = (v0[k]+v1[k]+v2[k]+v3[k]) / 4.0;
	    }

	  normals[0][0] = 1.0;
	  normals[0][1] = 0.0;
	  normals[1][0] = 0.0;
	  normals[1][1] = 1.0;
	  normals[2][0] = -1.0;
	  normals[2][1] = 0.0;
	  normals[3][0] = 0.0;
	  normals[3][1] = -1.0;
	  
	  for(int k=0; k<4; k++)
	    elengths[k] = sqrt((coord[k][0]-coord[4][0])*(coord[k][0]-coord[4][0]) + 
			       (coord[k][1]-coord[4][1])*(coord[k][1]-coord[4][1]));

	  areas[0] = fabs( (coord[4][1]-coord[0][1])*(coord[4][0]-coord[3][0]) );
	  areas[1] = areas[0];
	  areas[2] = areas[0];
	  areas[3] = areas[0];
	}
      else
	cout << "Under construction..." << endl;
    }  
}

//=============================================================================

void DualMesh::GetElemDualCoord(int i, Array<double *> &coord)
{
  int n = (mesh->GetMType()==Element::TRIANGLE) ? 3 : 4;  

  Array<int> ind(n);
  Array<double *> v(n);
  mesh->GetElementVertices (i, ind);
  for (int k=0; k<n; k++) { v[k] = mesh->GetVertex(ind[k]); }

  if(n==3)
    {
      double *v0 = v[0];
      double *v1 = v[1];
      double *v2 = v[2];
      if(dualmeshorder==1)
	{
	  for(int k=0; k<2; k++)
	    {
	      coord[0][k] = (v1[k]+v0[k]) / 2.0;
	      coord[1][k] = (v2[k]+v1[k]) / 2.0;
	      coord[2][k] = (v2[k]+v0[k]) / 2.0;
	      coord[3][k] = (v0[k]+v1[k]+v2[k]) / 3.0;
	    }
	}
      else if(dualmeshorder==2)
	{
	  if(dualmeshtype==0)
	    {
	      for(int k=0; k<2; k++)
		{
		  coord[0][k] = v0[k]*3.0/4.0 + v1[k]/4.0; 
		  coord[1][k] = v0[k]/2.0 + v1[k]/4.0 + v2[k]/4.0; 
		  coord[2][k] = v0[k]*3.0/4.0 + v2[k]/4.0;
		  coord[3][k] = v0[k]*2.0/3.0 + v1[k]/6.0 + v2[k]/6.0;
		  
		  coord[4][k] = v1[k]*3.0/4.0 + v2[k]/4.0; 
		  coord[5][k] = v1[k]/2.0 + v0[k]/4.0 + v2[k]/4.0; 
		  coord[6][k] = v1[k]*3.0/4.0 + v0[k]/4.0;
		  coord[7][k] = v1[k]*2.0/3.0 + v0[k]/6.0 + v2[k]/6.0;
		  
		  coord[8][k] = v2[k]*3.0/4.0 + v0[k]/4.0; 
		  coord[9][k] = v2[k]/2.0 + v1[k]/4.0 + v0[k]/4.0; 
		  coord[10][k] = v2[k]*3.0/4.0 + v1[k]/4.0;
		  coord[11][k] = v2[k]*2.0/3.0 + v1[k]/6.0 + v0[k]/6.0;
	  	  
		  coord[12][k] = (v0[k] + v1[k] + v2[k])/3.0;
		}      
	    }
	  else if(dualmeshtype==1)
	    {
	      for(int k=0; k<2; k++)
		{
		  coord[0][k] = v0[k]*5.0/6.0 + v1[k]/6.0; 
		  coord[1][k] = v0[k]/2.0 + v1[k]/4.0 + v2[k]/4.0; 
		  coord[2][k] = v0[k]*5.0/6.0 + v2[k]/6.0;
		  coord[3][k] = v0[k]*2.0/3.0 + v1[k]/6.0 + v2[k]/6.0;
		  
		  coord[4][k] = v1[k]*5.0/6.0 + v2[k]/6.0; 
		  coord[5][k] = v1[k]/2.0 + v0[k]/4.0 + v2[k]/4.0; 
		  coord[6][k] = v1[k]*5.0/6.0 + v0[k]/6.0;
		  coord[7][k] = v1[k]*2.0/3.0 + v0[k]/6.0 + v2[k]/6.0;
		  
		  coord[8][k] = v2[k]*5.0/6.0 + v0[k]/6.0; 
		  coord[9][k] = v2[k]/2.0 + v1[k]/4.0 + v0[k]/4.0; 
		  coord[10][k] = v2[k]*5.0/6.0 + v1[k]/6.0;
		  coord[11][k] = v2[k]*2.0/3.0 + v1[k]/6.0 + v0[k]/6.0;
	  	  
		  coord[12][k] = (v0[k] + v1[k] + v2[k])/3.0;
		}      
	    }
	  else
	    cout << "Quadratic DualMesh is not constructed for this type: " << dualmeshtype << endl;	  
	}
      else if(dualmeshorder==3)
	{
	  for(int k=0; k<2; k++)
	    {
	      coord[0][k] = v0[k]*5.0/6.0 + v1[k]/6.0;   
	      coord[1][k] = v0[k]/2.0 + v1[k]/2.0;
	      coord[2][k] = v1[k]*5.0/6.0 + v0[k]/6.0;
	      
	      coord[3][k] = v1[k]*5.0/6.0 + v2[k]/6.0;   
	      coord[4][k] = v1[k]/2.0 + v2[k]/2.0;
	      coord[5][k] = v2[k]*5.0/6.0 + v1[k]/6.0;
	      
	      coord[6][k] = v2[k]*5.0/6.0 + v0[k]/6.0;   
	      coord[7][k] = v0[k]/2.0 + v2[k]/2.0;
	      coord[8][k] = v0[k]*5.0/6.0 + v2[k]/6.0;
	      
	      coord[9][k] = v0[k]*7.0/9.0 + v1[k]/9.0 + v2[k]/9.0;
	      coord[10][k] = v1[k]*7.0/9.0 + v0[k]/9.0 + v2[k]/9.0;
	      coord[11][k] = v2[k]*7.0/9.0 + v1[k]/9.0 + v0[k]/9.0;
	      
	      coord[12][k] = v0[k]*5.0/9.0 + v1[k]*2.0/9.0 + v2[k]*2.0/9.0;
	      coord[14][k] = v1[k]*5.0/9.0 + v0[k]*2.0/9.0 + v2[k]*2.0/9.0;
	      coord[16][k] = v2[k]*5.0/9.0 + v1[k]*2.0/9.0 + v0[k]*2.0/9.0;
	      
	      coord[13][k] = v0[k]*4.0/9.0 + v1[k]*4.0/9.0 + v2[k]/9.0;
	      coord[15][k] = v2[k]*4.0/9.0 + v1[k]*4.0/9.0 + v0[k]/9.0;
	      coord[17][k] = v0[k]*4.0/9.0 + v2[k]*4.0/9.0 + v1[k]/9.0;
	      
	      coord[18][k] = v0[k]/3.0 + v1[k]/3.0 + v2[k]/3.0;
	    }
	}
      else
	{
	  cout << "DualMesh is under construction for order: " << dualmeshorder << endl;
	}
    }
  else
    {
      cout << "Under construction..." << endl;
    }
}

//=============================================================================

void DualMesh::GetElemDualNormals(int i, Array<double *> &normals)
{
  ;
}

//=============================================================================

void DualMesh::GetElemDualEdgeLengths(int i, Array<double> &elengths)
{
  ;
}

//=============================================================================

void DualMesh::GetElemDualAreas(int i, Array<double> &areas)
{
  int n = (mesh->GetMType()==Element::TRIANGLE) ? 3 : 4;  

  Array<int> ind(n);
  Array<double *> v(n);
  mesh->GetElementVertices (i, ind);
  for (int k=0; k<n; k++) { v[k] = mesh->GetVertex(ind[k]); }

  if(n==3)
    {
      double *v0 = v[0];
      double *v1 = v[1];
      double *v2 = v[2];
      if(dualmeshorder==1)
	{
	  double coord[4][2];
	  for(int k=0; k<2; k++)
	    {
	      coord[0][k] = (v1[k]+v0[k]) / 2.0;
	      coord[1][k] = (v2[k]+v1[k]) / 2.0;
	      coord[2][k] = (v2[k]+v0[k]) / 2.0;
	      coord[3][k] = (v0[k]+v1[k]+v2[k]) / 3.0;
	    }

	  double *p, *q;
	  p = new double[2];
	  q = new double[2];
	  p[0] = coord[2][0]-coord[0][0];
	  p[1] = coord[2][1]-coord[0][1];
	  q[0] = coord[3][0]-v0[0];
	  q[1] = coord[3][1]-v0[1];
	  areas[0] = 0.5*fabs(p[0]*q[1]-p[1]*q[0]);
	  areas[1] = areas[0]; 
	  areas[2] = areas[1];
	  delete []p;
	  delete []q;
	}
      else if(dualmeshorder==2)
	{
	  if(dualmeshtype==0)
	    {
	      double coord[4][2];
	      for(int k=0; k<2; k++)
		{
		  coord[0][k] = v0[k]*3.0/4.0 + v1[k]/4.0; 
		  coord[1][k] = v0[k]/2.0 + v1[k]/4.0 + v2[k]/4.0; 
		  coord[2][k] = v0[k]*3.0/4.0 + v2[k]/4.0;
		  coord[3][k] = v0[k]*2.0/3.0 + v1[k]/6.0 + v2[k]/6.0;
		}	      
	      double *p, *q;
	      p = new double[2];
	      q = new double[2];
	      p[0] = coord[2][0]-coord[0][0];
	      p[1] = coord[2][1]-coord[0][1];
	      q[0] = coord[3][0]-v0[0];
	      q[1] = coord[3][1]-v0[1];
	      areas[0] = 0.5*fabs(p[0]*q[1]-p[1]*q[0]);
	      delete []p;
	      delete []q;
	    }
	  
	  else if(dualmeshtype==1)
	    {
	      double coord[4][2];
	      for(int k=0; k<2; k++)
		{
		  coord[0][k] = v0[k]*5.0/6.0 + v1[k]/6.0; 
		  coord[1][k] = v0[k]/2.0 + v1[k]/4.0 + v2[k]/4.0; 
		  coord[2][k] = v0[k]*5.0/6.0 + v2[k]/6.0;
		  coord[3][k] = v0[k]*2.0/3.0 + v1[k]/6.0 + v2[k]/6.0;
		}
	      
	      double *p, *q;
	      p = new double[2];
	      q = new double[2];
	      p[0] = coord[2][0]-coord[0][0];
	      p[1] = coord[2][1]-coord[0][1];
	      q[0] = coord[3][0]-v0[0];
	      q[1] = coord[3][1]-v0[1];
	      areas[0] = 0.5*fabs(p[0]*q[1]-p[1]*q[0]);

	      delete []p;
	      delete []q;
	    }
	  else
	    cout << "Quadratic DualMesh is not constructed for this type: " << dualmeshtype << endl;	  
	}
      else if(dualmeshorder==3)
	{
	  double coord[4][2];
	  for(int k=0; k<2; k++)
	    {
	      coord[0][k] = v0[k]*5.0/6.0 + v1[k]/6.0;   
	      coord[1][k] = v0[k]*5.0/6.0 + v2[k]/6.0;	      
	      coord[2][k] = v0[k]*7.0/9.0 + v1[k]/9.0 + v2[k]/9.0;
	    }

	  double *p, *q;
	  p = new double[2];
	  q = new double[2];
	  p[0] = coord[1][0]-coord[0][0];
	  p[1] = coord[1][1]-coord[0][1];
	  q[0] = coord[2][0]-v0[0];
	  q[1] = coord[2][1]-v0[1];
	  areas[0] = 0.5*fabs(p[0]*q[1]-p[1]*q[0]);
	  areas[1] = areas[0];
	  areas[2] = areas[0];
	  areas[3] = 3 * areas[0];
	  areas[4] = 3 * areas[0];
	  areas[5] = 3 * areas[0];
	  areas[6] = 3 * areas[0];
	  areas[7] = 3 * areas[0];
	  areas[8] = 3 * areas[0];
	  areas[9] = 6 * areas[0];	  
	  delete []p;
	  delete []q;	  
	}
      else
	{
	  cout << "DualMesh is under construction for order: " << dualmeshorder << endl;
	}
    }
  else
    {
      double *v0 = v[0];
      //double *v1 = v[1];
      double *v2 = v[2];
      double *v3 = v[3];
      if(dualmeshorder==1)
	{
	  areas[0] = 0.25 * fabs( (v3[1]-v0[1])*(v2[0]-v0[0]) );
	  areas[1] = areas[0];
	  areas[2] = areas[0];
	  areas[3] = areas[0];
	}
      else
	cout << "Under construction..." << endl;
    }  
}

/*/=============================================================================

void DualMesh::PrintDualMesh()
{
  int pcount = 0;
  for(int k=0; k<dof; k++)
    {
      if(cvnumel[k]==1){ pcount += 4;}
      else{ pcount += 3*cvnumel[k];}
    }


  ofstream out;
  out.open("dualmesh.vtk");
  out << "# vtk DataFile Version 2.0" << endl;
  out << "Unstructured Grid" << endl;
  out << "ASCII" << endl;
  out << "DATASET UNSTRUCTURED_GRID" << endl;
  out << "POINTS " << pcount << " double" << endl;

  for(int k=0; k<dof; k++)
    {
      if(cvnumel[k]==1)
	{
	  int el = cvelindex[k][0];
	  int locvert = cvlocalindex[k][0];

	  Array<double *> coord;
	  cvolumes[el]->GetCoordinates(coord);
	  Array<double *> vcoord(2);
	  for(int i=0; i<2; i++){vcoord[i] = new double[numlocaldof];}
	  mesh->GetElementVerticesCoord(el, vcoord);

	  //cout <<  el << " " << locvert << " " << vcoord[0][locvert] << " " << vcoord[1][locvert] << endl;
	  for(int i=0; i<GetDim(); i++)
	    {
	      out << coord[locvert][i] << " ";
	    }
	  out << 0.0 << endl;

	  for(int i=0; i<GetDim(); i++)
	    {
	      out << coord[numlocaldof][i] << " ";
	    }
	  out << 0.0 << endl;

	  for(int i=0; i<GetDim(); i++)
	    {
	      out << coord[(locvert+numlocaldof-1)%numlocaldof][i] << " ";
	    }
	  out << 0.0 << endl;

	  for(int i=0; i<GetDim(); i++)
	    {
	      out << vcoord[i][locvert] << " ";
	    }
	  out << 0.0 << endl;

	  for(int i=0; i<2; i++){delete []vcoord[i]; }
	}
      else
	{
	  Array<int> tempelindex(cvnumel[k]);
	  Array<int> templocindexstart(cvnumel[k]);
	  Array<int> templocindexend(cvnumel[k]);
	  int locvert = cvlocalindex[k][0];
	  tempelindex[0] = cvelindex[k][0];
	  templocindexstart[0] = locvert;
	  templocindexend[0] = (locvert+numlocaldof-1)%numlocaldof;

	  Array<double *> coord;
	  cvolumes[tempelindex[0]]->GetCoordinates(coord);
	  double xtest = coord[templocindexend[0]][0];
	  double ytest = coord[templocindexend[0]][1];

	  tempelindex[cvnumel[k]-1] = -1;
	  tempelindex[1] = -1;
	  //cout << k << endl;
	  int indexcount = 1;
	  while(tempelindex[cvnumel[k]-1] == -1)
	    {
	      for(int p=1; p<cvnumel[k]; p++)
		{
		  int el = cvelindex[k][p];
		  int locvertex = cvlocalindex[k][p];
		  int locbdrind = (locvertex+numlocaldof-1)%numlocaldof; 
		  cvolumes[el]->GetCoordinates(coord);

		  if( fabs(xtest - coord[locvertex][0]) < 1e-7 && fabs(ytest - coord[locvertex][1]) < 1e-7 && el!=tempelindex[indexcount-1])
		    {
		      //cout << 1 << " " << el << " " <<  locvertex << " " << locbdrind << " " << indexcount << endl;
		      tempelindex[indexcount] = cvelindex[k][p];
		      templocindexstart[indexcount] = locvertex;
		      templocindexend[indexcount] = locbdrind;

		      cvolumes[tempelindex[indexcount]]->GetCoordinates(coord);
		      xtest = coord[templocindexend[indexcount]][0];
		      ytest = coord[templocindexend[indexcount]][1];
		      indexcount++;
		      break;
		    }
	      
		 
		  if( fabs(xtest - coord[locbdrind][0]) < 1e-7 && fabs(ytest - coord[locbdrind][1]) < 1e-7 && el!=tempelindex[indexcount-1])
		    {
		      //cout << 2 << " " << el << " " <<  locvertex << " " << locbdrind << " " << indexcount << endl;
		      tempelindex[indexcount] = cvelindex[k][p];
		      templocindexstart[indexcount] = locbdrind;
		      templocindexend[indexcount] = locvertex;
		  
		      cvolumes[tempelindex[indexcount]]->GetCoordinates(coord);
		      xtest = coord[templocindexend[indexcount]][0];
		      ytest = coord[templocindexend[indexcount]][1];
		      indexcount++;
		      break;
		    }
		}

	      if(tempelindex[1] == -1)
		{
		  int temp = templocindexstart[0];
		  templocindexstart[0] = templocindexend[0];
		  templocindexend[0] = temp;

		  // cout << templocindexend[0] << " " << templocindexstart[0] << endl;
		  cvolumes[tempelindex[0]]->GetCoordinates(coord);
		  xtest = coord[templocindexend[0]][0];
		  ytest = coord[templocindexend[0]][1];
		}	      
	    } 

	  for(int p=0; p<cvnumel[k]; p++)
	    {  
	      Array<double *> coord;
	      cvolumes[tempelindex[p]]->GetCoordinates(coord);
	      
	      for(int i=0; i<GetDim(); i++)
		{
		  out << coord[templocindexstart[p]][i] << " ";
		}
	      out << 0.0 << endl;
	      
	      for(int i=0; i<GetDim(); i++)
		{
		  out << coord[numlocaldof][i] << " ";
		}
	      out << 0.0 << endl;
	      
	      for(int i=0; i<GetDim(); i++)
		{
		  out << coord[templocindexend[p]][i] << " ";
		}
	      out << 0.0 << endl;
	    }
	}
    }
 
  out << endl;
  out << "CELLS " << dof << " " << dof + pcount << endl;

  int counter = 0;
  for(int k=0; k<dof; k++)
    {
      if(cvnumel[k]==1)
	{
	  out << 4 << "\t";
	  for(int i=0; i<4; i++)
	    {
	      out << counter << "\t";
	      counter++;
	    }
	  out << endl;
	}
      else
	{
	  out << 3*cvnumel[k] << "\t";
	  for(int i=0; i<3*cvnumel[k]; i++)
	    {
	      out << counter << "\t";
	      counter++;
	    }
	  out << endl;
	}
    }

  out << endl;


  out << "CELL_TYPES " << dof << endl;
  for(int j=0; j<dof; j++)
    {
      out << 7 << endl;
    }
  out.close();
}

//=============================================================================

void DualMesh::PrintSaturation(ofstream &out, Vector &Saturation, double scalex, double scaley)
{
  int pcount = 0;
  for(int k=0; k<dof; k++)
    {
      if(cvnumel[k]==1){ pcount += 4;}
      else{ pcount += 3*cvnumel[k];}
    }

  out << "# vtk DataFile Version 2.0" << endl;
  out << "Unstructured Grid" << endl;
  out << "ASCII" << endl;
  out << "DATASET UNSTRUCTURED_GRID" << endl;
  out << "POINTS " << pcount << " double" << endl;

  for(int k=0; k<dof; k++)
    {
      if(cvnumel[k]==1)
	{
	  int el = cvelindex[k][0];
	  int locvert = cvlocalindex[k][0];

	  Array<double *> coord;
	  cvolumes[el]->GetCoordinates(coord);

	  Array<double *> vcoord(2);
	  for(int i=0; i<2; i++){vcoord[i] = new double[numlocaldof];}
	  mesh->GetElementVerticesCoord(el, vcoord);

	  //cout <<  el << " " << locvert << " " << vcoord[0][locvert] << " " << vcoord[1][locvert] << endl;
	  out << scalex*coord[locvert][0] << " " << scaley*coord[locvert][1] << " " << 0.0 << endl;
	  out << scalex*coord[numlocaldof][0] << " " << scaley*coord[numlocaldof][1] << " " << 0.0 << endl;
	  out << scalex*coord[(locvert+numlocaldof-1)%numlocaldof][0] << " " << scaley*coord[(locvert+numlocaldof-1)%numlocaldof][1] << " " << 0.0 << endl;
	  out << scalex*vcoord[0][locvert] << " " << scaley*vcoord[1][locvert] << " " << 0.0 << endl;

	  for(int i=0; i<2; i++){delete []vcoord[i]; }
	}
      else
	{
	  Array<int> tempelindex(cvnumel[k]);
	  Array<int> templocindexstart(cvnumel[k]);
	  Array<int> templocindexend(cvnumel[k]);
	  int locvert = cvlocalindex[k][0];
	  tempelindex[0] = cvelindex[k][0];
	  templocindexstart[0] = locvert;
	  templocindexend[0] = (locvert+numlocaldof-1)%numlocaldof;

	  Array<double *> coord;
	  cvolumes[tempelindex[0]]->GetCoordinates(coord);
	  double xtest = coord[templocindexend[0]][0];
	  double ytest = coord[templocindexend[0]][1];

	  tempelindex[cvnumel[k]-1] = -1;
	  tempelindex[1] = -1;
	  //cout << k << endl;
	  int indexcount = 1;
	  while(tempelindex[cvnumel[k]-1] == -1)
	    {
	      for(int p=1; p<cvnumel[k]; p++)
		{
		  int el = cvelindex[k][p];
		  int locvertex = cvlocalindex[k][p];
		  int locbdrind = (locvertex+numlocaldof-1)%numlocaldof; 
		  cvolumes[el]->GetCoordinates(coord);

		  if( fabs(xtest - coord[locvertex][0]) < 1e-7 && fabs(ytest - coord[locvertex][1]) < 1e-7 && el!=tempelindex[indexcount-1])
		    {
		      //cout << 1 << " " << el << " " <<  locvertex << " " << locbdrind << " " << indexcount << endl;
		      tempelindex[indexcount] = cvelindex[k][p];
		      templocindexstart[indexcount] = locvertex;
		      templocindexend[indexcount] = locbdrind;

		      cvolumes[tempelindex[indexcount]]->GetCoordinates(coord);
		      xtest = coord[templocindexend[indexcount]][0];
		      ytest = coord[templocindexend[indexcount]][1];
		      indexcount++;
		      break;
		    }
	      
		 
		  if( fabs(xtest - coord[locbdrind][0]) < 1e-7 && fabs(ytest - coord[locbdrind][1]) < 1e-7 && el!=tempelindex[indexcount-1])
		    {
		      //cout << 2 << " " << el << " " <<  locvertex << " " << locbdrind << " " << indexcount << endl;
		      tempelindex[indexcount] = cvelindex[k][p];
		      templocindexstart[indexcount] = locbdrind;
		      templocindexend[indexcount] = locvertex;
		  
		      cvolumes[tempelindex[indexcount]]->GetCoordinates(coord);
		      xtest = coord[templocindexend[indexcount]][0];
		      ytest = coord[templocindexend[indexcount]][1];
		      indexcount++;
		      break;
		    }
		}

	      if(tempelindex[1] == -1)
		{
		  int temp = templocindexstart[0];
		  templocindexstart[0] = templocindexend[0];
		  templocindexend[0] = temp;

		  // cout << templocindexend[0] << " " << templocindexstart[0] << endl;
		  cvolumes[tempelindex[0]]->GetCoordinates(coord);
		  xtest = coord[templocindexend[0]][0];
		  ytest = coord[templocindexend[0]][1];
		}	      
	    } 

	  for(int p=0; p<cvnumel[k]; p++)
	    {  
	      Array<double *> coord;
	      cvolumes[tempelindex[p]]->GetCoordinates(coord);
	      
	      out << scalex*coord[templocindexstart[p]][0] << " " << scaley*coord[templocindexstart[p]][1] << " " << 0.0 << endl;
	      out << scalex*coord[numlocaldof][0] << " " << scaley*coord[numlocaldof][1] << " " << 0.0 << endl;
	      out << scalex*coord[templocindexend[p]][0] << " " << scaley*coord[templocindexend[p]][1] << " " << 0.0 << endl;
	  
	    }
	}
    }
 
  out << endl;
  out << "CELLS " << dof << " " << dof + pcount << endl;

  int counter = 0;
  for(int k=0; k<dof; k++)
    {
      if(cvnumel[k]==1)
	{
	  out << 4 << "\t";
	  for(int i=0; i<4; i++)
	    {
	      out << counter << "\t";
	      counter++;
	    }
	  out << endl;
	}
      else
	{
	  out << 3*cvnumel[k] << "\t";
	  for(int i=0; i<3*cvnumel[k]; i++)
	    {
	      out << counter << "\t";
	      counter++;
	    }
	  out << endl;
	}
    }

  out << endl;


  out << "CELL_TYPES " << dof << endl;
  for(int j=0; j<dof; j++)
    {
      out << 7 << endl;
    }

  out << "CELL_DATA " << dof << endl;
  out << "SCALARS pval DOUBLE" << endl;
  out << "LOOKUP_TABLE default" << endl;

  for(int j=0; j<dof; j++)
    {  
      out << Saturation(j) << endl;
    }
  
}

//=============================================================================

void DualMesh::PrintGLESaturation(ofstream &out, Vector &Saturation, double scalex, double scaley)
{

  out << dof << endl;
  for(int k=0; k<dof; k++)
    {
      
      if(cvnumel[k]==1)
	{
	  out << 8 << " " << Saturation(k) << endl;
	  int el = cvelindex[k][0];
	  int locvert = cvlocalindex[k][0];
	  Array<double *> coord;
	  cvolumes[el]->GetCoordinates(coord);

	  Array<double *> vcoord(2);
	  for(int i=0; i<2; i++){vcoord[i] = new double[numlocaldof];}
	  mesh->GetElementVerticesCoord(el, vcoord);

	  out << scalex*coord[locvert][0] << " " << scaley*coord[locvert][1] << " ";
	  out << scalex*coord[numlocaldof][0] << " " << scaley*coord[numlocaldof][1] << " ";
	  out << scalex*coord[(locvert+numlocaldof-1)%numlocaldof][0] << " " << scaley*coord[(locvert+numlocaldof-1)%numlocaldof][1] << " ";
	  out << scalex*vcoord[0][locvert] << " " << scaley*vcoord[1][locvert] << " ";
	  out << endl;
	  for(int i=0; i<2; i++){delete []vcoord[i]; }
	}
      else
	{
	  out << 6*cvnumel[k] << " " << Saturation(k) << endl;
	  Array<int> tempelindex(cvnumel[k]);
	  Array<int> templocindexstart(cvnumel[k]);
	  Array<int> templocindexend(cvnumel[k]);
	  int locvert = cvlocalindex[k][0];
	  tempelindex[0] = cvelindex[k][0];
	  templocindexstart[0] = locvert;
	  templocindexend[0] = (locvert+numlocaldof-1)%numlocaldof;

	  Array<double *> coord;
	  cvolumes[tempelindex[0]]->GetCoordinates(coord);
	  double xtest = coord[templocindexend[0]][0];
	  double ytest = coord[templocindexend[0]][1];

	  tempelindex[cvnumel[k]-1] = -1;
	  tempelindex[1] = -1;
	  //cout << k << endl;
	  int indexcount = 1;
	  while(tempelindex[cvnumel[k]-1] == -1)
	    {
	      for(int p=1; p<cvnumel[k]; p++)
		{
		  int el = cvelindex[k][p];
		  int locvertex = cvlocalindex[k][p];
		  int locbdrind = (locvertex+numlocaldof-1)%numlocaldof; 
		  cvolumes[el]->GetCoordinates(coord);

		  if( fabs(xtest - coord[locvertex][0]) < 1e-7 && fabs(ytest - coord[locvertex][1]) < 1e-7 && el!=tempelindex[indexcount-1])
		    {
		      //cout << 1 << " " << el << " " <<  locvertex << " " << locbdrind << " " << indexcount << endl;
		      tempelindex[indexcount] = cvelindex[k][p];
		      templocindexstart[indexcount] = locvertex;
		      templocindexend[indexcount] = locbdrind;

		      cvolumes[tempelindex[indexcount]]->GetCoordinates(coord);
		      xtest = coord[templocindexend[indexcount]][0];
		      ytest = coord[templocindexend[indexcount]][1];
		      indexcount++;
		      break;
		    }
	      
		 
		  if( fabs(xtest - coord[locbdrind][0]) < 1e-7 && fabs(ytest - coord[locbdrind][1]) < 1e-7 && el!=tempelindex[indexcount-1])
		    {
		      //cout << 2 << " " << el << " " <<  locvertex << " " << locbdrind << " " << indexcount << endl;
		      tempelindex[indexcount] = cvelindex[k][p];
		      templocindexstart[indexcount] = locbdrind;
		      templocindexend[indexcount] = locvertex;
		  
		      cvolumes[tempelindex[indexcount]]->GetCoordinates(coord);
		      xtest = coord[templocindexend[indexcount]][0];
		      ytest = coord[templocindexend[indexcount]][1];
		      indexcount++;
		      break;
		    }
		}

	      if(tempelindex[1] == -1)
		{
		  int temp = templocindexstart[0];
		  templocindexstart[0] = templocindexend[0];
		  templocindexend[0] = temp;

		  // cout << templocindexend[0] << " " << templocindexstart[0] << endl;
		  cvolumes[tempelindex[0]]->GetCoordinates(coord);
		  xtest = coord[templocindexend[0]][0];
		  ytest = coord[templocindexend[0]][1];
		}	      
	    } 

	  for(int p=0; p<cvnumel[k]; p++)
	    {  
	      Array<double *> coord;
	      cvolumes[tempelindex[p]]->GetCoordinates(coord);
	      
	      out << scalex*coord[templocindexstart[p]][0] << " " << scaley*coord[templocindexstart[p]][1] << " ";
	      out << scalex*coord[numlocaldof][0] << " " << scaley*coord[numlocaldof][1] << " ";
	      out << scalex*coord[templocindexend[p]][0] << " " << scaley*coord[templocindexend[p]][1] << " ";	  
	    }
	  out << endl;
	}
    }
  
}
*/

//=============================================================================

DualMesh::~DualMesh()
{
  ;
}
