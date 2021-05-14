#include "mesh_header.h"


//=============================================================================

ControlVolume::ControlVolume(const Array<double *> &v, const Array<int> &dof)
{
  numdof = dof.Size();  
  vcoord.SetSize(numdof);
  for (int i=0; i<numdof; i++) 
    { 
      vcoord[i] = new double[2];
      for(int j=0; j<2; j++){ vcoord[i][j] = v[i][j]; }
    }

  vindex.SetSize(dof.Size());
  for (int i=0; i<dof.Size(); i++)
    {
      vindex[i] = dof[i]; 
    }

  coord.SetSize(numdof+1);
  for (int i=0; i<=numdof; i++) 
    { 
      coord[i] = new double[2];
      for(int j=0; j<2; j++){ coord[i][j] = 0.0; }
    }

  for(int j=0; j<2; j++)
    {
      for(int i=0; i<numdof; i++)
	{
	  coord[i][j] = 0.5*(v[(i+1)%numdof][j] + v[i][j]);
	  coord[numdof][j] += v[i][j];
	}
      coord[numdof][j] /= (double) numdof;
    }
  
  elengths.SetSize(numdof);
  for(int i=0; i<numdof; i++)
    {
      double diffx = coord[i][0]-coord[numdof][0];
      double diffy = coord[i][1]-coord[numdof][1];
      elengths[i] = sqrt(diffx*diffx+diffy*diffy);
    }
 
 
  normals.SetSize(numdof);
  for (int i=0; i<numdof; i++) { normals[i] = new double[2]; }
  for(int i=0; i<numdof; i++)
    {
      double x = coord[numdof][0]-coord[i][0];
      double y = coord[numdof][1]-coord[i][1];
      
      normals[i][0] = -y/elengths[i];
      normals[i][1] = x/elengths[i];
      
      double vx = v[(i+1)%numdof][0]-v[i][0];
      double vy = v[(i+1)%numdof][1]-v[i][1];
      if (normals[i][0]*vx+normals[i][1]*vy < 0)
	{
	  normals[i][0] = -normals[i][0];
	  normals[i][1] = -normals[i][1];
	}   

    }
  
  
  areas.SetSize(numdof);
  for(int k=0; k<numdof; k++)
    {
      int kk = (k+numdof-1)%numdof;

      double diffx = (v[k][0]-coord[k][0]);
      double diffy = (v[k][1]-coord[k][1]);
      double a = sqrt(diffx*diffx+diffy*diffy);

      diffx = (coord[numdof][0]-coord[k][0]);
      diffy = (coord[numdof][1]-coord[k][1]);
      double b = sqrt(diffx*diffx+diffy*diffy);
      
      diffx = (coord[numdof][0]-coord[kk][0]);
      diffy = (coord[numdof][1]-coord[kk][1]);
      double c = sqrt(diffx*diffx+diffy*diffy);
      
      diffx = (v[k][0]-coord[kk][0]);
      diffy = (v[k][1]-coord[kk][1]);
      double d = sqrt(diffx*diffx+diffy*diffy);

      diffx = (v[k][0]-coord[numdof][0]);
      diffy = (v[k][1]-coord[numdof][1]);
      double q = sqrt(diffx*diffx+diffy*diffy);

      diffx = (coord[kk][0]-coord[k][0]);
      diffy = (coord[kk][1]-coord[k][1]);
      double p = sqrt(diffx*diffx+diffy*diffy);

      areas[k] = 0.25*sqrt(4.0*p*p*q*q-(b*b+d*d-a*a-c*c)*(b*b+d*d-a*a-c*c));
    }

    

  
}

//=============================================================================

void ControlVolume::GetEdgeLengths(Array<double> &_elengths)
{
  _elengths.SetSize(numdof);
  for(int i=0; i<numdof; i++)
    {
      _elengths[i] = elengths[i];
    }
}

//=============================================================================

void ControlVolume::GetElementVerticesCoord(Array<double *> &_vcoord)
{
  _vcoord.SetSize(numdof);
  for(int i=0; i<numdof; i++)
    {
      _vcoord[i] = vcoord[i];
    }

}

//=============================================================================

void ControlVolume::GetCoordinates(Array<double *> &_coord)
{
  _coord.SetSize(numdof+1);
  for(int i=0; i<=numdof; i++)
    {
      _coord[i] = coord[i];
    }
}

//=============================================================================

void ControlVolume::GetNormals(Array<double *> &_normals)
{
  _normals.SetSize(numdof); 
 for(int i=0; i<numdof; i++)
    {
      _normals[i] = normals[i];
    }
}

//=============================================================================

void ControlVolume::GetAreas(Array<double> &_areas)
{
  _areas.SetSize(numdof);
  for(int i=0; i<numdof; i++)
    {
      _areas[i] = areas[i];
    }
} 

//==============================================================================

void ControlVolume::GetDOF(Array<int> &_dof)
{
  _dof.SetSize(vindex.Size());
  for(int i=0; i<vindex.Size(); i++)
    {
      _dof[i] = vindex[i];
    }
} 

//==============================================================================

ControlVolume::~ControlVolume()
{
  for (int i=0; i<normals.Size(); i++) {delete []normals[i]; }
  for (int i=0; i<coord.Size(); i++) { delete []coord[i]; }
  for (int i=0; i<vcoord.Size(); i++) { delete []vcoord[i]; }

}
