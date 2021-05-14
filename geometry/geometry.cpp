#include "geometry_header.h"

Geometry::Geometry(ifstream *file)
{
  int temp;
  *file >> temp;
  NumOfPatches = temp;
  patches.SetSize(temp);
  faces.SetSize(temp);
  *file >> temp;
  PatchDim = temp;
  *file >> temp;
  SpaceDim = temp;
  for (int i=0; i<NumOfPatches; i++)
    {
      patches[i] = new Patch(file,PatchDim,SpaceDim);
    }
  /////
  int counter = faces.Size()*patches[0]->GetFaceSize();
  double *temp1 = new double[faces.Size()*patches[0]->GetFaceSize()];
  for (int i=0; i<faces.Size(); i++)
    {
      faces[i] = new Vector(patches[0]->GetFaceSize());
      for (int j=0; j<faces[0]->Size(); j++)
	{
	  temp1[i*faces[0]->Size()+j] = patches[i]->GetFaceVal(j);
	}
    }
  connect.SetSize(faces.Size()*patches[0]->GetFaceSize());
  for (int i=0; i<connect.Size(); i++)
    {
      connect[i] = new Vector(4); 
      for(int j=0; j<4; j++)
	{
	  connect[i]->Elem(j) = -1;
	}
      faces[i/patches[0]->GetFaceSize()]->Elem(i % patches[0]->GetFaceSize()) = 0;
      for (int j=0; j<i; j++)
	{
	  if (temp1[i] - temp1[j] == 0)
	    {
	      connect[i]->Elem(0) = i/patches[0]->GetFaceSize();
	      connect[i]->Elem(1) = j/patches[0]->GetFaceSize();
	      connect[i]->Elem(2) = i % patches[0]->GetFaceSize();
	      connect[i]->Elem(3) = j % patches[0]->GetFaceSize();
	      faces[i/patches[0]->GetFaceSize()]->Elem(i % patches[0]->GetFaceSize()) = 1;
	      counter -= 1;
	      break;
	    }
	}
      if ( faces[i/patches[0]->GetFaceSize()]->Elem(i % patches[0]->GetFaceSize()) < 1)
	{
	  for (int j=i+1; j<faces.Size()*patches[0]->GetFaceSize(); j++)
	    {
	      if (temp1[i] - temp1[j] == 0)
		{
		  connect[i]->Elem(0) = i/patches[0]->GetFaceSize();
		  connect[i]->Elem(1) = j/patches[0]->GetFaceSize();
		  connect[i]->Elem(2) = i % patches[0]->GetFaceSize();
		  connect[i]->Elem(3) = j % patches[0]->GetFaceSize();
		  faces[i/patches[0]->GetFaceSize()]->Elem(i % patches[0]->GetFaceSize()) = 1;
		  counter -= 1;
		  break;
		}
	    }
	}
    }
  NumOfFaces = counter;
  delete []temp1;

   
  for (int i=0; i<faces.Size(); i++)
    {
      for (int j=0; j<faces[0]->Size(); j++)
	{
	  cout << faces[i]->Elem(j) << " ";
	}
    }
   cout << endl << " wokkka wokka " << endl;
  for (int i=0; i<4; i++)
    {
      for (int j=0; j<connect.Size(); j++)
	{
	  cout << connect[j]->Elem(i) << " " ;
	}
      cout << endl;
    }
}

//============================================================================

void Geometry::PrintNURBS(int *n)
{
  if (PatchDim==2 && SpaceDim==2)
    {
      for (int i=0; i<NumOfPatches; i++)
	{
	  char *filename = new char[256];
	  sprintf(filename, "2D-NURBSsurface_P%d.out",i);
	  ofstream fileout(filename);
	  fileout << setprecision(5) << setiosflags(ios::scientific | ios::showpos);
	  double h = 1.0/n[i];
	  double x;
	  double y;
	  gsl_vector *B1, *B2;
	  gsl_bspline_workspace *bw1, *bw2;
	  gsl_vector *cknots1, *cknots2;
	  int xsize = patches[i]->GetKnotSize(0)+GetPatchDegs(i,0)-2;
	  int ysize = patches[i]->GetKnotSize(1)+GetPatchDegs(i,1)-2;
  	  B1 = gsl_vector_alloc(xsize);
	  B2 = gsl_vector_alloc(ysize);
	  cknots1 = gsl_vector_alloc(patches[i]->GetKnotSize(0));
	  cknots2 = gsl_vector_alloc(patches[i]->GetKnotSize(1));
	  for (int j=0; j<patches[i]->GetKnotSize(0); j++)
	    {
	      gsl_vector_set(cknots1,j,patches[i]->GetKnot(0,j));
	    }
	  for (int j=0; j<patches[i]->GetKnotSize(1); j++)
	    {
	      gsl_vector_set(cknots2,j,patches[i]->GetKnot(1,j));
	    }
	  bw1 = gsl_bspline_alloc(GetPatchDegs(i,0),patches[i]->GetKnotSize(0)); 
	  bw2 = gsl_bspline_alloc(GetPatchDegs(i,1),patches[i]->GetKnotSize(1));
	  gsl_bspline_knots(cknots1,bw1);
	  gsl_bspline_knots(cknots2,bw2);
	  int outcount = 0;
	  for (int j=0; j<n[i]+1; j++)
	    {
	      gsl_bspline_eval(j*h,B2,bw2);
	      for (int kk=0; kk<n[i]+1; kk++)
		{ 
		  gsl_bspline_eval(kk*h,B1,bw1);
		  x = 0;
		  y = 0;
		  int counter = 0;
		  double temp1 = 0;
		  for (int k=0; k<xsize; k++)
		    {
		      for (int l=0; l<ysize; l++)
			{
			  double temp = gsl_vector_get(B1,k) * gsl_vector_get(B2,l)*GetWeight(i,counter);
			  temp1 += temp;
			  x += temp * GetCtrlPt(i,0,counter);
			  y += temp * GetCtrlPt(i,1,counter);
			  counter += 1;
			}
		    }
		  fileout << kk*h  << "  " << j*h << "  " << x/temp1 << "  " << y/temp1 << endl; 
		  outcount += 1;
		}
	    }
	  gsl_vector_free(B1);
	  gsl_bspline_free(bw1);
	  gsl_vector_free(cknots1);
	  gsl_vector_free(B2);
	  gsl_bspline_free(bw2);
	  gsl_vector_free(cknots2);
	  fileout.close();
	  delete []filename;
	}
    }
 else if(PatchDim==2 && SpaceDim==3)
    {
      Array<double *> val(2);
      val[0] = new double[(n[0]+1)*(n[0]+1)];
      val[1] = new double[(n[0]+1)*(n[0]+1)];
      for (int i=0; i<NumOfPatches; i++)
	{
	  char *filename = new char[256];
	  sprintf(filename, "3D-NURBSsurface_P%d.vtk",i);
	  ofstream out(filename);
	  out << "# vtk DataFile Version 2.0" << endl;
	  out << "Unstructured Grid" << endl;
	  out << "ASCII" << endl;
	  out << "DATASET UNSTRUCTURED_GRID" << endl;
	  out << "POINTS " << (n[i]+1)*(n[i]+1) << " double" << endl;
	  double h = 1.0/n[i];
	  double x;
	  double y;
	  double z;
	  gsl_vector *B1, *B2;
	  gsl_bspline_workspace *bw1, *bw2;
	  gsl_vector *cknots1, *cknots2;
	  int xsize = patches[i]->GetKnotSize(0)+GetPatchDegs(i,0)-2;
	  int ysize = patches[i]->GetKnotSize(1)+GetPatchDegs(i,1)-2;
	  cout << "nas" << i << "  " << xsize << "  " << ysize << endl;
  	  B1 = gsl_vector_alloc(xsize);
	  B2 = gsl_vector_alloc(ysize);
	  cknots1 = gsl_vector_alloc(patches[i]->GetKnotSize(0));
	  cknots2 = gsl_vector_alloc(patches[i]->GetKnotSize(1));
	  for (int j=0; j<patches[i]->GetKnotSize(0); j++)
	    {
	      gsl_vector_set(cknots1,j,patches[i]->GetKnot(0,j));
	    }
	  for (int j=0; j<patches[i]->GetKnotSize(1); j++)
	    {
	      gsl_vector_set(cknots2,j,patches[i]->GetKnot(1,j));
	    }
	  bw1 = gsl_bspline_alloc(GetPatchDegs(i,0),patches[i]->GetKnotSize(0)); 
	  bw2 = gsl_bspline_alloc(GetPatchDegs(i,1),patches[i]->GetKnotSize(1));
	  gsl_bspline_knots(cknots1,bw1);
	  gsl_bspline_knots(cknots2,bw2);
	  int outcount = 0;
	  for (int j=0; j<n[i]+1; j++)
	    {
	     gsl_bspline_eval(j*h,B2,bw2);
	      for (int kk=0; kk<n[i]+1; kk++)
		{ 
		  gsl_bspline_eval(kk*h,B1,bw1);
		  x = 0;
		  y = 0;
		  z = 0;
		  int counter = 0;
		  double temp1 = 0;
		  for (int k=0; k<xsize; k++)
		    {
		      for (int l=0; l<ysize; l++)
			{
			  //cout << "weight " << GetWeight(i,counter) << endl;
			  double temp = gsl_vector_get(B1,k) * gsl_vector_get(B2,l);// * GetWeight(i,counter);
			  double temp2 = temp * GetWeight(i,counter);
			  //cout << "weight " << GetCtrlPt(i,0,counter) << "  "  << GetCtrlPt(i,1,counter) << "  " << GetCtrlPt(i,2,counter) << "  " << "  " << GetWeight(i,counter)  << endl;
			  temp1 += temp2;
			  x += temp * GetCtrlPt(i,0,counter);
			  y += temp * GetCtrlPt(i,1,counter);
			  z += temp * GetCtrlPt(i,2,counter);
			  counter += 1;
			}
		    }
		  out << x/temp1 << "  " << y/temp1 << "  " << z/temp1 << endl; 
		  val[i][outcount] = temp1;		  
		  outcount += 1;
		  //cout << "bspline " << x << "  " << y << " " << z << endl;
		}
	    }
	  gsl_vector_free(B1);
	  gsl_bspline_free(bw1);
	  gsl_vector_free(cknots1);
	  gsl_vector_free(B2);
	  gsl_bspline_free(bw2);
	  gsl_vector_free(cknots2);

	  out << "CELLS " << n[i]*n[i] << " " << 5*n[i]*n[i] << endl;
	  for (int j=0; j<n[i]; j++)
	    {
	      for (int k=0; k<n[i]; k++)
		{
		  int place = (n[i]+1)*j+k;
		  out << 4 << " " << place << "  " << place+1 << "  " << place+n[i]+2 << "  " << place+n[i]+1 << endl;
		}
	    }
	  out << "CELL_TYPES " << n[i]*n[i] << endl;
	  for (int j=0; j<n[i]*n[i]; j++)
	    {
	      out << 9 << endl;
	    } 
	  out.close();
	  delete []filename;
	}
      double diff=0;
      for (int i=0; i <(n[0]+1)*(n[0]+1);i++)
	{
	  // cout << " ahsdf " <<  << endl;
	  diff += sqrt((val[0][i]- val[1][i])*(val[0][i]- val[1][i]));
	}
      cout << "diff is " << diff << endl;
      delete []val[0];
      delete []val[1];
    }
 else if (PatchDim == 3 && SpaceDim == 3)
    {
      for (int i=0; i<NumOfPatches; i++)
	{
	  char *filename = new char[256];
	  sprintf(filename, "3D-NURBSvolume_P%d.vtk",i);
	  ofstream out(filename);
	  out << "# vtk DataFile Version 2.0" << endl;
	  out << "Unstructured Grid" << endl;
	  out << "ASCII" << endl;
	  out << "DATASET UNSTRUCTURED_GRID" << endl;
	  out << "POINTS " << (n[i]+1)*(n[i]+1)*(n[i]+1) << " double" << endl;
	  double h = 1.0/n[i];
	  double x;
	  double y;
	  double z;
	  gsl_vector *B1, *B2, *B3;
	  gsl_bspline_workspace *bw1, *bw2, *bw3;
	  gsl_vector *cknots1, *cknots2, *cknots3;
	  int xsize = patches[i]->GetKnotSize(0)+GetPatchDegs(i,0)-2;
	  int ysize = patches[i]->GetKnotSize(1)+GetPatchDegs(i,1)-2;
	  int zsize = patches[i]->GetKnotSize(2)+GetPatchDegs(i,2)-2;
  	  B1 = gsl_vector_alloc(xsize);
	  B2 = gsl_vector_alloc(ysize);
	  B3 = gsl_vector_alloc(zsize);
	  cknots1 = gsl_vector_alloc(patches[i]->GetKnotSize(0));
	  cknots2 = gsl_vector_alloc(patches[i]->GetKnotSize(1));
	  cknots3 = gsl_vector_alloc(patches[i]->GetKnotSize(2));
	  for (int j=0; j<patches[i]->GetKnotSize(0); j++)
	    {
	      gsl_vector_set(cknots1,j,patches[i]->GetKnot(0,j));
	    }
	  for (int j=0; j<patches[i]->GetKnotSize(1); j++)
	    {
	      gsl_vector_set(cknots2,j,patches[i]->GetKnot(1,j));
	    }
	  for (int j=0; j<patches[i]->GetKnotSize(2); j++)
	    {
	      gsl_vector_set(cknots3,j,patches[i]->GetKnot(2,j));
	    }
	  bw1 = gsl_bspline_alloc(GetPatchDegs(i,0),patches[i]->GetKnotSize(0)); 
	  bw2 = gsl_bspline_alloc(GetPatchDegs(i,1),patches[i]->GetKnotSize(1));
	  bw3 = gsl_bspline_alloc(GetPatchDegs(i,2),patches[i]->GetKnotSize(2));
	  gsl_bspline_knots(cknots1,bw1);
	  gsl_bspline_knots(cknots2,bw2);
	  gsl_bspline_knots(cknots3,bw3);
	  int outcount = 0;
	  for (int j=0; j<n[i]+1; j++)
	    {
	     gsl_bspline_eval(j*h,B3,bw3);
	      for (int kk=0; kk<n[i]+1; kk++)
		{ 
		  gsl_bspline_eval(kk*h,B2,bw2);
		  for (int ll=0; ll<n[i]+1; ll++)
		    {
		      gsl_bspline_eval(ll*h,B1,bw1);
		      x = 0;
		      y = 0;
		      z = 0;
		      int counter = 0;
		      double temp1 = 0;
		      for (int k=0; k<zsize; k++)
			{
			  for (int l=0; l<xsize; l++)
			    {
			      for (int jj=0; jj<ysize; jj++)
				{
				  double temp = gsl_vector_get(B2,jj)*gsl_vector_get(B3,k) * gsl_vector_get(B1,l) * GetWeight(i,counter);
				  temp1 += temp;
				  x += temp * GetCtrlPt(i,0,counter);
				  y += temp * GetCtrlPt(i,1,counter);
				  z += temp * GetCtrlPt(i,2,counter); 
				  counter += 1;
				}
			    }
			}
		      out << x/temp1 << "  " << y/temp1 << "  " << z/temp1 << endl; 
		      outcount += 1;
		    }
		}
	    }
	  gsl_vector_free(B1);
	  gsl_bspline_free(bw1);
	  gsl_vector_free(cknots1);
	  gsl_vector_free(B2);
	  gsl_bspline_free(bw2);
	  gsl_vector_free(cknots2);
	  gsl_vector_free(B3);
	  gsl_bspline_free(bw3);
	  gsl_vector_free(cknots3);

	  out << "CELLS " << n[i]*n[i]*n[i] << " " << 9*n[i]*n[i]*n[i] << endl;
	  for (int j=0; j<n[i]; j++)
	    {
	      for (int k=0; k<n[i]; k++)
		{
		  for (int l=0; l<n[i]; l++)
		    {
		      int place = (n[i]+1)*k+l + (n[i]+1)*(n[i]+1)*j;
		      out << 8 << " " << place << "  " << place+1 << "  " << place+n[i]+2 << "  " << place+n[i]+1 << "  " << place+(n[i]+1)*(n[i]+1) << "  " << place+1+(n[i]+1)*(n[i]+1) << "  " << place+n[i]+2+(n[i]+1)*(n[i]+1) << "  " << place+n[i]+1+(n[i]+1)*(n[i]+1) << endl;
		    }
		}
	    }
	  out << "CELL_TYPES " << n[i]*n[i]*n[i] << endl;
	  for (int j=0; j<n[i]*n[i]*n[i]; j++)
	    {
	      out << 12 << endl;
	    } 
	  out.close();
	  delete []filename;
	}

    }
  else
    {
      cout << " Batman says: Recheck input file...THEN SWEAR TO ME!!! " << endl;
      cout << "          _==/           |     |          \\==_       " << endl; 
      cout << "        /XX/             |\\___/|            \\XX\\	 " << endl;
      cout << "       /XXXX\\            |XXXXX|           /XXXX\\ " << endl;
      cout << "      |XXXXXX\\_         _XXXXXXX_         _/XXXXXX| " << endl;
      cout << "     XXXXXXXXXXXxxxxxxxXXXXXXXXXXXxxxxxxxXXXXXXXXXXX " << endl; 
      cout << "    |XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX| " << endl; 
      cout << "    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX " << endl; 
      cout << "    |XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX| " << endl; 
      cout << "     XXXXXX/^^^^\"\\XXXXXXXXXXXXXXXXXXXXX/^^^^^\\XXXXXX " << endl; 
      cout << "      |XXX|       \\XXX/^^\\XXXXX/^^\\XXX/       |XXX| " << endl;
      cout << "        \\XX\\       \\X/    \\XXX/    \\X/       /XX/ " << endl; 
      cout << "           \"\\       \"      \\X/      \"      /\" "      << endl;
    }
}

//============================================================================

void Geometry::PrintBSSurface(int *n)
{
  if (PatchDim==2 && SpaceDim==2)
    {
      for (int i=0; i<NumOfPatches; i++)
	{
	  char *filename = new char[256];
	  sprintf(filename, "2D-bsplinesurface_P%d.out",i);
	  ofstream fileout(filename);
	  fileout << setprecision(5) << setiosflags(ios::scientific | ios::showpos);
	  double h = 1.0/n[i];
	  double *x = new double[(n[i]+1)*(n[i]+1)];
	  double *y = new double[(n[i]+1)*(n[i]+1)];
	  gsl_vector *B1, *B2;
	  gsl_bspline_workspace *bw1, *bw2;
	  gsl_vector *cknots1, *cknots2;
	  int xsize = patches[i]->GetKnotSize(0)+GetPatchDegs(i,0)-2;
	  int ysize = patches[i]->GetKnotSize(1)+GetPatchDegs(i,1)-2;
  	  B1 = gsl_vector_alloc(xsize);
	  B2 = gsl_vector_alloc(ysize);
	  cknots1 = gsl_vector_alloc(patches[i]->GetKnotSize(0));
	  cknots2 = gsl_vector_alloc(patches[i]->GetKnotSize(1));
	  for (int j=0; j<patches[i]->GetKnotSize(0); j++)
	    {
	      gsl_vector_set(cknots1,j,patches[i]->GetKnot(0,j));
	    }
	  for (int j=0; j<patches[i]->GetKnotSize(1); j++)
	    {
	      gsl_vector_set(cknots2,j,patches[i]->GetKnot(1,j));
	    }
	  bw1 = gsl_bspline_alloc(GetPatchDegs(i,0),patches[i]->GetKnotSize(0)); 
	  bw2 = gsl_bspline_alloc(GetPatchDegs(i,1),patches[i]->GetKnotSize(1));
	  gsl_bspline_knots(cknots1,bw1);
	  gsl_bspline_knots(cknots2,bw2);
	  int outcount = 0;
	  for (int j=0; j<n[i]+1; j++)
	    {
	      for (int kk=0; kk<n[i]+1; kk++)
		{ 
		  gsl_bspline_eval(j*h,B1,bw1);
		  gsl_bspline_eval(kk*h,B2,bw2);
		  x[outcount] = 0;
		  y[outcount] = 0;
		  int counter = 0;
		  for (int k=0; k<xsize; k++)
		    {
		      for (int l=0; l<ysize; l++)
			{
			  double temp = gsl_vector_get(B1,k) * gsl_vector_get(B2,l);
			  x[outcount] += temp * GetCtrlPt(i,0,counter);
			  y[outcount] += temp * GetCtrlPt(i,1,counter);
			  counter += 1;
			}
		    }
		  fileout << j*h  << "  " << kk*h << "  " << x[outcount] << "  " << y[outcount] << endl; 
		  outcount += 1;
		}
	    }
	  delete []x;
	  delete []y;
	  gsl_vector_free(B1);
	  gsl_bspline_free(bw1);
	  gsl_vector_free(cknots1);
	  gsl_vector_free(B2);
	  gsl_bspline_free(bw2);
	  gsl_vector_free(cknots2);
	  fileout.close();
	  delete []filename;
	}
    }
}

//============================================================================

void Geometry::Refine(Array<int *> &factor)
{
  for (int w=0; w<NumOfPatches; w++)
    {
      Array<Vector *> temp(PatchDim); // array to house the positions of the new knots 
      Array<Vector *> rfd(PatchDim); // array to house the refined knots vector
      Array<Vector *> ctrlpts(SpaceDim);
      Array<Vector *> ctemp(PatchDim);
      Array<Vector *> weights(1);
      Array<Vector *> tknots(PatchDim);
      for (int i=0; i<PatchDim; i++)
	{
	  int counter = 1;
	  for (int j=0; j<patches[w]->GetKnotSize(i)-1; j++)
	    {
	      double h = (patches[w]->GetKnot(i,j+1)-patches[w]->GetKnot(i,j));
	      if(h > 0)
		{
		  for (int k=1; k<factor[w][i]; k++)
		    {
		      counter += 1;
		    }
		}
	      counter += 1 ;
	    }
	  rfd[i] = new Vector(counter);
	  temp[i] = new Vector(counter - patches[w]->GetKnotSize(i));
	  cout << "increase is " << counter - patches[w]->GetKnotSize(i) << endl;
	}
      for (int i=0; i<PatchDim; i++)
	{
	  int tc = 0;
	  int counter = 0;
	  for (int j=0; j<patches[w]->GetKnotSize(i)-1; j++)
	    {
	      rfd[i]->Elem(counter) = patches[w]->GetKnot(i,j);
	      double h = (patches[w]->GetKnot(i,j+1)-patches[w]->GetKnot(i,j));
	      if(h > 0)
		{
		  for (int k=1; k<factor[w][i]; k++)
		    {
		      counter += 1;
		      temp[i]->Elem(tc) = counter;
		      tc += 1;
		      rfd[i]->Elem(counter) = patches[w]->GetKnot(i,j) + h/factor[w][i]*k;
		    }
		}
	      counter += 1 ;
	    }
	  rfd[i]->Elem(counter) = patches[w]->GetKnot(i,patches[w]->GetKnotSize(i)-1);
	}
      int numgridpts = 1;
      for (int i=0; i<PatchDim; i++) 
	{numgridpts *= rfd[i]->Size() + patches[w]->GetDegs(i) - 2;}
      Array<Vector *> temppts(numgridpts); // each individual vector is a control point coupled with a vectorweight of dimension SpaceDim+1 
      cout << "numgrid pts " << numgridpts << endl;
      for (int i=0; i<numgridpts; i++)
	{
	  temppts[i] = new Vector(SpaceDim+1);
	}
      for (int i=0; i<patches[w]->GetCtrlPtSize(); i++)
	{
	  for (int j=0; j<SpaceDim; j++)
	    {
	      temppts[i]->Elem(j) = patches[w]->GetCtrlPt(j,i)*patches[w]->GetWeight(i);
	    }
	  temppts[i]->Elem(SpaceDim) = patches[w]->GetWeight(i);
	}
      //-------------------------------------------
      //refine y-direction first
      //-------------------------------------------
      int sz = patches[w]->GetCtrlPtSize()/(patches[w]->GetKnotSize(1)+patches[w]->GetDegs(1)-2);
      Vector *knots = new Vector(patches[w]->GetKnotSize(1)); 
      for (int i=0; i<knots->Size(); i++)
	{
	  double a  = patches[w]->GetKnot(1,i);
	  knots->Elem(i) = a;
	}
      for (int i=0; i<temp[1]->Size(); i++)
	{
	  Array<Vector *> oldpts(patches[w]->GetCtrlPtSize() + i*sz);
	  for (int j=0; j<oldpts.Size(); j++)
	    {
	      oldpts[j] = new Vector(SpaceDim+1);
	      for (int k=0; k<SpaceDim+1; k++)
		{
		  double a = temppts[j]->Elem(k);
		  oldpts[j]->Elem(k) = a;
		}
	    }
	  Vector *oldknots = new Vector(knots->Size());
	  for (int j=0; j<knots->Size(); j++)
	    {
	      double a = knots->Elem(j);
	      oldknots->Elem(j) = a;
	    }
	  Vector *alpha = new Vector;
	  alphafill(w,1,temp[1]->Elem(i)-1,oldknots,rfd[1]->Elem(temp[1]->Elem(i)),alpha);
	  knots->SetSize(oldknots->Size()+1);
	  //update knots
	  for (int j=0; j<temp[1]->Elem(i); j++)
	    {
	      double a = oldknots->Elem(j);
	      knots->Elem(j) = a;
	    }
	  double tpt = rfd[1]->Elem(temp[1]->Elem(i));
	  knots->Elem(temp[1]->Elem(i)) = tpt;
	  for (int j=temp[1]->Elem(i)+1; j<knots->Size(); j++)
	    {
	      double a = oldknots->Elem(j-1);
	      knots->Elem(j) = a;
	    }
	  cout << "knots is [";
	  for (int j=0; j<knots->Size(); j++)
	    {
	      cout << knots->Elem(j) << " " ;
	    }
		 cout << "]" << endl;
	  //
	  //update ctrlpoints
	  int counter = 0;
	  for (int j=0; j<sz; j++)
	    {
	      add(*oldpts[counter-j],0,*oldpts[counter-j],*temppts[counter]);
	      for (int k=1; k<alpha->Size()-1; k++)
		{
		  counter += 1;
		  Vector *tempvec = new Vector(SpaceDim+1);
		  for (int l=0; l<SpaceDim+1;l++)
		    {
		      double a = oldpts[counter-j]->Elem(l) * alpha->Elem(k);
		      tempvec->Elem(l) = a;
		    }
		  add(*tempvec,1.0-alpha->Elem(k),*oldpts[counter-j-1],*temppts[counter]);
		  delete tempvec;
		}
	      counter += 1;
	      add(*oldpts[counter-j-1],0,*oldpts[counter-j-1],*temppts[counter]);
	      counter += 1;
	    }		    
	  delete alpha;
	  delete oldknots;
	  for (int j=0; j<oldpts.Size(); j++)
	    {
	      delete oldpts[j];
	    }
	}
      int offset = knots->Size()+patches[w]->GetDegs(1)-2;
            
      //-------------------------------------------
      //refine x-direction second
      //-------------------------------------------
      knots->SetSize(patches[w]->GetKnotSize(0));
      sz = (patches[w]->GetCtrlPtSize()+temp[1]->Size()*sz)/(knots->Size()+patches[w]->GetDegs(0)-2);
      cout << sz << endl;
      for (int i=0; i<knots->Size(); i++)
	{
	  double a  = patches[w]->GetKnot(0,i);
	  knots->Elem(i) = a;
	}
      int increase = temp[1]->Size()*(patches[w]->GetKnotSize(0)+patches[w]->GetDegs(0)-2);
      for (int i=0; i<temp[0]->Size(); i++)
	{
	  Array<Vector *> oldpts(patches[w]->GetCtrlPtSize() + i*sz + increase);
	  for (int j=0; j<oldpts.Size(); j++)
	    {
	      oldpts[j] = new Vector(SpaceDim+1);
	      for (int k=0; k<SpaceDim+1; k++)
		{
		  double a = temppts[j]->Elem(k);
		  oldpts[j]->Elem(k) = a;
		}
	    }
	  Vector *oldknots = new Vector(knots->Size());
	  for (int j=0; j<knots->Size(); j++)
	    {
	      double a = knots->Elem(j);
	      oldknots->Elem(j) = a;
	    }
	  Vector *alpha = new Vector;
	  cout << "the position " << temp[0]->Elem(i)-1 << endl;
	  alphafill(w,0,temp[0]->Elem(i)-1,oldknots,rfd[0]->Elem(temp[0]->Elem(i)),alpha);
	  knots->SetSize(oldknots->Size()+1);
	  //update knots
	  for (int j=0; j<temp[0]->Elem(i); j++)
	    {
	      double a = oldknots->Elem(j);
	      knots->Elem(j) = a;
	    }
	  double tpt = rfd[0]->Elem(temp[0]->Elem(i));
	  knots->Elem(temp[0]->Elem(i)) = tpt;
	  for (int j=temp[0]->Elem(i)+1; j<knots->Size(); j++)
	    {
	      double a =oldknots->Elem(j-1);
	      knots->Elem(j) = a;
	    }
	  //int counter = 0;
	  for (int j=0; j<sz; j++)
	    {
	      add(*oldpts[j],0,*oldpts[0],*temppts[j]);
	      for (int k=1; k<alpha->Size()-1; k++)
		{
		  //counter += 1;
		  Vector *tempvec = new Vector(SpaceDim+1);
		  for (int l=0; l<SpaceDim+1; l++)
		    {
		      double a = oldpts[j+k*offset]->Elem(l) * alpha->Elem(k);
		      tempvec->Elem(l) = a;
		    }
		  add(*tempvec,1.0-alpha->Elem(k),*oldpts[j+(k-1)*offset],*temppts[j+k*offset]);
		  delete tempvec;
		}
	      //counter += 1;
	      add(*oldpts[j+(alpha->Size()-2)*offset],0,*oldpts[0],*temppts[j+(alpha->Size()-1)*offset]);
	      // counter += 1;
	    }
	  for (int j=0; j< alpha->Size(); j++)
	    {
	        cout << "alpha " << alpha->Elem(j) << endl;
	    }	 
	  delete alpha;
	  delete oldknots;
	  for (int j=0; j<oldpts.Size(); j++)
	    {
	      delete oldpts[j];
	    }
	}
      for (int i=0; i<temppts.Size(); i++)
	{
	    cout << "x ref [" << temppts[i]->Elem(0) << " " <<  temppts[i]->Elem(1) << " " <<  temppts[i]->Elem(2) << " " <<  temppts[i]->Elem(3) << "]" << endl;
	}
      delete knots;
      /* if (w==0)
	{
	  for (int i=1; i<temppts.Size(); i++)
	    {
	      if (temppts[i]->Elem(SpaceDim)<1)
		{
		  temppts[i]->Elem(SpaceDim) *= 0.25;
		}
	    }
	    }*/
      patches[w]->Refine(rfd,temppts);
      for (int i=0; i<temppts.Size(); i++)
	{
	  delete temppts[i];
	}
      for (int i=0; i<PatchDim; i++)
	{
	  delete rfd[i];
	  delete temp[i];
	}
    }
}
	  
//============================================================================
void Geometry::alphafill(int w, int i, int k, Vector *knots, double newknot, Vector *alpha)
{
  // w -> patch
  // i -> dimension x, y, or z
  // k -> knot index preceding newknot
  // knots -> knots vector to be updated
  // newknot -> value of knot value used to update knots (must be => kth entry of knots)
  alpha->SetSize(knots->Size()+patches[w]->GetDegs(i)-1);
  Vector *tknots = new Vector(knots->Size()+2*patches[i]->GetDegs(i)-2);
  for (int j=0; j<patches[w]->GetDegs(i)-1; j++)
    {
      tknots->Elem(j) = 0;
      tknots->Elem(tknots->Size()-(j+1)) = 1.0;
    }
  for (int j=0; j<knots->Size(); j++)
    {
      double a = knots->Elem(j);
      tknots->Elem(j+patches[w]->GetDegs(i)-1) = a;
    }
  for (int j=0; j<tknots->Size(); j++)
    {
      // cout << tknots->Elem(j) << "  " ;
    }
  for (int j=0;  j<alpha->Size(); j++)
    {
      if (j <= k )
	{
	  // cout << " ppants " << newknot << endl;
	 alpha->Elem(j) = 1.0;
	}
      else if (j <= k+patches[w]->GetDegs(i)-1 )
	{
	  double a = (newknot - tknots->Elem(j))/(tknots->Elem(j+patches[w]->GetDegs(i)-1) - tknots->Elem(j));
	  alpha->Elem(j) = a;
	  //cout << "opants " << tknots->Elem(j+patches[w]->GetDegs(i)-1) - tknots->Elem(j) << " , " << newknot - tknots->Elem(j) << endl;
	}
      else
	{
	  //cout <<"qpants" << endl;
	  alpha->Elem(j) = 0;
	}
    }
  delete tknots;
 }

//============================================================================

Geometry::~Geometry()
{
  for (int i=0; i<NumOfPatches; i++)
    {
      delete patches[i];
      delete faces[i]; 
    }
  for (int i=0; i<connect.Size(); i++)
    {
      delete connect[i];
    }
}
