#include "timedependent_header.h"

//============================================================================

double MinMod(double a, double b)
{
  double ma = fabs(a);
  double mb = fabs(b);
  double min_a_b = (ma<mb) ? ma : mb;
  double signa = ((a<0.0) ? -1 : 1);
  double signb = ((b<0.0) ? -1 : 1);
  double fact = (double) (signa + signb);
  return (0.5 * fact * min_a_b);
} 

double SlopeInfo_0(int i, Mesh *mesh, Vector &sol)
{
  
  int nx = mesh->GetNx();
  int l, r;
  l = (i%(nx+1)==0) ? i : i-1;
  r = (i%(nx+1)==nx) ? i : i+1;
  double a = sol(r) - sol(i);
  double b = sol(i) - sol(l);
  return MinMod(a, b);
  
  /*
  int nx = mesh->GetNx();
  int ny = mesh->GetNE()/nx; ny = ny/2;
  int l, r;
  l = (i%(nx+1)==0) ? i : i-1;
  r = (i%(nx+1)==nx) ? i : i+1;
  int t, b;
  t = (i>=ny*(nx+1)) ? i : i+(nx+1);
  b = (i<=nx) ? i : i-(nx+1);

  if(i==-1) //((i>=ny*(nx+1))||(i<=nx))
    {
      double aa = sol(r) - sol(i);
      double bb = sol(i) - sol(l);
      return MinMod(aa, bb);
    }
  else
    {
      double aa = (2.0*sol(r)+sol(t))/3.0 - sol(i);
      double bb = sol(i) - (2.0*sol(l)+sol(b))/3.0;
      return 6.0*MinMod(aa, bb)/5.0;
    }
  */
}

double SlopeInfo_1(int i, Mesh *mesh, Vector &sol)
{
  
  int nx = mesh->GetNx();
  int ny = mesh->GetNE()/nx; ny = ny/2;
  int t, b;
  t = (i>=ny*(nx+1)) ? i : i+(nx+1);
  b = (i<=nx) ? i : i-(nx+1);
  double aa = sol(t) - sol(i);
  double bb = sol(i) - sol(b);
  return MinMod(aa, bb);
  
  /*
  int nx = mesh->GetNx();
  int ny = mesh->GetNE()/nx; ny = ny/2;
  int l, r;
  l = (i%(nx+1)==0) ? i : i-1;
  r = (i%(nx+1)==nx) ? i : i+1;
  int t, b;
  t = (i>=ny*(nx+1)) ? i : i+(nx+1);
  b = (i<=nx) ? i : i-(nx+1);

  if(i==-1) //((i%(nx+1)==0)||(i%(nx+1)==nx))
    {
      double aa = sol(t) - sol(i);
      double bb = sol(i) - sol(b);
      return MinMod(aa, bb);
    }
  else
    {
      double aa = (2.0*sol(t)+sol(r))/3.0 - sol(i);
      double bb = sol(i) - (2.0*sol(b)+sol(l))/3.0;
      return 6.0*MinMod(aa, bb)/5.0;
    }
  */
}

double SlopeInfo_2(int i, Mesh *mesh, Vector &sol)
{
  int nx = mesh->GetNx();
  int ny = mesh->GetNE()/nx; ny = ny/2;
  int lt, br;
  lt = ((i%(nx+1)==0)||(i>=ny*(nx+1))) ? i : i+nx;
  br = ((i%(nx+1)==nx)||(i<=nx)) ? i : i-nx;
  double a = sol(br) - sol(i);
  double b = sol(i) - sol(lt);
  return MinMod(a, b);
}

//============================================================================

void LocalLinearTriangularEUF(Array<double> &locflux, Array<double> &locsat, 
			      Function *fff, Array<double> &locupwindflux)
{
  Array<double> tem(3);
  double co[3][1];
  for(int i=0; i<3; i++) co[i][0] = locsat[i];
  
  if(locflux[0] > 0)
    tem[0] = fff->Eval(co[0]);
  else
    tem[0] = fff->Eval(co[1]);

  if(locflux[1] > 0)
    tem[1] = fff->Eval(co[1]);
  else
    tem[1] = fff->Eval(co[2]);

  if(locflux[2] > 0)
    tem[2] = fff->Eval(co[2]);
  else
    tem[2] = fff->Eval(co[0]);

  locupwindflux[0] = tem[0] * locflux[0] - tem[2] * locflux[2];
  locupwindflux[1] = tem[1] * locflux[1] - tem[0] * locflux[0];
  locupwindflux[2] = tem[2] * locflux[2] - tem[1] * locflux[1];
}

void LocalLinearTriangular2ndOrderEUF(int i, Mesh *mesh, Array<double> &locflux, Vector &sat, 
				      Function *fff, Array<double> &locupwindflux)
{
  Array<int> ind;
  mesh->GetElementVertices(i, ind);

  Array<double> tem(3);
  double val[1], sig[1];
  if(i%2==0)
    {
      if(locflux[0] > 0)
	{
	  val[0] = sat(ind[0]);
	  sig[0] = 0.5 * SlopeInfo_0(ind[0], mesh, sat);
	  val[0] += sig[0];
	}
      else
	{
	  val[0] = sat(ind[1]);
	  sig[0] = 0.5 * SlopeInfo_0(ind[1], mesh, sat);
	  val[0] -= sig[0];
	}
      tem[0] = locflux[0] * fff->Eval(val);
  
      if(locflux[1] > 0)
	{
	  val[0] = sat(ind[1]);
	  sig[0] = 0.5 * SlopeInfo_2(ind[1], mesh, sat);
	  val[0] -= sig[0];
	}
      else
	{
	  val[0] = sat(ind[2]);
	  sig[0] = 0.5 * SlopeInfo_2(ind[2], mesh, sat);
	  val[0] += sig[0];
	}
      tem[1] = locflux[1] * fff->Eval(val);
  
      if(locflux[2] > 0)
	{
	  val[0] = sat(ind[2]);
	  sig[0] = 0.5 * SlopeInfo_1(ind[2], mesh, sat);
	  val[0] -= sig[0];
	}
      else
	{
	  val[0] = sat(ind[0]);
	  sig[0] = 0.5 * SlopeInfo_1(ind[0], mesh, sat);
	  val[0] += sig[0];
	}
      tem[2] = locflux[2] * fff->Eval(val);
    }
  else
    {
      if(locflux[0] > 0)
	{
	  val[0] = sat(ind[0]);
	  sig[0] = 0.5 * SlopeInfo_1(ind[0], mesh, sat);
	  val[0] += sig[0];
	}
      else
	{
	  val[0] = sat(ind[1]);
	  sig[0] = 0.5 * SlopeInfo_1(ind[1], mesh, sat);
	  val[0] -= sig[0];
	}
      tem[0] = locflux[0] * fff->Eval(val);
  
      if(locflux[1] > 0)
	{
	  val[0] = sat(ind[1]);
	  sig[0] = 0.5 * SlopeInfo_0(ind[1], mesh, sat);
	  val[0] -= sig[0];
	}
      else
	{
	  val[0] = sat(ind[2]);
	  sig[0] = 0.5 * SlopeInfo_0(ind[2], mesh, sat);
	  val[0] += sig[0];
	}
      tem[1] = locflux[1] * fff->Eval(val);
  
      if(locflux[2] > 0)
	{
	  val[0] = sat(ind[2]);
	  sig[0] = 0.5 * SlopeInfo_2(ind[2], mesh, sat);
	  val[0] += sig[0];
	}
      else
	{
	  val[0] = sat(ind[0]);
	  sig[0] = 0.5 * SlopeInfo_2(ind[0], mesh, sat);
	  val[0] -= sig[0];
	}
      tem[2] = locflux[2] * fff->Eval(val);      
    }
  
  locupwindflux[0] = tem[0] - tem[2];
  locupwindflux[1] = tem[1] - tem[0];
  locupwindflux[2] = tem[2] - tem[1];
}


//============================================================================

void LocalQuadraticTriangularEUF(Array<double> &locflux, Array<double> &locsat, 
				 Function *fff, Array<double> &locupwindflux)
{
  Array<double> tem(9);
  double co[6][1];
  for(int i=0; i<6; i++) co[i][0] = locsat[i];

  if(locflux[0] > 0)
    tem[0] = fff->Eval(co[0]);
  else
    tem[0] = fff->Eval(co[3]);
  if(locflux[1] > 0)
    tem[1] = fff->Eval(co[3]);
  else
    tem[1] = fff->Eval(co[5]);
  if(locflux[2] > 0)
    tem[2] = fff->Eval(co[5]);
  else
    tem[2] = fff->Eval(co[0]);

  if(locflux[3] > 0)
    tem[3] = fff->Eval(co[1]);
  else
    tem[3] = fff->Eval(co[4]);
  if(locflux[4] > 0)
    tem[4] = fff->Eval(co[4]);
  else
    tem[4] = fff->Eval(co[3]);
  if(locflux[5] > 0)
    tem[5] = fff->Eval(co[3]);
  else
    tem[5] = fff->Eval(co[1]);

  if(locflux[6] > 0)
    tem[6] = fff->Eval(co[2]);
  else
    tem[6] = fff->Eval(co[5]);
  if(locflux[7] > 0)
    tem[7] = fff->Eval(co[5]);
  else
    tem[7] = fff->Eval(co[4]);
  if(locflux[8] > 0)
    tem[8] = fff->Eval(co[4]);
  else
    tem[8] = fff->Eval(co[2]);
  

  locupwindflux[0] = tem[0] * locflux[0] - tem[2] * locflux[2];
  locupwindflux[1] = tem[3] * locflux[3] - tem[5] * locflux[5];
  locupwindflux[2] = tem[6] * locflux[6] - tem[8] * locflux[8];
  locupwindflux[3] = tem[1] * locflux[1] - tem[0] * locflux[0] + tem[5] * locflux[5] - tem[4] * locflux[4];
  locupwindflux[4] = tem[4] * locflux[4] - tem[3] * locflux[3] + tem[8] * locflux[8] - tem[7] * locflux[7];
  locupwindflux[5] = tem[2] * locflux[2] - tem[1] * locflux[1] + tem[7] * locflux[7] - tem[6] * locflux[6];
}

void LocalQuadraticTriangular2ndOrderEUF(int i, Mesh *mesh, Array<double> &locflux, Vector &sat,
					 Vector &SlopeInfo_0, Vector &SlopeInfo_1, Vector &SlopeInfo_2,
					 Function *fff, Array<double> &locupwindflux)
{
  Array<int> ind;
  mesh->GetElementVertices(i, ind);
  Array<int> edges;
  mesh->GetElementEdges(i, edges);
  for (int k = 0; k < edges.Size(); k++)
    ind.Append (mesh->GetNV() + edges[k]);

  Array<double> tem(9);
  double val[1];
  if(i%2==0)
    {
      if(locflux[0] > 0)
	val[0] = sat(ind[0]) + 0.5 * SlopeInfo_0(ind[0]);
      else
	val[0] = sat(ind[3]) - 0.5 * SlopeInfo_0(ind[3]);
      tem[0] = locflux[0] * fff->Eval(val);
  
      if(locflux[1] > 0)
	val[0] = sat(ind[3]) - 0.5 * SlopeInfo_2(ind[3]);
      else
	val[0] = sat(ind[5]) + 0.5 * SlopeInfo_2(ind[5]);
      tem[1] = locflux[1] * fff->Eval(val);
  
      if(locflux[2] > 0)
	val[0] = sat(ind[5]) - 0.5 * SlopeInfo_1(ind[5]);
      else
	val[0] = sat(ind[0]) + 0.5 * SlopeInfo_1(ind[0]);
      tem[2] = locflux[2] * fff->Eval(val);

      if(locflux[3] > 0)
	val[0] = sat(ind[1]) - 0.5 * SlopeInfo_2(ind[1]);
      else
	val[0] = sat(ind[4]) + 0.5 * SlopeInfo_2(ind[4]);
      tem[3] = locflux[3] * fff->Eval(val);
  
      if(locflux[4] > 0)
	val[0] = sat(ind[4]) - 0.5 * SlopeInfo_1(ind[4]);
      else
	val[0] = sat(ind[3]) + 0.5 * SlopeInfo_1(ind[3]);
      tem[4] = locflux[4] * fff->Eval(val);
  
      if(locflux[5] > 0)
	val[0] = sat(ind[3]) + 0.5 * SlopeInfo_0(ind[3]);
      else
	val[0] = sat(ind[1]) - 0.5 * SlopeInfo_0(ind[1]);
      tem[5] = locflux[5] * fff->Eval(val);

      if(locflux[6] > 0)
	val[0] = sat(ind[2]) - 0.5 * SlopeInfo_1(ind[2]);
      else
	val[0] = sat(ind[5]) + 0.5 * SlopeInfo_1(ind[5]);
      tem[6] = locflux[6] * fff->Eval(val);
  
      if(locflux[7] > 0)
	val[0] = sat(ind[5]) + 0.5 * SlopeInfo_0(ind[5]);
      else
	val[0] = sat(ind[4]) - 0.5 * SlopeInfo_0(ind[4]);
      tem[7] = locflux[7] * fff->Eval(val);
  
      if(locflux[8] > 0)
	val[0] = sat(ind[4]) - 0.5 * SlopeInfo_2(ind[4]);
      else
	val[0] = sat(ind[2]) + 0.5 * SlopeInfo_2(ind[2]);
      tem[8] = locflux[8] * fff->Eval(val);      
    }
  else
    {
      if(locflux[0] > 0)
	val[0] = sat(ind[0]) + 0.5 * SlopeInfo_1(ind[0]);
      else
	val[0] = sat(ind[3]) - 0.5 * SlopeInfo_1(ind[3]);
      tem[0] = locflux[0] * fff->Eval(val);
  
      if(locflux[1] > 0)
	val[0] = sat(ind[3]) - 0.5 * SlopeInfo_0(ind[3]);
      else
	val[0] = sat(ind[5]) + 0.5 * SlopeInfo_0(ind[5]);
      tem[1] = locflux[1] * fff->Eval(val);
  
      if(locflux[2] > 0)
	val[0] = sat(ind[5]) + 0.5 * SlopeInfo_2(ind[5]);
      else
	val[0] = sat(ind[0]) - 0.5 * SlopeInfo_2(ind[0]);
      tem[2] = locflux[2] * fff->Eval(val);

      if(locflux[3] > 0)
	val[0] = sat(ind[1]) - 0.5 * SlopeInfo_0(ind[1]);
      else
	val[0] = sat(ind[4]) + 0.5 * SlopeInfo_0(ind[4]);
      tem[3] = locflux[3] * fff->Eval(val);
  
      if(locflux[4] > 0)
	val[0] = sat(ind[4]) + 0.5 * SlopeInfo_2(ind[4]);
      else
	val[0] = sat(ind[3]) - 0.5 * SlopeInfo_2(ind[3]);
      tem[4] = locflux[4] * fff->Eval(val);
  
      if(locflux[5] > 0)
	val[0] = sat(ind[3]) + 0.5 * SlopeInfo_1(ind[3]);
      else
	val[0] = sat(ind[1]) - 0.5 * SlopeInfo_1(ind[1]);
      tem[5] = locflux[5] * fff->Eval(val);

      if(locflux[6] > 0)
	val[0] = sat(ind[2]) + 0.5 * SlopeInfo_2(ind[2]);
      else
	val[0] = sat(ind[5]) - 0.5 * SlopeInfo_2(ind[5]);
      tem[6] = locflux[6] * fff->Eval(val);
  
      if(locflux[7] > 0)
	val[0] = sat(ind[5]) + 0.5 * SlopeInfo_1(ind[5]);
      else
	val[0] = sat(ind[4]) - 0.5 * SlopeInfo_1(ind[4]);
      tem[7] = locflux[7] * fff->Eval(val);
  
      if(locflux[8] > 0)
	val[0] = sat(ind[4]) - 0.5 * SlopeInfo_0(ind[4]);
      else
	val[0] = sat(ind[2]) + 0.5 * SlopeInfo_0(ind[2]);
      tem[8] = locflux[8] * fff->Eval(val);      
    }
  
  locupwindflux[0] = tem[0] - tem[2];
  locupwindflux[1] = tem[3] - tem[5];
  locupwindflux[2] = tem[6] - tem[8];
  locupwindflux[3] = tem[1] - tem[0] + tem[5] - tem[4];
  locupwindflux[4] = tem[4] - tem[3] + tem[8] - tem[7];
  locupwindflux[5] = tem[2] - tem[1] + tem[7] - tem[6];
}

//============================================================================

void LocalCubicTriangularEUF(Array<double> &locflux, Array<double> &locsat, 
			     Function *fff, Array<double> &locupwindflux)
{
  Array<double> tem(18);
  double co[10][1];
  for(int i=0; i<10; i++) co[i][0] = locsat[i];
  
  if(locflux[0] > 0)
    tem[0] = fff->Eval(co[0]);
  else
    tem[0] = fff->Eval(co[3]);
  if(locflux[1] > 0)
    tem[1] = fff->Eval(co[8]);
  else
    tem[1] = fff->Eval(co[0]);
  if(locflux[2] > 0)
    tem[2] = fff->Eval(co[1]);
  else
    tem[2] = fff->Eval(co[5]);
  if(locflux[3] > 0)
    tem[3] = fff->Eval(co[4]);
  else
    tem[3] = fff->Eval(co[1]);
  if(locflux[4] > 0)
    tem[4] = fff->Eval(co[2]);
  else
    tem[4] = fff->Eval(co[7]);
  if(locflux[5] > 0)
    tem[5] = fff->Eval(co[6]);
  else
    tem[5] = fff->Eval(co[2]);

  if(locflux[6] > 0)
    tem[6] = fff->Eval(co[3]);
  else
    tem[6] = fff->Eval(co[4]);
  if(locflux[7] > 0)
    tem[7] = fff->Eval(co[5]);
  else
    tem[7] = fff->Eval(co[6]);
  if(locflux[8] > 0)
    tem[8] = fff->Eval(co[7]);
  else
    tem[8] = fff->Eval(co[8]);
  if(locflux[9] > 0)
    tem[9] = fff->Eval(co[3]);
  else
    tem[9] = fff->Eval(co[8]);
  if(locflux[10] > 0)
    tem[10] = fff->Eval(co[5]);
  else
    tem[10] = fff->Eval(co[4]);
  if(locflux[11] > 0)
    tem[11] = fff->Eval(co[7]);
  else
    tem[11] = fff->Eval(co[6]);

  if(locflux[12] > 0)
    tem[12] = fff->Eval(co[9]);
  else
    tem[12] = fff->Eval(co[3]);
  if(locflux[13] > 0)
    tem[13] = fff->Eval(co[4]);
  else
    tem[13] = fff->Eval(co[9]);
  if(locflux[14] > 0)
    tem[14] = fff->Eval(co[9]);
  else
    tem[14] = fff->Eval(co[5]);
  if(locflux[15] > 0)
    tem[15] = fff->Eval(co[6]);
  else
    tem[15] = fff->Eval(co[9]);
  if(locflux[16] > 0)
    tem[16] = fff->Eval(co[9]);
  else
    tem[16] = fff->Eval(co[7]);
  if(locflux[17] > 0)
    tem[17] = fff->Eval(co[8]);
  else
    tem[17] = fff->Eval(co[9]);
  
  locupwindflux[0] = tem[0] * locflux[0] - tem[1] * locflux[1];
  locupwindflux[1] = tem[2] * locflux[2] - tem[3] * locflux[3];
  locupwindflux[2] = tem[4] * locflux[4] - tem[5] * locflux[5];

  locupwindflux[3] = tem[9] * locflux[9] - tem[0] * locflux[0] + tem[6] * locflux[6] - tem[12] * locflux[12];
  locupwindflux[4] = tem[13] * locflux[13] - tem[6] * locflux[6] + tem[3] * locflux[3] - tem[10] * locflux[10];
  locupwindflux[5] = tem[10] * locflux[10] - tem[2] * locflux[2] + tem[7] * locflux[7] - tem[14] * locflux[14];
  locupwindflux[6] = tem[15] * locflux[15] - tem[7] * locflux[7] + tem[5] * locflux[5] - tem[11] * locflux[11];
  locupwindflux[7] = tem[11] * locflux[11] - tem[4] * locflux[4] + tem[8] * locflux[8] - tem[16] * locflux[16];
  locupwindflux[8] = tem[17] * locflux[17] - tem[8] * locflux[8] + tem[1] * locflux[1] - tem[9] * locflux[9];

  locupwindflux[9] = tem[12] * locflux[12] - tem[13] * locflux[13] + tem[14] * locflux[14]
    - tem[15] * locflux[15] + tem[16] * locflux[16] - tem[17] * locflux[17];
}

//============================================================================

void LocalLinearRectangularEUF(Array<double> &locflux, Array<double> &locsat, 
			       Function *fff, Array<double> &locupwindflux)
{
  ;
}

//============================================================================

void LocalQuadraticRectangularEUF(Array<double> &locflux, Array<double> &locsat, 
			       Function *fff, Array<double> &locupwindflux)
{
  ;
}

//============================================================================

void LocalCubicRectangularEUF(Array<double> &locflux, Array<double> &locsat, 
			       Function *fff, Array<double> &locupwindflux)
{
  ;
}

//============================================================================

void ElementUpwindFlux(DualMesh *dualmesh, Array<double*> &flux, Vector &sat, 
		       Function *fff, Array<double> &upwindflux)
{
  Mesh *mesh = dualmesh->GetMesh();
  int NE = mesh->GetNE();
  int order = dualmesh->GetDualMeshOrder();
  int dof = 0;
  if(order==1) dof = mesh->GetNV();
  else if(order==2) dof = mesh->GetNV() + mesh->GetNEdges();
  else dof = mesh->GetNV() + 2 * mesh->GetNEdges() + mesh->GetNE();
  
  int n;
  Array<int> ind;
  for(int i=0; i<NE; i++)
    {
      if(order==1)
	{
	  mesh->GetElementVertices(i, ind);
	  n = (mesh->GetMType()==Element::TRIANGLE) ? 3 : 4;
	  Array<double> locflux(n);
	  Array<double> locsat(n);
	  Array<double> locupwindflux(n);
	  for(int k=0; k<n; k++)
	    {
	      locflux[k] = flux[k][i];
	      locsat[k] = sat(ind[k]);
	      locupwindflux[k] = 0.0;
	    }
	  
	  if(mesh->GetMType()==Element::TRIANGLE)
	    {
	      //LocalLinearTriangularEUF(locflux, locsat, fff, locupwindflux);
	      LocalLinearTriangular2ndOrderEUF(i, mesh, locflux, sat, fff, locupwindflux);
	    }
	  else
	    LocalLinearRectangularEUF(locflux, locsat, fff, locupwindflux);

	  for(int k=0; k<n; k++)
	    {
	      upwindflux[ind[k]] += locupwindflux[k];
	    }
	}
      else if(order==2)
	{
	  mesh->GetElementVertices(i, ind);
	  Array<int> edges;
	  mesh->GetElementEdges(i, edges);
	  for (int k = 0; k < edges.Size(); k++)
	    ind.Append (mesh->GetNV() + edges[k]);
	  
	  n = (mesh->GetMType()==Element::TRIANGLE) ? 6 : 9;
	  Array<double> locflux(9);
	  Array<double> locsat(n);
	  Array<double> locupwindflux(n);
	  for(int k=0; k<n; k++)
	    {
	      locsat[k] = sat(ind[k]);
	      locupwindflux[k] = 0.0;
	    }
	  for(int k=0; k<9; k++)
	    {
	      locflux[k] = flux[k][i];
	    }
	  
	  if(mesh->GetMType()==Element::TRIANGLE)
	    {
	      LocalQuadraticTriangularEUF(locflux, locsat, fff, locupwindflux);
	      //LocalQuadraticTriangular2ndOrderEUF(i, mesh, locflux, sat, SI_0, SI_1, SI_2, fff, locupwindflux);
	    }
	  else
	    LocalQuadraticRectangularEUF(locflux, locsat, fff, locupwindflux);

	  for(int k=0; k<n; k++)
	    {
	      upwindflux[ind[k]] += locupwindflux[k];
	    }
	}
      else
	{
	  mesh->GetElementVertices(i, ind);
	  // === need more for ind ===
	  n = (mesh->GetMType()==Element::TRIANGLE) ? 10 : 16;
	  Array<double> locflux(n);
	  Array<double> locsat(n);
	  Array<double> locupwindflux(n);
	  for(int k=0; k<n; k++)
	    {
	      locflux[k] = flux[k][i];
	      locsat[k] = sat(ind[k]);
	      locupwindflux[k] = 0.0;
	    }
	  
	  if(mesh->GetMType()==Element::TRIANGLE)
	    LocalCubicTriangularEUF(locflux, locsat, fff, locupwindflux);
	  else
	    LocalCubicRectangularEUF(locflux, locsat, fff, locupwindflux);

	  for(int k=0; k<n; k++)
	    {
	      upwindflux[ind[k]] += locupwindflux[k];
	    }
	}
    }  
}

//============================================================================

void BdrElementUpwindFlux(Array<double> &bdrflux, Vector &sat, 
			  Function *fff, Array<double> &upwindbdrflux)
{ // no upwind at the global boundary
  double coord[1];
  
  for(int i=0; i<bdrflux.Size(); i++)
    {
      coord[0] = sat(i);
      upwindbdrflux[i] = bdrflux[i] * fff->Eval(coord);
      //cout << setprecision(9) << i << "\t" << upwindbdrflux[i] << "\t" << bdrflux[i] << endl;
    }
}

//============================================================================

void TPFTimeMarchSolve(double dt, int numfinetimestep, DualMesh *dualmesh, Array<double> &areas,
		       Array<double*> &flux, Array<double> &bdrflux,
		       Function *fff, Vector &satold, Vector &satnew)
{
  Mesh *mesh = dualmesh->GetMesh();
  int order = dualmesh->GetDualMeshOrder();
  int dof = 0;
  if(order==1) dof = mesh->GetNV();
  else if(order==2) dof = mesh->GetNV() + mesh->GetNEdges();
  else dof = mesh->GetNV() + 2 * mesh->GetNEdges() + mesh->GetNE();

  Array<double> upwindflux(dof);
  Array<double> upwindbdrflux(dof);
  
  for(int k=0; k<numfinetimestep; k++)
    {
      for(int j=0; j<upwindflux.Size(); j++)
	{
	  upwindflux[j] = 0.0;
	  upwindbdrflux[j] = 0.0;
	}

      ElementUpwindFlux(dualmesh, flux, satold, fff, upwindflux);
      BdrElementUpwindFlux(bdrflux, satold, fff, upwindbdrflux);
      for(int i=0; i<dof; i++)
	{
	  satnew(i) = satold(i) - dt * (upwindflux[i] + upwindbdrflux[i])  / areas[i];

	  satold = satnew;
	  
	  /*	  
	  if(satold(i)<-1.0e-10)
	    {
	      cout << setprecision(10) << i << "\t" << satold(i) << endl;
	      exit(1);
	    }

	  if(satold(i)>1+1.0e-10)
	    {
	      cout << i << "\t" << satold(i) << endl;
	      exit(1);
	    }
	  */	 
	}
    }
}


/*
void SlopeInfo(int i, Mesh *mesh, Vector &sol, Vector &SlopeInfo_0, Vector &SlopeInfo_1, Vector &SlopeInfo_2)
{
  SlopeInfo_0 = 0.0;
  SlopeInfo_1 = 0.0;
  SlopeInfo_2 = 0.0;
  
  Array<int> ind;
  mesh->GetElementVertices(i, ind);
  Array<int> edges;
  mesh->GetElementEdges(i, edges);
  for (int k = 0; k < edges.Size(); k++)
    ind.Append (mesh->GetNV() + edges[k]);

  Array<int> inda;
  Array<int> indb;
  int nx = mesh->GetNx();
  int ne = mesh->GetNE()/2;
  int nv = mesh->GetNV();
  int d = i/2;
  double a, b;

  if(i%2==0)
    {
      if(d==0)
	{
	  mesh->GetElementVertices(2, inda);
	  mesh->GetElementEdges(2, edges);
	  for (int k = 0; k <3; k++) inda.Append (nv + edges[k]);
	  a = sol(ind[1]) - sol(ind[3]);
	  b = sol(inda[3]) - sol(ind[1]);
	  SlopeInfo_0(1) = MinMod(a, b);
	  
	  mesh->GetElementVertices(2*nx, inda);
	  mesh->GetElementEdges(2*nx, edges);
	  for (int k = 0; k <3; k++) inda.Append (nv + edges[k]);
	  a = sol(ind[2]) - sol(ind[5]);
	  b = sol(inda[5]) - sol(ind[2]);
	  SlopeInfo_1(2) = MinMod(a, b);

	  a = sol(ind[3]) - sol(ind[0]);
 	  b = sol(ind[1]) - sol(ind[3]);
	  SlopeInfo_0(3) = MinMod(a, b);

	  a = sol(ind[2]) - sol(ind[5]);
	  b = sol(ind[5]) - sol(ind[0]);
	  SlopeInfo_1(5) = MinMod(a, b);

	  a = sol(ind[4]) - sol(ind[2]);
	  b = sol(ind[1]) - sol(ind[4]);
	  SlopeInfo_2(4) = MinMod(a, b);	  
	}
      else if(d==(nx-1))
	{

	}
      else if(d==(ne-nx))
	{

	}
      else if(d = (ne-1))
	{

	}
      else if(d<nx)
	{

	}
      else if(d>ne-nx && d<ne)
	{

	}
      else if(d%nx==0)
	{

	}
      else if(d%nx==nx-1)
	{

	}
      else
	{

	}
    }
  else
    {


    }
  
  
  int nx = mesh->GetNx();
  int l, r;
  l = (i%(nx+1)==0) ? i : i-1;
  r = (i%(nx+1)==nx) ? i : i+1;
  double a = sol(r) - sol(i);
  double b = sol(i) - sol(l);
  return MinMod(a, b);
  
}


 */
