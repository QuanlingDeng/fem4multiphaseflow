#include "msrobin_header.h"

//===========================================================================

directMsRM::directMsRM(Array<msd*> &_subdomains, int &_dim, int _numbdrcon)
{
  //set dimension of global system
  dim = _dim;
  nbdrcon = _numbdrcon;
  Nsd = _subdomains.Size();

  //Gather Subdomain Data
  subdomains.SetSize(Nsd);
  for(int k=0; k<Nsd; k++) { subdomains[k] = _subdomains[k]; }

  //Set memory for solution storage
  gsol.SetSize(Nsd);
  alpha.SetSize(Nsd);

  //Build Solution Array
  for(int k=0; k<Nsd; k++)
    {
      int tbf = subdomains[k]->getnumbasis();
      int totvert = subdomains[k]->gettotvert();   
      gsol[k] = new Vector(nbdrcon*subdomains[k]->gettotvert()); 
      alpha[k] = alpha[k] = new double[tbf];
      for (int l=0; l<tbf; l++){ alpha[k][l] = 1.0; }
    }
}

//====================================================================

void directMsRM::solve()
{
  cout << "Noniterative Multiscale Robin Method Start" << endl;
  mat = new SparseMatrix(dim, dim);
  solalpha.SetSize(dim);
  rhs.SetSize(dim);
  rhs = 0.0;
  solalpha = 1.0;

 //fill the matrix
  cout << "filling matrix...";
  for(int k=0; k<Nsd; k++)
    {
      int nbdry = subdomains[k]->getnbdrs();
      for(int i=0; i<nbdrcon; i++)
	{
	  for(int b=0; b<nbdry; b++)
	    {
	      int bc = b + i*nbdry;
	      if(subdomains[k]->getGBdry(bc) == false)
		{
		  addlocalsystem(b, bc, k);
		}
	    }
	}
    }

  mat->Finalize();

  ofstream out;
  out.open("msmatrix.out");
  mat->MPrint(out);
  out.close();
  //mat->Print();
  //rhs.Print();

  cout << "done" << endl << "solving...";

  //solve
  char *linalgsolver = new char[256];
  strcpy(linalgsolver, "multifrontal");
  solve(solalpha, linalgsolver, 2*dim, 0, 1e-22, 1e-20); 
  delete []linalgsolver;

  cout << "done" << endl;

  //solalpha.Print();
  //set alphas
  for(int k=0; k<Nsd; k++)
    {
      int nbdry = subdomains[k]->getnbdrs();
      for(int i=0; i<nbdrcon; i++)
	{
	  for(int b=0; b<nbdry; b++)
	    {
	      int bc = b + i*nbdry;
	      int locind = subdomains[k]->getlocindex(bc);
	      if(subdomains[k]->getGBdry(bc) == false)
		{
		  int mns = subdomains[k]->getnumseg(bc);
		  int gbind = subdomains[k]->getgbindex(bc);
		  
		  for(int v=0; v<=mns; v++)
		    {
		      alpha[k][locind+v] = solalpha(gbind+v);
		    }
		}
	      else{ alpha[k][locind] = 1.0; }
	    }
	}
	  
      subdomains[k]->Build(*gsol[k], alpha[k]); //bulid subdomain solution
    }

  cout << "Noniterative Multiscale Robin Method End" << endl;
 
}
//=============================================================================

void directMsRM::addlocalsystem(int b, int bc, int k)
{
  int kn = subdomains[k]->getnbinfo(0, b);
  int bn = subdomains[k]->getnbinfo(1, b);
  int bcn = subdomains[k]->getnbinfo(1, bc);

  //construct the local system governing the interface condition
  int mns = subdomains[k]->getnumseg(bc);
  int nns = subdomains[kn]->getnumseg(bcn);
  int mind = subdomains[k]->getgbindex(bc); 
  int nind = subdomains[kn]->getgbindex(bcn); 
  int ne = subdomains[k]->getnbdryelem(b);

  //double gamma = subdomains[kn]->getgamma(bn);
  double gamma = subdomains[kn]->getrobincoeff(bn);

  int bdryelemindex = 0;
  for(int bb=0; bb<b; bb++)
    {
      bdryelemindex += subdomains[k]->getnbdryelem(bb);
    }

  int nbdryelemindex = 0;
  for(int bb=0; bb<bn; bb++)
    {
      nbdryelemindex += subdomains[kn]->getnbdryelem(bb);
    }

  Array<double *> bdryfcts(2);
  for(int l=0; l<2; l++){ bdryfcts[l] = new double[ne+1]; }

  Array<double> h(ne);
  for(int p=0; p<ne; p++){h[p] = subdomains[k]->getbdryelemlength(b, bdryelemindex+p); } 

  for(int v=0; v<=mns; v++)
    {
      for(int p=0; p<=ne; p++){bdryfcts[0][p] = subdomains[k]->getbdryfctval(bc, v, p); }

      //Integrate hat/phi functions against hat functions and add to matrix
      for(int u=0; u<=mns; u++)
	{
	  for(int p=0; p<=ne; p++){ bdryfcts[1][p] = subdomains[k]->getbdryfctval(bc, u, p); }
	  mat->Elem(mind+v, mind+u) = intbdryfcts(bdryfcts, ne, h);
	}
	      
      for(int u=0; u<=nns; u++)
	{
	  for(int p=0; p<=ne; p++){ bdryfcts[1][p] = subdomains[kn]->getbdryfctval(bcn, u, p); } 
	  mat->Elem(mind+v, nind+u) = intbdryfcts(bdryfcts, ne, h);
	}

      
      //Integrate Hat/phi function against basis functions and add to matrix
      int nbdry = subdomains[kn]->getnbdrs();
      for(int i=0; i<nbdrcon; i++)
	{
	  for(int bb=0; bb<nbdry; bb++)
	    {  
	      int bbc = bb + i*nbdry;
	      int locind = subdomains[kn]->getlocindex(bbc);
	      
	      if(subdomains[kn]->getGBdry(bbc) == false)
		{
		  int gbind = subdomains[kn]->getgbindex(bbc);
		  int nb = subdomains[kn]->getnumseg(bbc);

		  for(int kk=0; kk<=nb; kk++)
		    {
		      for(int p=0; p<=ne; p++)
			{
			  int index = subdomains[kn]->getbdrelemgindex(bcn, p);
			  bdryfcts[1][p] = subdomains[kn]->getbasisval(locind+kk, index);
			} 
		      mat->Elem(mind+v, gbind+kk) += 2.0*gamma*intbdryfcts(bdryfcts, ne, h);		  
		    }
		}
	      else
		{
		  for(int p=0; p<=ne; p++)
		    {
		      int index = subdomains[kn]->getbdrelemgindex(bcn, p);
		      bdryfcts[1][p] = subdomains[kn]->getbasisval(locind, index); 
		    } 
		  rhs(mind+v) -= 2.0*gamma*intbdryfcts(bdryfcts, ne, h);
		} 
	      
	    }
	}
      
      int hatindex = subdomains[kn]->getnumbasis()-1;
      for(int p=0; p<=ne; p++)
	{
	  int index = subdomains[kn]->getbdrelemgindex(bcn, p);
	  bdryfcts[1][p] = subdomains[kn]->getbasisval(hatindex, index);
	} 
      rhs(mind+v) -=  2.0*gamma*intbdryfcts(bdryfcts, ne, h);
    }	 
	  
  for(int l=0; l<2; l++){delete []bdryfcts[l]; }
}

//=============================================================================

double directMsRM::intbdryfcts(Array<double *> &bdryfcts, int n, Array<double> &h)
{
  //integrate two boundary fcts that are peicwise linear
  double val = 0.0;
  for(int p=0; p<n; p++)
    {
      double psi0 = bdryfcts[0][p];
      double phi0 = bdryfcts[1][p];
      
      double psi1 = bdryfcts[0][p+1];
      double phi1 = bdryfcts[1][p+1];
      
      val += h[p]*((psi0*phi0 + psi1*phi1)/3.0 + (psi0*phi1 + psi1*phi0)/6.0);
    }  
  return val;
}

//=============================================================================

void directMsRM::buildsol(Vector &sol, int Nx, int Ny, int nx, int ny)
{
  int totGvert = (Nx*nx+1)*(Ny*ny+1);
  int dof = nbdrcon*totGvert;
  sol.SetSize(dof);

  Array<int> count(dof);
  for(int j=0; j<dof; j++)
    {
      count[j]=0;
      sol(j) = 0.0;  
    }
 
  for(int u=0; u<nbdrcon; u++)
    {
      int offset = u*totGvert;
      for (int j=0; j<Ny; j++)
	{
	  for (int i=0; i<Nx; i++)
	    {
	      int Gind = j*Nx + i;
	      int totvert = subdomains[Gind]->gettotvert();
	      for (int k=0; k<=ny; k++)
		{
		  for (int l=0; l<=nx; l++)
		    {
		      int FGind = offset + l + i*nx + (k + j*ny)*(Nx*nx + 1);
		      int ind = u*totvert + k*(nx+1) + l;
		      sol(FGind) += gsol[Gind]->Elem(ind);
		      count[FGind]++;
		    }
		}
	      
	    }
	}
    }
  for(int j=0; j<dof; j++)
    {
      sol(j) /= ((double) count[j]);
    }
  
}

//=============================================================================

void directMsRM::solve(Vector &sol, char *linalgsolv, int maxit, int printit, 
		double rtol, double atol)
{

  sol = 0.0;
  if (!strcmp(linalgsolv, "multifrontal")) 
    {
      int mem_fact = 100;
      MultiFrontalMatrixInverse *InvMat = 
	new MultiFrontalMatrixInverse(*mat, mem_fact);
      InvMat->Mult(rhs, sol);
      delete InvMat;
    }
  else 
    {      
      MatrixInverse *prec = new GSSmoother(*mat);
      MatrixInverse *invmat;
      
      if (!strcmp(linalgsolv, "pcg"))
	invmat = new PCGMatrixInverse(*mat, *prec, printit, maxit, rtol, atol);
      if (!strcmp(linalgsolv, "bicgstab"))
	invmat = new BICGSTABMatrixInverse(*mat, *prec, printit, maxit, 
					   rtol, atol);
      if (!strcmp(linalgsolv, "gmres"))
	invmat = new GMRESMatrixInverse(*mat, *prec, printit, maxit, 40, 
					rtol, atol);
      
      invmat->Mult(rhs, sol);
      delete invmat;
      delete prec;
    }
}

//=============================================================================

double directMsRM::computel2error(Function *exactsol)
{  
  double l2error = 0.0;
  for(int k=0; k<Nsd; k++)
    {
      int mtype = subdomains[k]->getmeshtype();
      Mesh *tempmesh = subdomains[k]->getmesh();
      Data *tempdata = subdomains[k]->getdata();
      FEM *tempfem;
      if(mtype==Element::QUADRILATERAL){ tempfem = new BilinearFEM2D(tempmesh, tempdata); }
      if(mtype==Element::TRIANGLE){ tempfem = new LinearFEM2D(tempmesh, tempdata); }

      double temperr = tempfem->ComputeL2Error(exactsol, *gsol[k]);
      l2error += temperr*temperr;

      delete tempfem;
    }
  
  return sqrt(l2error);
}

//=============================================================================

double directMsRM::computel2error(Function *exactsolx, Function *exactsoly)
{  
  double l2errorx = 0.0;
  double l2errory = 0.0;
  for(int k=0; k<Nsd; k++)
    {
      int mtype = subdomains[k]->getmeshtype();
      int meshnv = subdomains[k]->gettotvert();
      Mesh *tempmesh = subdomains[k]->getmesh();
      Data *tempdata = subdomains[k]->getdata();
      FEM *tempfem;
      if(mtype==Element::QUADRILATERAL){ tempfem = new BilinearFEM2D(tempmesh, tempdata); }
      if(mtype==Element::TRIANGLE){ tempfem = new LinearFEM2D(tempmesh, tempdata); }

      Vector tempu1(meshnv);
      Vector tempu2(meshnv);
      for(int i=0; i<meshnv; i++){tempu1(i) = gsol[k]->Elem(i); tempu2(i) = gsol[k]->Elem(i+meshnv);}
      double temperrx = tempfem->ComputeL2Error(exactsolx, tempu1);
      double temperry = tempfem->ComputeL2Error(exactsoly, tempu2);
      l2errorx += temperrx*temperrx;
      l2errory += temperry*temperry;

      delete tempfem;
    }
  
  cout << "L2error u1: " << sqrt(l2errorx) << endl;
  cout << "L2error u2: " << sqrt(l2errory) << endl << endl;

  return sqrt(l2errorx+l2errory);
}

//=============================================================================

double directMsRM::computeh1error(Array<Function *> &exactderivative)
{  
  double h1error = 0.0;
  for(int k=0; k<Nsd; k++)
    {
      int mtype = subdomains[k]->getmeshtype();
      Mesh *tempmesh = subdomains[k]->getmesh();
      Data *tempdata = subdomains[k]->getdata();
      FEM *tempfem;
      if(mtype==Element::QUADRILATERAL){ tempfem = new BilinearFEM2D(tempmesh, tempdata); }
      if(mtype==Element::TRIANGLE){ tempfem = new LinearFEM2D(tempmesh, tempdata); }

      double temperr = tempfem->ComputeH1Error(exactderivative, *gsol[k]);
      h1error += temperr*temperr;

      delete tempfem;
    }
  
  return sqrt(h1error);
}

//=============================================================================

void directMsRM::paraviewprintmesh()
{
  int nv = 0;
  int ne = 0;
  int ncdp = 0;
  for(int k=0; k<Nsd; k++)
    {
      int numel = subdomains[k]->getne();
      nv += subdomains[k]->gettotvert();
      ne += numel;
      for(int i=0; i<numel; i++)
	{
	  ncdp += subdomains[k]->getelementnvertices(i) + 1;
	}
    }

  ofstream out;
  out.open("msmesh.vtk");
  out << "# vtk DataFile Version 2.0" << endl;
  out << "Unstructured Grid" << endl;
  out << "ASCII" << endl;
  out << "DATASET UNSTRUCTURED_GRID" << endl;
  out << "POINTS " << nv << " DOUBLE" << endl;

  for(int k=0; k<Nsd; k++)
    {
      int numv = subdomains[k]->gettotvert();
      for(int j=0; j<numv; j++)
	{
	  for(int i=0; i<2; i++)
	    {
	      double coord = subdomains[k]->getvertex(j, i);
	      out << coord << " ";
	    }
	  out << 0.0 << endl;
	}
    }

  
  out << endl;
  out << "CELLS " << ne << " " << ncdp << endl;

  for(int k=0; k<Nsd; k++)
    {
      int numel = subdomains[k]->getne();
      for(int j=0; j<numel; j++)
	{
	  int numv = subdomains[k]->getelementnvertices(j);
	  out << numv << "     ";
	  Array<int> ind;
	  subdomains[k]->getelementvertices(j, ind);
	  for(int i=0; i<numv; i++)
	    {
	      out << subdomains[k]->getgbvindex(ind[i]) << "    ";
	    }
	  out << endl;
	}
    }
  
  out << endl;
  out << "CELL_TYPES " << ne << endl;
  for(int k=0; k<Nsd; k++)
    {
      int type = subdomains[k]->getmeshtype();
      int ctype;
      if(type==Element::QUADRILATERAL){ctype=9;  }
      if(type==Element::TRIANGLE){ ctype=5; }

      int numel = subdomains[k]->getne();
      for(int j=0; j<numel; j++)
	{
	  out << ctype << endl;
	}
    }
      
  out.close();
  
}

//=============================================================================

//=============================================================================

void directMsRM::paraviewprintdeformedmesh()
{
  int nv = 0;
  int ne = 0;
  int ncdp = 0;
  for(int k=0; k<Nsd; k++)
    {
      int numel = subdomains[k]->getne();
      nv += subdomains[k]->gettotvert();
      ne += numel;
      for(int i=0; i<numel; i++)
	{
	  ncdp += subdomains[k]->getelementnvertices(i) + 1;
	}
    }

  ofstream out;
  out.open("msdefmesh.vtk");
  out << "# vtk DataFile Version 2.0" << endl;
  out << "Unstructured Grid" << endl;
  out << "ASCII" << endl;
  out << "DATASET UNSTRUCTURED_GRID" << endl;
  out << "POINTS " << nv << " DOUBLE" << endl;

  for(int k=0; k<Nsd; k++)
    {
      int numv = subdomains[k]->gettotvert();
      for(int j=0; j<numv; j++)
	{
	  for(int i=0; i<2; i++)
	    {
	      double coord = subdomains[k]->getvertex(j, i) + gsol[k]->Elem(j+i*numv) ;
	      out << coord << " ";
	    }
	  out << 0.0 << endl;
	}
    }

  
  out << endl;
  out << "CELLS " << ne << " " << ncdp << endl;

  for(int k=0; k<Nsd; k++)
    {
      int numel = subdomains[k]->getne();
      for(int j=0; j<numel; j++)
	{
	  int numv = subdomains[k]->getelementnvertices(j);
	  out << numv << "     ";
	  Array<int> ind;
	  subdomains[k]->getelementvertices(j, ind);
	  for(int i=0; i<numv; i++)
	    {
	      out << subdomains[k]->getgbvindex(ind[i]) << "    ";
	    }
	  out << endl;
	}
    }
  
  out << endl;
  out << "CELL_TYPES " << ne << endl;
  for(int k=0; k<Nsd; k++)
    {
      int type = subdomains[k]->getmeshtype();
      int ctype;
      if(type==Element::QUADRILATERAL){ctype=9;  }
      if(type==Element::TRIANGLE){ ctype=5; }

      int numel = subdomains[k]->getne();
      for(int j=0; j<numel; j++)
	{
	  out << ctype << endl;
	}
    }
      
  out.close();
  
}

//=============================================================================

void directMsRM::paraviewprintsol() 
{
  int nv = 0;
  int ne = 0;
  int ncdp = 0;
  for(int k=0; k<Nsd; k++)
    {
      int numel = subdomains[k]->getne();
      nv += subdomains[k]->gettotvert();
      ne += numel;
      for(int i=0; i<numel; i++)
	{
	  ncdp += subdomains[k]->getelementnvertices(i) + 1;
	}
    }

  ofstream out;
  out.open("mssol.vtk");
  out << "# vtk DataFile Version 2.0" << endl;
  out << "Unstructured Grid" << endl;
  out << "ASCII" << endl;
  out << "DATASET UNSTRUCTURED_GRID" << endl;
  out << "POINTS " << nv << " DOUBLE" << endl;

  for(int k=0; k<Nsd; k++)
    {
      int numv = subdomains[k]->gettotvert();
      for(int j=0; j<numv; j++)
	{
	  for(int i=0; i<2; i++)
	    {
	      double coord = subdomains[k]->getvertex(j, i);
	      out << coord << " ";
	    }
	  out << 0.0 << endl;
	}
    }

  
  out << endl;
  out << "CELLS " << ne << " " << ncdp << endl;

  for(int k=0; k<Nsd; k++)
    {
      int numel = subdomains[k]->getne();
      for(int j=0; j<numel; j++)
	{
	  int numv = subdomains[k]->getelementnvertices(j);
	  out << numv << "     ";
	  Array<int> ind;
	  subdomains[k]->getelementvertices(j, ind);
	  for(int i=0; i<numv; i++)
	    {
	      out << subdomains[k]->getgbvindex(ind[i]) << "    ";
	    }
	  out << endl;
	}
    }
  
  out << endl;
  out << "CELL_TYPES " << ne << endl;
  for(int k=0; k<Nsd; k++)
    {
      int type = subdomains[k]->getmeshtype();
      int ctype;
      if(type==Element::QUADRILATERAL){ctype=9;  }
      if(type==Element::TRIANGLE){ ctype=5; }

      int numel = subdomains[k]->getne();
      for(int j=0; j<numel; j++)
	{
	  out << ctype << endl;
	}
    }

  
  out << endl;
  out << "POINT_DATA " << nv << endl;
  out << "SCALARS pval DOUBLE" << endl;
  out << "LOOKUP_TABLE default" << endl;

 
  for(int k=0; k<Nsd; k++)
    {
      int numv = subdomains[k]->gettotvert();
      for(int j=0; j<numv; j++)
	{
	  out << gsol[k]->Elem(j) << endl;
	}	
    }
      
  out.close();
  
}

//=============================================================================

void directMsRM::paraviewprintellipticdata() 
{
  int nv = 0;
  int ne = 0;
  int ncdp = 0;
  for(int k=0; k<Nsd; k++)
    {
      int numel = subdomains[k]->getne();
      nv += subdomains[k]->gettotvert();
      ne += numel;
      for(int i=0; i<numel; i++)
	{
	  ncdp += subdomains[k]->getelementnvertices(i) + 1;
	}
    }
  for(int b=0; b<nbdrcon; b++)
    {
      char fname[256];
      sprintf(fname, "msellpticdata_%d.vtk", b);

      ofstream out;
      out.open(fname);
      out << "# vtk DataFile Version 2.0" << endl;
      out << "Unstructured Grid" << endl;
      out << "ASCII" << endl;
      out << "DATASET UNSTRUCTURED_GRID" << endl;
      out << "POINTS " << nv << " DOUBLE" << endl;

      for(int k=0; k<Nsd; k++)
	{
	  int numv = subdomains[k]->gettotvert();
	  for(int j=0; j<numv; j++)
	    {
	      for(int i=0; i<2; i++)
		{
		  double coord = subdomains[k]->getvertex(j, i);
		  out << coord << " ";
		}
	      out << 0.0 << endl;
	    }
	}

  
      out << endl;
      out << "CELLS " << ne << " " << ncdp << endl;

      for(int k=0; k<Nsd; k++)
	{
	  int numel = subdomains[k]->getne();
	  for(int j=0; j<numel; j++)
	    {
	      int numv = subdomains[k]->getelementnvertices(j);
	      out << numv << "     ";
	      Array<int> ind;
	      subdomains[k]->getelementvertices(j, ind);
	      for(int i=0; i<numv; i++)
		{
		  out << subdomains[k]->getgbvindex(ind[i]) << "    ";
		}
	      out << endl;
	    }
	}
  
      out << endl;
      out << "CELL_TYPES " << ne << endl;
      for(int k=0; k<Nsd; k++)
	{
	  int type = subdomains[k]->getmeshtype();
	  int ctype;
	  if(type==Element::QUADRILATERAL){ctype=9;  }
	  if(type==Element::TRIANGLE){ ctype=5; }

	  int numel = subdomains[k]->getne();
	  for(int j=0; j<numel; j++)
	    {
	      out << ctype << endl;
	    }
	}

  
      out << endl;
      out << "POINT_DATA " << nv << endl;
      out << "SCALARS pval DOUBLE" << endl;
      out << "LOOKUP_TABLE default" << endl;

 
      for(int k=0; k<Nsd; k++)
	{
	  int numv = subdomains[k]->gettotvert();
	  for(int j=0; j<numv; j++)
	    {
	      out << subdomains[k]->getellipticdata(j, b) << endl;
	    }	
	}
      
      out.close();
    }
}


//=============================================================================
directMsRM::~directMsRM()
{
  delete mat;
  for (int i=0; i<gsol.Size(); i++)
    {
      delete gsol[i];
    }
  for (int i=0; i<alpha.Size(); i++)
    {
      delete []alpha[i];
    }
}
