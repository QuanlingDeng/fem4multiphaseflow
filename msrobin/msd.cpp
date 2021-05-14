#include "msrobin_header.h"

//===========================================================================

msd::msd(Array<int> &_n, Array<double> &_l, int _mtype, int _nbdrcon)
{
  mtype = _mtype;
  nbdrcon = _nbdrcon;
  locdata = new Data();
  locdata->SetNumberofBoundaryConditions(nbdrcon);
 
  if(_mtype==Element::TRIANGLE){locmesh = new TMesh<2>(_n, _l, Element::TRIANGLE);}
  if(_mtype==Element::QUADRILATERAL){ locmesh = new TMesh<2>(_n, _l, Element::QUADRILATERAL);}

  locpointsourcevalues.SetSize(0);
  locpointsourceindex.SetSize(0);

 totvert = locmesh->GetNV();
 dof = nbdrcon*totvert;
}

//=============================================================================
void msd::setlocellipticdata(Function *permfunc)
{
  locellipticdata.SetSize(1);
  locellipticdata[0] = new double[totvert];

  for(int k=0; k<totvert; k++)
    {
      double* coord = locmesh->GetVertex(k);
      locellipticdata[0][k] = permfunc->Eval(coord);
    }

  locdata->SetEllipticData(locellipticdata);
}

//=============================================================================
void msd::setmaterialdata(Function *modfunc, double v)
{
  locellipticdata.SetSize(2);
  locellipticdata[0] = new double[totvert];
  locellipticdata[1] = new double[totvert];

  for(int k=0; k<totvert; k++)
    {
      double* coord = locmesh->GetVertex(k);
      double data = modfunc->Eval(coord);
      locellipticdata[0][k] = data*v/((1.0+v)*(1.0-2.0*v));
      locellipticdata[1][k] = data/(1+v);
    }

   locdata->SetEllipticData(locellipticdata);
}

//=============================================================================
void msd::setmaterialdata(Array<double *> &lamecoeff)
{
  locellipticdata.SetSize(2);
  locellipticdata[0] = new double[totvert];
  locellipticdata[1] = new double[totvert];

  for(int k=0; k<totvert; k++)
    {
      locellipticdata[0][k] = lamecoeff[0][k];
      locellipticdata[1][k] = lamecoeff[1][k];
    }

   locdata->SetEllipticData(locellipticdata);
}

//=============================================================================
void msd::setlocforcedata(Function *forcefunc)
{
  locforcedata.SetSize(1);
  locforcedata[0] = new double[totvert];

  for(int k=0; k<totvert; k++)
    {
      double* coord = locmesh->GetVertex(k);
      locforcedata[0][k] = forcefunc->Eval(coord);
    }

  locdata->SetForceData(locforcedata);
}

//=============================================================================
void msd::setlocforcedata(Function *forcefuncx, Function *forcefuncy)
{
  locforcedata.SetSize(2); 
  locforcedata[0] = new double[totvert];
  locforcedata[1] = new double[totvert];

  for(int k=0; k<totvert; k++)
    {
      double* coord = locmesh->GetVertex(k);
      locforcedata[0][k] = forcefuncx->Eval(coord); 
      locforcedata[1][k] = forcefuncy->Eval(coord);
    }

  locdata->SetForceData(locforcedata);
}

//=============================================================================
void msd::setpointsourcedata(const Array<double> &_fvalue, const Array<int> &_nodeindex)
{
  locpointsourcevalues.SetSize(_fvalue.Size());
  for(int i=0; i<locpointsourcevalues.Size(); i++) { locpointsourcevalues[i] = _fvalue[i]; }

  locpointsourceindex.SetSize(_nodeindex.Size());
  for(int i=0; i<locpointsourceindex.Size(); i++) { locpointsourceindex[i] = _nodeindex[i]; }
}

//=============================================================================
void msd::setlocbndrydata(Array<double *> &_dbdryval, Array<double *> &_nbdryval,
			  Array<double *> &_rbdryval, Array<double *> &_rcoeff,
			  Array<bool> &_BdrDirichlet, Array<bool> &_BdrNeumann,
			  Array<bool> &_BdrRobin, Array<bool> &_GBdry)
{
  //build 
  GBdry.SetSize(_GBdry.Size());
  for(int b=0; b<GBdry.Size(); b++)
    {
      GBdry[b] = _GBdry[b];
    }

  BdrDirichlet.SetSize(_BdrDirichlet.Size()); 
  for(int b=0; b<BdrDirichlet.Size(); b++)
    {
      BdrDirichlet[b] = _BdrDirichlet[b];
    }

  BdrNeumann.SetSize(_BdrNeumann.Size());
  for(int b=0; b<BdrNeumann.Size(); b++)
    {
      BdrNeumann[b] = _BdrNeumann[b];
    }

  BdrRobin.SetSize(_BdrRobin.Size());
  for(int b=0; b<BdrRobin.Size(); b++)
    {
      BdrRobin[b] = _BdrRobin[b];
    }

  int nbdrs = getnbdrs();
  dbdryval.SetSize(_dbdryval.Size());
  for(int b=0; b<dbdryval.Size()/nbdrcon; b++)
    {
      for(int i=0; i<nbdrcon; i++)
	{ 
	  int bb = b + i*nbdrs;
	  int n = getnbdryelem(b)+1;
	  dbdryval[bb] = new double[n];
	  for(int k=0; k<n; k++)
	    {
	      dbdryval[bb][k] = _dbdryval[bb][k];
	    }
	}
    }
  
  nbdryval.SetSize(_nbdryval.Size());
  for(int b=0; b<nbdryval.Size()/nbdrcon; b++)
    {
      for(int i=0; i<nbdrcon; i++)
	{ 
	  int bb = b + i*nbdrs;
	  int n = getnbdryelem(b)+1;
	  nbdryval[bb] = new double[n];
	  for(int k=0; k<n; k++)
	    {
	      nbdryval[bb][k] = _nbdryval[bb][k];
	    }
	}
    }

  rbdryval.SetSize(_rbdryval.Size());
  for(int b=0; b<rbdryval.Size()/nbdrcon; b++)
    {
      for(int i=0; i<nbdrcon; i++)
	{ 
	  int bb = b + i*nbdrs;
	  int n = getnbdryelem(b)+1;
	  rbdryval[bb] = new double[n];
	  for(int k=0; k<n; k++)
	    {
	      rbdryval[bb][k] = _rbdryval[bb][k];
	    }
	}
    }

  rcoeff.SetSize(_rcoeff.Size());
  for(int b=0; b<rcoeff.Size()/nbdrcon; b++)
    {  
      for(int i=0; i<nbdrcon; i++)
	{ 
	  int bb = b + i*nbdrs;
	  int n = getnbdryelem(b)+1;
	  rcoeff[bb] = new double[n];
	  for(int k=0; k<n; k++)
	    {
	      rcoeff[bb][k] = _rcoeff[bb][k];
	    }
	}
    }
}

//=============================================================================

void msd::setgbvindex(int startindex)
{
  //global vertex index information
  gbvindex.SetSize(totvert);

  for(int k=0; k<totvert; k++)
    {
      gbvindex[k] = startindex+k;
    }
}

//=============================================================================

void msd::setnbinfo(Array<int*> &_nbinfo)
{
  //boundary information for neighboring subdomain boundaries
  nbinfo.SetSize(_nbinfo.Size());
  for(int b=0; b<nbinfo.Size(); b++)
    {
      int nb = getnbdrs();
      nbinfo[b] = new int[nbdrcon*nb];
      for(int k=0; k<nbdrcon*nb; k++)
	{
	  nbinfo[b][k] = _nbinfo[b][k];
	}
    }
}

//=============================================================================

int msd::getbdrelemgindex(int b, int p)
{
  int nbdrys = getnbdrs();
  if(b<nbdrys)
    {
      return locmesh->GetBdryGindex(b, p);
    }
  else{
    int bb = b - nbdrys;
    return totvert + locmesh->GetBdryGindex(bb, p); 
  }
}

//=============================================================================

void msd::buildmsbasisfcts(Array<int> &_ns)
{
  int nbdry = getnbdrs();
  boundary.SetSize(nbdrcon*nbdry);
  int index = 0;

  for(int i=0; i<nbdrcon; i++)
    {
      for(int b=0; b<nbdry; b++)
	{
	  int bb = b + i*nbdry;
	  int ne = getnbdryelem(b);

	  if(GBdry[bb] == false)
	    {
	      boundary[bb] = new msbdry(_ns[bb], ne, index);
	      index += _ns[bb]+1;
	    }
	  else
	    {
	      boundary[bb] = new msbdry(ne, index); 
 	      index++;
	    }

	}
    }
  tbf = index+1;

  //---------------------------------------------------------
  // Compute Basis Functions
  //--------------------------------------------------------- 
  char *linalgsolver = new char[256];
  strcpy(linalgsolver, "multifrontal");

  basisfcts.SetSize(tbf);
  for (int k=0; k<tbf; k++) 
    {
      basisfcts[k] = new Vector(nbdrcon*totvert);
      *basisfcts[k] = 0.0;
    }

  Array<double *> tempdbdryval(dbdryval.Size());
  Array<double *> tempnbdryval(nbdryval.Size());
  Array<double *> temprbdryval(rbdryval.Size());

  for(int i=0; i<nbdrcon; i++)
    {
      for(int b=0; b<nbdry; b++)
	{
	  int bb = b + i*nbdry;
	  int n = getnbdryelem(b)+1;
	  tempdbdryval[bb] = new double[n];
	  tempnbdryval[bb] = new double[n];
	  temprbdryval[bb] = new double[n];
	  for(int k=0; k<n; k++)
	    {
	      tempdbdryval[bb][k] = 0.0;
	      tempnbdryval[bb][k] = 0.0;
	      temprbdryval[bb][k] = 0.0;
	    }
	}
    }  
  locdata->SetNeumannData(tempnbdryval, BdrNeumann);
  locdata->SetDirichletData(tempdbdryval, BdrDirichlet);
  locdata->SetRobinData(rcoeff, temprbdryval, BdrRobin);

  Array<double *> templocforcedata(nbdrcon);
  for(int k=0; k<nbdrcon; k++){ templocforcedata[k] = new double[totvert]; }
  for(int i=0; i<nbdrcon; i++)
    {
      for(int k=0; k<totvert; k++)
	{
	  templocforcedata[i][k] = 0.0;
	}
    }
  locdata->SetForceData(templocforcedata);

  for(int i=0; i<nbdrcon; i++)
    {
      for(int b=0; b<nbdry; b++)
	{
	  int bb = b + i*nbdry;
	  index = boundary[bb]->getlocindex();
	  if(GBdry[bb]==false)
	    {
	      int numbf = boundary[bb]->getnumseg()+1;
	      int n = getnbdryelem(b)+1;
	      
	      for(int k=0; k<numbf; k++)
		{
		  
		  for(int p=0; p<n; p++){ temprbdryval[bb][p] = boundary[bb]->getphifctval(k, p); }		  
		  locdata->SetRobinData(rcoeff, temprbdryval, BdrRobin);
		  
		  FEM *tempfem;
		  if(mtype==Element::TRIANGLE && nbdrcon==1){ tempfem = new LinearFEM2D(locmesh, locdata); } 	  
		  if(mtype==Element::QUADRILATERAL && nbdrcon==1){ tempfem = new BilinearFEM2D(locmesh, locdata); }
		  if(mtype==Element::TRIANGLE && nbdrcon==2){ tempfem = new LinearFEMELAST2D(locmesh, locdata); }
		  if(mtype==Element::QUADRILATERAL && nbdrcon==2){ tempfem = new BilinearFEMELAST2D(locmesh, locdata); }
		  
		  tempfem->Assemble();
		  tempfem->FinalizeMatrix();
		  tempfem->Solve(*basisfcts[index], linalgsolver);

		  
		  delete tempfem;
		  index++;
	      
		}
	      for(int p=0; p<n; p++){ temprbdryval[bb][p] = 0.0; }
	      locdata->SetRobinData(rcoeff, temprbdryval, BdrRobin);
	    }
	  
	  if(GBdry[bb]==true)
	    {
	      int n = getnbdryelem(b)+1;
	      for(int p=0; p<n; p++)
		{
		  tempdbdryval[bb][p] = dbdryval[bb][p];
		  tempnbdryval[bb][p] = nbdryval[bb][p];
		  temprbdryval[bb][p] = rbdryval[bb][p];
		 }
	      
	      locdata->SetNeumannData(tempnbdryval, BdrNeumann);
	      locdata->SetDirichletData(tempdbdryval, BdrDirichlet);
	      locdata->SetRobinData(rcoeff, temprbdryval, BdrRobin);
	      
	      FEM *tempfem;
	      if(mtype==Element::TRIANGLE && nbdrcon==1){ tempfem = new LinearFEM2D(locmesh, locdata); } 	  
	      if(mtype==Element::QUADRILATERAL && nbdrcon==1){ tempfem = new BilinearFEM2D(locmesh, locdata); }	      
	      if(mtype==Element::TRIANGLE && nbdrcon==2){ tempfem = new LinearFEMELAST2D(locmesh, locdata); }
	      if(mtype==Element::QUADRILATERAL && nbdrcon==2){ tempfem = new BilinearFEMELAST2D(locmesh, locdata); }
	      
	      tempfem->Assemble();
	      tempfem->FinalizeMatrix();
	      tempfem->Solve(*basisfcts[index], linalgsolver);
	       
	      delete tempfem;	     
	      for(int p=0; p<n; p++)
		{
		  tempdbdryval[bb][p] = 0.0;
		  tempnbdryval[bb][p] = 0.0;
		  temprbdryval[bb][p] = 0.0;
		}
	      
	      locdata->SetNeumannData(tempnbdryval, BdrNeumann);
	      locdata->SetDirichletData(tempdbdryval, BdrDirichlet);
	      locdata->SetRobinData(rcoeff, temprbdryval, BdrRobin);
	    } 
	}
    }


  locdata->SetForceData(locforcedata);
  locdata->SetNeumannData(tempnbdryval, BdrNeumann);
  locdata->SetDirichletData(tempdbdryval, BdrDirichlet);
  locdata->SetRobinData(rcoeff, temprbdryval, BdrRobin);
   
  if(locpointsourceindex.Size()>0)
    {
      locdata->SetPointSourceData(locpointsourcevalues, locpointsourceindex);	
    }
  
  FEM *tempfem;	 
  if(mtype==Element::TRIANGLE && nbdrcon==1){ tempfem = new LinearFEM2D(locmesh, locdata); } 	  
  if(mtype==Element::QUADRILATERAL && nbdrcon==1){ tempfem = new BilinearFEM2D(locmesh, locdata); }
  if(mtype==Element::TRIANGLE && nbdrcon==2){ tempfem = new LinearFEMELAST2D(locmesh, locdata); }
  if(mtype==Element::QUADRILATERAL && nbdrcon==2){ tempfem = new BilinearFEMELAST2D(locmesh, locdata); }
  tempfem->Assemble();
  tempfem->FinalizeMatrix();
  tempfem->Solve(*basisfcts[tbf-1], linalgsolver);
 
  delete tempfem;
  for(int k=0; k<nbdryval.Size(); k++)
    delete []tempnbdryval[k];

  for(int k=0; k<dbdryval.Size(); k++)
    delete []tempdbdryval[k];

  for(int k=0; k<rbdryval.Size(); k++)
    delete []temprbdryval[k];

  for(int k=0; k<locforcedata.Size(); k++)
    delete []templocforcedata[k];

  delete []linalgsolver;

}

//=============================================================================

void msd::Build(Vector &sol, double* alpha)
{
  sol.SetSize(dof);
  sol = 0.0;

  for (int k=0; k<tbf; k++)
    {
      for (int l=0; l<dof; l++)
	{
	  sol(l) += alpha[k]*basisfcts[k]->Elem(l); 
	}
    }
}

//=============================================================================

double msd::getbdryelemlength(int b, int p)
{
  Array<double *> coord(2);
  coord[0] = new double[2];
  coord[1] = new double[2];
  
  locmesh->GetBdrElementVerticesCoord(p, coord);

  double dx, dy, leng;
  dx = coord[0][0] - coord[0][1];
  dy = coord[1][0] - coord[1][1];

  leng = sqrt(pow(dx,2) + pow(dy,2));

  delete []coord[0];
  delete []coord[1];

  return leng;
}

//=============================================================================

msd::~msd()
{
  delete locmesh;
  delete locdata;

 for(int k=0; k<locellipticdata.Size(); k++)
    delete []locellipticdata[k];

 for(int k=0; k<locforcedata.Size(); k++)
    delete []locforcedata[k];

  for(int k=0; k<nbdryval.Size(); k++)
    delete []nbdryval[k];

  for(int k=0; k<dbdryval.Size(); k++)
    delete []dbdryval[k];

  for(int k=0; k<rbdryval.Size(); k++)
    delete []rbdryval[k];

  for(int k=0; k<rcoeff.Size(); k++)
    delete []rcoeff[k];

  for(int k=0; k<boundary.Size(); k++)
    delete boundary[k];

  for(int k=0; k<basisfcts.Size(); k++)
    delete basisfcts[k];

  for(int k=0; k<nbinfo.Size(); k++)
    delete []nbinfo[k];




}

