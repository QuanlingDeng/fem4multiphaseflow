#include "fem_header.h"

BasisFunctions::BasisFunctions(int _order)
{
  if(_order>3)
    {
      cout << "Basis Functions are not constructed for polynomial with order " << _order << endl;
      exit(1);
    }
  
  order = _order; 
}


//=========================================================

double BasisFunctions::BF1D(int ind, double x)
{
  //if(x>1 || x<0) cout << "Input Error:" << x << endl;
  
  if(order==0)
    return 1.0;
  else if(order==1)
    {
      if(ind==0)
	return 1.0 - x;
      else
	return x;
    }
  else if(order==2)
    {
      if(ind==0)
	return (1.0 - x) * (1.0 - 2.0*x);
      else if(ind==1)
	return x * (2.0*x - 1.0);
      else
	return 4.0 * x * (1.0 - x);
    }
  else
    {
      if(ind==0)
	return (1.0 - x) * (1.0 - 3.0*x) * (1.0 - 1.5*x);
      else if(ind==1)
	return x * (1.0 - 3.0*x) * (1.0 - 1.5*x);
      else if(ind==2)
	return 9.0 * x * (1.0 - x) * (1.0 - 1.5*x);
      else
	return -2.25 * x * (1.0 - x) * (1.0 - 3.0*x) ;
    }
}

double BasisFunctions::DerBF1D(int ind, double x)
{
  //if(x>1 || x<0) cout << "Input Error:" << x << endl;
  
  if(order==0)
    return 0.0;
  else if(order==1)
    {
      if(ind==0)
	return -1.0;
      else
	return 1.0;
    }
  else if(order==2)
    {
      if(ind==0)
	return 4.0 * x - 3.0; 
      else if(ind==1)
	return 4.0 * x - 1.0; 
      else
	return 4.0 - 8.0 * x;
    }
  else
    {
      if(ind==0)
	return -13.5 * x * x + 18.0 *x - 5.5;
      else if(ind==1)
	return 13.5 * x * x - 9.0 *x + 1.0;
      else if(ind==2)
	return 40.5 * x * x - 45.0 * x - 9.0;
      else
	return -20.25 * x * x + 18.0 * x - 2.25; 
    }
}

//=========================================================

double BasisFunctions::BF2D(int ind, double x, double y)
{
  //if(x>1 || x<0 || y>1 || y<0) cout << "Input Error: x = " << x << " and y = " << y << endl;
  
  if(order==0)
    return 1.0;
  else if(order==1)
    {
      if(ind==0)
	return 1.0 - x - y;
      else if(ind==1)
	return x;
      else
	return y;
    }
  else if(order==2)
    {
      if(ind==0)
	return (1.0 - x -y) * (1.0 - 2.0*x - 2.0*y);
      else if(ind==1)
	return x * (2.0*x - 1.0);
      else if(ind==2)
	return y * (2.0*y - 1.0);
      else if(ind==3)
	return 4.0 * x * (1.0 - x -y);
      else if(ind==4)
	return 4.0 * x * y;
      else
	return 4.0 * y * (1.0 - x - y);
    }
  else
    {
      if(ind==0)
	return (1.0 - x - y) * (1.0 - 3.0*x - 3.0*y) * (1.0 - 1.5*x - 1.5*y);
      else if(ind==1)
	return x * (3.0*x - 1.0) * (1.5*x - 1.0);
      else if(ind==2)
	return y * (3.0*y - 1.0) * (1.5*y - 1.0);
      else if(ind==3)
	return 9.0 * x * (1.0 - x - y) * (1.0 - 1.5*x - 1.5*y);
      else if(ind==4)
	return 4.5 * x * (1.0 - x - y) * (3.0*x - 1.0);
      else if(ind==5)
	return 4.5 * x * y * (3.0 * x - 1.0);
      else if(ind==6)
	return 4.5 * x * y * (3.0 * y - 1.0);
      else if(ind==7)
	return 4.5 * y * (1.0 - x - y) * (3.0*y - 1.0);
      else if(ind==8)
	return 9.0 * y * (1.0 - x - y) * (1.0 - 1.5*x - 1.5*y);
      else
	return 27.0 * x * y * (1.0 - x - y);
    }
}

double BasisFunctions::GradBF2D(int ind, int dimflag, double x, double y)
{
  //if(x>1 || x<0 || y>1 || y<0) cout << "Input Error: x = " << x << " and y = " << y << endl;
  
  if(order==0) return 0.0;
  
  if(dimflag==0)
    {
      if(order==1)
	{
	  if(ind==0)
	    return -1.0;
	  else if(ind==1)
	    return 1.0;
	  else
	    return 0.0;
	}
      else if(order==2)
	{ 
	  if(ind==0)
	    return 4.0*x + 4.0*y - 3.0;
	  else if(ind==1)
	    return 4.0*x - 1.0;
	  else if(ind==2)
	    return 0.0;
	  else if(ind==3)
	    return 4.0 * ( 1.0 - y - 2.0*x );
	  else if(ind==4)
	    return 4.0 * y;
	  else
	    return -4.0 * y;
	}
      else
	{
	  if(ind==0)
	    return 0.5 * ( -27.0*x*x - 27.0*y*y - 54.0*x*y + 36.0*x + 36.0*y - 11.0 );
	  else if(ind==1)
	    return 1.0 - 9.0*x + 13.5*x*x;
	  else if(ind==2)
	    return 0.0;
	  else if(ind==3)
	    return 4.5 * ( 2.0 + 9.0*x*x - 5.0*y + 3.0*y*y - 10.0*x + 12.0*x*y );
	  else if(ind==4)
	    return -4.5 * ( 1.0 + 9.0*x*x - y - 8.0*x + 6.0*x*y);
	  else if(ind==5)
	    return 4.5 * y * ( -1.0 + 6.0*x );
	  else if(ind==6)
	    return 4.5 * y * ( -1.0 + 3.0*y );
	  else if(ind==7)
	    return -4.5 * y * ( -1.0 + 3.0*y );
	  else if(ind==8)
	    return 4.5 * y * ( -5.0 + 6.0*x + 6.0*y );
	  else
	    return 27.0 * y * ( 1.0 - 2.0*x - y ); 
	}
    }
  else
    {
      if(order==1)
	{
	  if(ind==0)
	    return -1.0;
	  else if(ind==1)
	    return 0.0;
	  else
	    return 1.0;
	}
      else if(order==2)
	{
	  if(ind==0)
	    return 4.0*x + 4.0*y - 3.0;
	  else if(ind==1)
	    return 0.0;
	  else if(ind==2)
	    return 4.0*y - 1.0;
	  else if(ind==3)
	    return -4.0 * x;
	  else if(ind==4)
	    return 4.0 * x;
	  else
	    return 4.0 * ( 1.0 - x - 2.0*y );
	}
      else
	{
	  if(ind==0)
	    return 0.5 * ( -27.0*x*x - 27.0*y*y - 54.0*x*y + 36.0*x + 36.0*y - 11.0 );
	  else if(ind==1)
	    return 0.0;
	  else if(ind==2)
	    return 1.0 - 9.0*y + 13.5*y*y;
	  else if(ind==3)
	    return 4.5 * x * ( -5.0 + 6.0*x + 6.0*y );
	  else if(ind==4)
	    return -4.5 * x * ( -1.0 + 3.0*x );
	  else if(ind==5)
	    return 4.5 * x * ( -1.0 + 3.0*x );
	  else if(ind==6)
	    return 4.5 * x * ( -1.0 + 6.0*y );
	  else if(ind==7)
	    return -4.5 * ( 1.0 + 9.0*y*y - x - 8.0*y + 6.0*x*y);
	  else if(ind==8)
	    return 4.5 * ( 2.0 + 9.0*y*y - 5.0*x + 3.0*x*x - 10.0*y + 12.0*x*y );
	  else
	    return 27.0 * x * ( 1.0 - 2.0*y - x ); 
	}
    }
}

//=========================================================

BasisFunctions::~BasisFunctions()
{
  ;
}
