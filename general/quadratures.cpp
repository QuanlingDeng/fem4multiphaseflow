#include "fem.h"


/*  
  Author: Q. Deng
  Date: 05/24/2015
*/

QuadratureRule::QuadratureRule( int NP )
{
  NPoints   = NP;
  QuadraturePoints = new QuadraturePoint[NP];
}

//============================================================================

void QuadratureRule::GaussianRule()
{
  int n = NPoints;
  int m = (n+1)/2;
  int i, j;
  double p1, p2, p3;
  double pp, z, z1;
  for (i = 1; i <= m; i++)
    {
      z = cos ( M_PI * (i - 0.25) / (n + 0.5));
      
      while(1)
        {
          p1 = 1;
          p2 = 0;
          for (j = 1; j <= n; j++)
            {
              p3 = p2;
              p2 = p1;
              p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
            }
          // p1 is Legendre polynomial
	  
          pp = n * (z*p1-p2) / (z*z - 1);
          z1 = z;
          z = z1-p1/pp;
	  
          if (fabs (z - z1) < 1e-14) break;
        }
      
      QuadraturePoints[i-1].x  = 0.5 * (1 - z);
      QuadraturePoints[n-i].x  = 0.5 * (1 + z);
      QuadraturePoints[i-1].weight = QuadraturePoints[n-i].weight = 
        1.0 / ( (1  - z * z) * pp * pp);
    }
}

//============================================================================

QuadratureRule::~QuadratureRule()
{
  delete []QuadraturePoints;
}

//============================================================================

QuadratureRules QuadRules;

QuadratureRules::QuadratureRules()
{
  cout << "Initializing IntRules ..." << endl;

  SegmentQuadratureRules();
  TriangleQuadratureRules();
  SquareQuadratureRules();
}

//============================================================================

QuadratureRule &QuadratureRules::Get( int GeomType, int Order )
{
  switch( GeomType )
    {
    case Geometry::SEGMENT:  return *SegmentQuadRules[Order];
    case Geometry::TRIANGLE: return *TriangleQuadRules[Order];
    case Geometry::SQUARE:   return *SquareQuadRules[Order];
    default:
#ifdef DEBUG
      cerr << "QuadratureRules::Get";
      return *TriangleQuadRules[Order];
#endif
      ;
    }
  return *TriangleQuadRules[Order]; 
}

//============================================================================

QuadratureRules::~QuadratureRules ()
{
  int i;

  for( i=0; i<TriangleQuadRules.Size(); i++ )
    if( TriangleQuadRules[i]!=NULL )
      delete TriangleQuadRules[i];
  
  for( i=0; i<SquareQuadRules.Size(); i++ )
    if( SquareQuadRules[i]!=NULL )
      delete SquareQuadRules[i];
}

//===========================================================================

// Quadrature rules for line segment [0,1]
void QuadratureRules::SegmentQuadratureRules()
{
  int i;

  SegmentQuadRules.SetSize(32);

  for( i=0; i<SegmentQuadRules.Size(); i++ )
    SegmentQuadRules[i] = NULL;

  // 1 point - 1 degree
  SegmentQuadRules[0] = new QuadratureRule(1);

  SegmentQuadRules[0] -> QuadPoint(0).x = .5;
  SegmentQuadRules[0] -> QuadPoint(0).weight = 1.;

  // 1 point - 1 degree
  SegmentQuadRules[1] = new QuadratureRule(1);

  SegmentQuadRules[1] -> QuadPoint(0).x = .5;
  SegmentQuadRules[1] -> QuadPoint(0).weight = 1.;   

  // 2 point - 3 degree
  SegmentQuadRules[2] = new QuadratureRule(2);

  SegmentQuadRules[2] -> QuadPoint(0).x = 0.211324865405187;
  SegmentQuadRules[2] -> QuadPoint(0).weight = .5;   
  SegmentQuadRules[2] -> QuadPoint(1).x = 0.788675134594812;
  SegmentQuadRules[2] -> QuadPoint(1).weight = .5;   

  // 2 point - 3 degree
  SegmentQuadRules[3] = new QuadratureRule(2);

  SegmentQuadRules[3] -> QuadPoint(0).x = 0.211324865405187;
  SegmentQuadRules[3] -> QuadPoint(0).weight = .5;   
  SegmentQuadRules[3] -> QuadPoint(1).x = 0.788675134594812;
  SegmentQuadRules[3] -> QuadPoint(1).weight = .5;   

  // 3 point - 5 degree
  SegmentQuadRules[4] = new QuadratureRule(3);
  SegmentQuadRules[4] -> GaussianRule();

  SegmentQuadRules[4] -> QuadPoint(0).x = 0.11270166537925831148;
  SegmentQuadRules[4] -> QuadPoint(0).weight = 0.2777777777777777777777778;   
  SegmentQuadRules[4] -> QuadPoint(1).x = 0.5;
  SegmentQuadRules[4] -> QuadPoint(1).weight = 0.4444444444444444444444444;
  SegmentQuadRules[4] -> QuadPoint(2).x = 0.88729833462074168852;
  SegmentQuadRules[4] -> QuadPoint(2).weight = 0.2777777777777777777777778;

  // 3 point - 5 degree
  SegmentQuadRules[5] = new QuadratureRule (3);

  SegmentQuadRules[5] -> QuadPoint(0).x = 0.11270166537925831148;
  SegmentQuadRules[5] -> QuadPoint(0).weight = 0.2777777777777777777777778;   
  SegmentQuadRules[5] -> QuadPoint(1).x = 0.5;
  SegmentQuadRules[5] -> QuadPoint(1).weight = 0.4444444444444444444444444;
  SegmentQuadRules[5] -> QuadPoint(2).x = 0.88729833462074168852;
  SegmentQuadRules[5] -> QuadPoint(2).weight = 0.2777777777777777777777778;

  for( i=6; i<SegmentQuadRules.Size(); i++ )
    {
      SegmentQuadRules[i] = new QuadratureRule(i/2+1);
      SegmentQuadRules[i] -> GaussianRule();
    }
}

//===========================================================================

// Quadrature rules for reference triangle {[0,0],[1,0],[0,1]}
void QuadratureRules::TriangleQuadratureRules()
{
  TriangleQuadRules.SetSize(9);

   for( int i=0; i<TriangleQuadRules.Size(); i++ )
     TriangleQuadRules[i] = NULL;

  // 1 point - 0 degree 
  TriangleQuadRules[0] = new QuadratureRule(1);

  TriangleQuadRules[0] -> QuadPoint(0).x = 0.33333333333333;
  TriangleQuadRules[0] -> QuadPoint(0).y = 0.33333333333333;
  TriangleQuadRules[0] -> QuadPoint(0).weight = 0.5;

  // 3 point - 2 degree (vertices)
  TriangleQuadRules[1] = new QuadratureRule (3);

  TriangleQuadRules[1] -> QuadPoint(0).x      = 0.;
  TriangleQuadRules[1] -> QuadPoint(0).y      = 0.;
  TriangleQuadRules[1] -> QuadPoint(0).weight = 0.16666666666667;

  TriangleQuadRules[1] -> QuadPoint(1).x      = 1.;
  TriangleQuadRules[1] -> QuadPoint(1).y      = 0.;
  TriangleQuadRules[1] -> QuadPoint(1).weight = 0.16666666666667;

  TriangleQuadRules[1] -> QuadPoint(2).x      = 0.;
  TriangleQuadRules[1] -> QuadPoint(2).y      = 1.;
  TriangleQuadRules[1] -> QuadPoint(2).weight = 0.16666666666667;

  // 3 point - 2 degree (midpoints)
  TriangleQuadRules[2] = new QuadratureRule (3);

  TriangleQuadRules[2] -> QuadPoint(0).x      = 0.5;
  TriangleQuadRules[2] -> QuadPoint(0).y      = 0;
  TriangleQuadRules[2] -> QuadPoint(0).weight = 0.16666666666667;

  TriangleQuadRules[2] -> QuadPoint(1).x      = 0.5;
  TriangleQuadRules[2] -> QuadPoint(1).y      = 0.5;
  TriangleQuadRules[2] -> QuadPoint(1).weight = 0.16666666666667;

  TriangleQuadRules[2] -> QuadPoint(2).x      = 0;
  TriangleQuadRules[2] -> QuadPoint(2).y      = 0.5;
  TriangleQuadRules[2] -> QuadPoint(2).weight = 0.16666666666667;

  // 4 point - 3 degree
  TriangleQuadRules[3] = new QuadratureRule (4);

  TriangleQuadRules[3] -> QuadPoint(0).x      = 0.33333333333333;
  TriangleQuadRules[3] -> QuadPoint(0).y      = 0.33333333333333;
  TriangleQuadRules[3] -> QuadPoint(0).weight = -0.28125;

  TriangleQuadRules[3] -> QuadPoint(1).x      = 0.2;
  TriangleQuadRules[3] -> QuadPoint(1).y      = 0.2;
  TriangleQuadRules[3] -> QuadPoint(1).weight = 0.26041666666665;

  TriangleQuadRules[3] -> QuadPoint(2).x      = 0.6;
  TriangleQuadRules[3] -> QuadPoint(2).y      = 0.2;
  TriangleQuadRules[3] -> QuadPoint(2).weight = 0.26041666666665;

  TriangleQuadRules[3] -> QuadPoint(3).x      = 0.2;
  TriangleQuadRules[3] -> QuadPoint(3).y      = 0.6;
  TriangleQuadRules[3] -> QuadPoint(3).weight = 0.26041666666665;

  // 6 point - 4 degree
  TriangleQuadRules[4] = new QuadratureRule (6);

  TriangleQuadRules[4] -> QuadPoint(0).x      = 0.091576213509771;
  TriangleQuadRules[4] -> QuadPoint(0).y      = 0.091576213509771;
  TriangleQuadRules[4] -> QuadPoint(0).weight = 0.054975871827661;

  TriangleQuadRules[4] -> QuadPoint(1).x      = 0.091576213509771;
  TriangleQuadRules[4] -> QuadPoint(1).y      = 0.816847572980459;
  TriangleQuadRules[4] -> QuadPoint(1).weight = 0.054975871827661;

  TriangleQuadRules[4] -> QuadPoint(2).x      = 0.816847572980459;
  TriangleQuadRules[4] -> QuadPoint(2).y      = 0.091576213509771;
  TriangleQuadRules[4] -> QuadPoint(2).weight = 0.054975871827661;

  TriangleQuadRules[4] -> QuadPoint(3).x      = 0.445948490915965;
  TriangleQuadRules[4] -> QuadPoint(3).y      = 0.445948490915965;
  TriangleQuadRules[4] -> QuadPoint(3).weight = 0.1116907948390055;

  TriangleQuadRules[4] -> QuadPoint(4).x      = 0.445948490915965;
  TriangleQuadRules[4] -> QuadPoint(4).y      = 0.108103018168070;
  TriangleQuadRules[4] -> QuadPoint(4).weight = 0.1116907948390055;

  TriangleQuadRules[4] -> QuadPoint(5).x      = 0.108103018168070;
  TriangleQuadRules[4] -> QuadPoint(5).y      = 0.445948490915965;
  TriangleQuadRules[4] -> QuadPoint(5).weight = 0.1116907948390055;

  // 7 point - 5 degree
  TriangleQuadRules[5] = new QuadratureRule (7);

  TriangleQuadRules[5] -> QuadPoint(0).x      = 0.3333333333333333333333333333333;
  TriangleQuadRules[5] -> QuadPoint(0).y      = 0.3333333333333333333333333333333;
  TriangleQuadRules[5] -> QuadPoint(0).weight = 0.1125;

  TriangleQuadRules[5] -> QuadPoint(1).x      = 0.1012865073234563388009873619151;
  TriangleQuadRules[5] -> QuadPoint(1).y      = 0.1012865073234563388009873619151;
  TriangleQuadRules[5] -> QuadPoint(1).weight = 0.06296959027241357629784197275009;

  TriangleQuadRules[5] -> QuadPoint(2).x      = 0.1012865073234563388009873619151;
  TriangleQuadRules[5] -> QuadPoint(2).y      = 0.7974269853530873223980252761698;
  TriangleQuadRules[5] -> QuadPoint(2).weight = 0.06296959027241357629784197275009;

  TriangleQuadRules[5] -> QuadPoint(3).x      = 0.7974269853530873223980252761698;
  TriangleQuadRules[5] -> QuadPoint(3).y      = 0.1012865073234563388009873619151;
  TriangleQuadRules[5] -> QuadPoint(3).weight = 0.06296959027241357629784197275009;

  TriangleQuadRules[5] -> QuadPoint(4).x      = 0.4701420641051150897704412095134;
  TriangleQuadRules[5] -> QuadPoint(4).y      = 0.4701420641051150897704412095134;
  TriangleQuadRules[5] -> QuadPoint(4).weight = 0.06619707639425309036882469391658;

  TriangleQuadRules[5] -> QuadPoint(5).x      = 0.4701420641051150897704412095134;
  TriangleQuadRules[5] -> QuadPoint(5).y      = 0.0597158717897698204591175809731;
  TriangleQuadRules[5] -> QuadPoint(5).weight = 0.06619707639425309036882469391658;

  TriangleQuadRules[5] -> QuadPoint(6).x      = 0.0597158717897698204591175809731;
  TriangleQuadRules[5] -> QuadPoint(6).y      = 0.4701420641051150897704412095134;
  TriangleQuadRules[5] -> QuadPoint(6).weight = 0.06619707639425309036882469391658;

  // 12 point - 6 degree
  TriangleQuadRules[6] = new QuadratureRule (12);

  TriangleQuadRules[6] -> QuadPoint(0).x      = 0.063089014491502;
  TriangleQuadRules[6] -> QuadPoint(0).y      = 0.063089014491502;
  TriangleQuadRules[6] -> QuadPoint(0).weight = 0.0254224531851035;

  TriangleQuadRules[6] -> QuadPoint(1).x      = 0.063089014491502;
  TriangleQuadRules[6] -> QuadPoint(1).y      = 0.873821971016996;
  TriangleQuadRules[6] -> QuadPoint(1).weight = 0.0254224531851035;

  TriangleQuadRules[6] -> QuadPoint(2).x      = 0.873821971016996;
  TriangleQuadRules[6] -> QuadPoint(2).y      = 0.063089014491502;
  TriangleQuadRules[6] -> QuadPoint(2).weight = 0.0254224531851035;

  TriangleQuadRules[6] -> QuadPoint(3).x      = 0.249286745170911;
  TriangleQuadRules[6] -> QuadPoint(3).y      = 0.249286745170911;
  TriangleQuadRules[6] -> QuadPoint(3).weight = 0.0583931378631895;

  TriangleQuadRules[6] -> QuadPoint(4).x      = 0.249286745170911;
  TriangleQuadRules[6] -> QuadPoint(4).y      = 0.501426509658179;
  TriangleQuadRules[6] -> QuadPoint(4).weight = 0.0583931378631895;

  TriangleQuadRules[6] -> QuadPoint(5).x      = 0.501426509658179;
  TriangleQuadRules[6] -> QuadPoint(5).y      = 0.249286745170911;
  TriangleQuadRules[6] -> QuadPoint(5).weight = 0.0583931378631895;

  TriangleQuadRules[6] -> QuadPoint(6).x      = 0.310352451033785;
  TriangleQuadRules[6] -> QuadPoint(6).y      = 0.053145049844816;
  TriangleQuadRules[6] -> QuadPoint(6).weight = 0.041425537809187;

  TriangleQuadRules[6] -> QuadPoint(7).x      = 0.310352451033785;
  TriangleQuadRules[6] -> QuadPoint(7).y      = 0.636502499121399;
  TriangleQuadRules[6] -> QuadPoint(7).weight = 0.041425537809187;

  TriangleQuadRules[6] -> QuadPoint(8).x      = 0.053145049844816;
  TriangleQuadRules[6] -> QuadPoint(8).y      = 0.310352451033785;
  TriangleQuadRules[6] -> QuadPoint(8).weight = 0.041425537809187;

  TriangleQuadRules[6] -> QuadPoint(9).x      = 0.053145049844816;
  TriangleQuadRules[6] -> QuadPoint(9).y      = 0.636502499121399;
  TriangleQuadRules[6] -> QuadPoint(9).weight = 0.041425537809187;

  TriangleQuadRules[6] -> QuadPoint(10).x      = 0.636502499121399;
  TriangleQuadRules[6] -> QuadPoint(10).y      = 0.310352451033785;
  TriangleQuadRules[6] -> QuadPoint(10).weight = 0.041425537809187;

  TriangleQuadRules[6] -> QuadPoint(11).x      = 0.636502499121399;
  TriangleQuadRules[6] -> QuadPoint(11).y      = 0.053145049844816;
  TriangleQuadRules[6] -> QuadPoint(11).weight = 0.041425537809187;

  // 13 point - 7 degree 
  TriangleQuadRules[7] = new QuadratureRule(13);

  TriangleQuadRules[7] -> QuadPoint(0).x      = 0.33333333333333;
  TriangleQuadRules[7] -> QuadPoint(0).y      = 0.33333333333333;
  TriangleQuadRules[7] -> QuadPoint(0).weight = -0.074785022233835;

  TriangleQuadRules[7] -> QuadPoint(1).x      = 0.2603459661;
  TriangleQuadRules[7] -> QuadPoint(1).y      = 0.2603459661;
  TriangleQuadRules[7] -> QuadPoint(1).weight = 0.087807628716602;

  TriangleQuadRules[7] -> QuadPoint(2).x      = 0.4793080678;
  TriangleQuadRules[7] -> QuadPoint(2).y      = 0.2603459661;
  TriangleQuadRules[7] -> QuadPoint(2).weight = 0.087807628716602;

  TriangleQuadRules[7] -> QuadPoint(3).x      = 0.2603459661;
  TriangleQuadRules[7] -> QuadPoint(3).y      = 0.4793080678;
  TriangleQuadRules[7] -> QuadPoint(3).weight = 0.087807628716602;

  TriangleQuadRules[7] -> QuadPoint(4).x      = 0.0651301029;
  TriangleQuadRules[7] -> QuadPoint(4).y      = 0.0651301029;
  TriangleQuadRules[7] -> QuadPoint(4).weight = 0.026673617804419;

  TriangleQuadRules[7] -> QuadPoint(5).x      = 0.8697397942;
  TriangleQuadRules[7] -> QuadPoint(5).y      = 0.0651301029;
  TriangleQuadRules[7] -> QuadPoint(5).weight = 0.026673617804419;

  TriangleQuadRules[7] -> QuadPoint(6).x      = 0.0651301029;
  TriangleQuadRules[7] -> QuadPoint(6).y      = 0.8697397942;
  TriangleQuadRules[7] -> QuadPoint(6).weight = 0.026673617804419;

  TriangleQuadRules[7] -> QuadPoint(7).x      = 0.31286549600487;
  TriangleQuadRules[7] -> QuadPoint(7).y      = 0.04869031542532;
  TriangleQuadRules[7] -> QuadPoint(7).weight = 0.038556880445128;

  TriangleQuadRules[7] -> QuadPoint(8).x      = 0.04869031542532;
  TriangleQuadRules[7] -> QuadPoint(8).y      = 0.31286549600487;
  TriangleQuadRules[7] -> QuadPoint(8).weight = 0.038556880445128;

  TriangleQuadRules[7] -> QuadPoint(9).x      = 0.63844418856981;
  TriangleQuadRules[7] -> QuadPoint(9).y      = 0.04869031542532;
  TriangleQuadRules[7] -> QuadPoint(9).weight = 0.038556880445128;

  TriangleQuadRules[7] -> QuadPoint(10).x      = 0.63844418856981;
  TriangleQuadRules[7] -> QuadPoint(10).y      = 0.31286549600487;
  TriangleQuadRules[7] -> QuadPoint(10).weight = 0.038556880445128;

  TriangleQuadRules[7] -> QuadPoint(11).x      = 0.04869031542532;
  TriangleQuadRules[7] -> QuadPoint(11).y      = 0.63844418856981;
  TriangleQuadRules[7] -> QuadPoint(11).weight = 0.038556880445128;

  TriangleQuadRules[7] -> QuadPoint(12).x      = 0.31286549600487;
  TriangleQuadRules[7] -> QuadPoint(12).y      = 0.63844418856981;
  TriangleQuadRules[7] -> QuadPoint(12).weight = 0.038556880445128;

  // Used only for computing max norm.
  int n = 20;
  int npp = (n+1)*(n+1) - (n+1)*n/2;
  TriangleQuadRules[8] = new QuadratureRule(npp);
  double dx, dy;
  dx = dy = 1.0 / (double) n;
  int k = 0;
  for (int j=0; j<=n; j++)
    {
      double yy = dy * (double) j;
      for (int i=0; i<=(n-j); i++)
	{
	  double xx = dx * (double) i;
	  TriangleQuadRules[8] -> QuadPoint(k).x      = xx;
	  TriangleQuadRules[8] -> QuadPoint(k).y      = yy;
	}
    }
}

//============================================================================

// Quadrature rules for unit square
void QuadratureRules::SquareQuadratureRules()
{
  SquareQuadRules.SetSize(14);

  int i,k,s,np;
  
  for( i=0; i<SquareQuadRules.Size(); i++ )
    {
      np = SegmentQuadRules[i]->GetNPoints();
      SquareQuadRules[i] = new QuadratureRule(np*np);
      for( k=0; k<np; k++ )
	for( s=0; s<np; s++ )
	  {
	    SquareQuadRules[i]->QuadPoint(k*np+s).x 
	      = SegmentQuadRules[i]->QuadPoint(k).x; 

	    SquareQuadRules[i]->QuadPoint(k*np+s).y 
	      = SegmentQuadRules[i]->QuadPoint(s).x; 

	    SquareQuadRules[i]->QuadPoint(k*np+s).weight 
	      = SegmentQuadRules[i]->QuadPoint(k).weight 
	      * SegmentQuadRules[i]->QuadPoint(s).weight;
	  }	        
    }
}

