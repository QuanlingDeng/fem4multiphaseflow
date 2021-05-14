#ifndef FILE_QUADRATURES
#define FILE_QUADRATURES

/// Class for integration point with weight
class QuadraturePoint
{
 public:
  double x, y, weight;
};

// Class for integration rule
class QuadratureRule
{
private:
  int NPoints;
  QuadraturePoint *QuadraturePoints;
  
public:
  QuadratureRule() {}

  /// Construct a quadrature rule with given number of points
  QuadratureRule(int NP);

  /// Computes Gaussian quadrature rule on (0,1) with NumPoints.
  void GaussianRule();

  /// Return the number of the points in the integration rule
  int GetNPoints() const { return NPoints; };
  
  /// Return a reference to the i-th quadrature point
  QuadraturePoint & QuadPoint( int i ) { return QuadraturePoints[i]; };
  
  /// Returns a const reference to the i-th quadrature point
  const QuadraturePoint &QuadPoint(int i) const { return QuadraturePoints[i];};  

  /// Destroys a QuadratureRule object
  ~QuadratureRule();
};

/// Container class for quadrature rules
class QuadratureRules
{
private:
  Array<QuadratureRule *> SegmentQuadRules;
  Array<QuadratureRule *> TriangleQuadRules;
  Array<QuadratureRule *> SquareQuadRules;
  
  void SegmentQuadratureRules();
  void TriangleQuadratureRules();
  void SquareQuadratureRules();

public:
  /// Defines all quadrature rules
  QuadratureRules();
  
  /// Returns an quadrature rule for given GeomType and Order.
  QuadratureRule &Get( int GeomType, int Order );
  
  /// Destroys an QuadratureRules object
  ~QuadratureRules();
};

/// A global object with all quadrature rules (defined in quadratures.cpp)
extern QuadratureRules QuadRules;

#endif
