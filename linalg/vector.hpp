#ifndef FILE_VECTOR
#define FILE_VECTOR

/*
   Data type vector

   Author: Tzanio, Joachim
   Date:   01/18/2001
*/

using namespace std;

#include "../general/array.hpp"

/// Vector data type. 
class Vector
{
protected:

  int size;
  double * data;
  
public:
  
  /// Default constructor for Vector. Sets size = 0 and data = NULL
  Vector () { size = 0; data = 0; };
  
  /// Creates vector of size s. 
  Vector (int s);

  /// Reads a vector from multpile files
  void Load (istream ** in, int np, int * dim);
  
  /// Load a vector from an input stream.
  void Load(istream &in, int Size);
  
  /// Load a vector from an input stream.
  void Load(istream &in) { int s; in >> s; Load (in, s); };

  /// Resizes the vector if the new size is different
  void SetSize(int s);
  
  /// Returns the size of the vector. 
  inline int Size() const {return size;};

//  double *GetData() { return data; }

  double *GetData() const { return data; }

  /// Sets value in vector. Index i = 0 .. size-1 
  double & Elem (int i);

  /// Sets value in vector. Index i = 0 .. size-1 
  const double & Elem (int i) const;

  /// Sets value in vector. Index i = 0 .. size-1 
  inline double & operator() (int i);

  /// Sets value in vector. Index i = 0 .. size-1 
  inline const double & operator() (int i) const;

  /// Return the inner-product.
  double operator*(const Vector & v) const;

  /// Return the l-infinity norm of a vector.
  double InfinityNorm() const;

  /// Redefine '=' for vector = vector. 
  Vector & operator=(const Vector &v);

  /// Redefine '=' for vector = constant.
  Vector & operator=(double value);

  Vector & operator*=(double c);

  Vector & operator-=(Vector &v);

  Vector & operator+=(Vector &v);

  /// (*this) += a * Va
  Vector & Add(double a, Vector &Va);

  /// (*this) = -(*this)
  void Neg();

  /// Do v = v1 - v2.
  friend void subtract(const Vector &v1, const Vector &v2, Vector &v);
  
  /// Do v = v1 + v2.
  friend void add(const Vector &v1, const Vector &v2, Vector &v);

  /// Do v = v1 + alpha * v2.
  friend void add(const Vector &v1,double alpha,const Vector &v2,Vector &v);

  void GetSubVector(const Array<int> &dofs, Vector &elemvect) const;

  void SetSubVector(const Array<int> &dofs, const Vector &elemvect);

  /// Add (element) subvector to the vector.
  void AddElementVector(const Array<int> & dofs, const Vector & elemvect);

  /// Prints vector to stream out.
  void Print (ostream & out = cout, int width =8);

  /// Prints vector to stream out in HYPRE_Vector format.
  void Print_HYPRE (ostream &out);

  /// Destroys vector. 
  ~Vector ();
};

class SubVector : public Vector
{
   public:
      SubVector (Vector &v, int beg, int _size)
         { data = v.GetData() + beg; size = _size; };
      SubVector (const Vector &v, int beg, int _size)
         { data = v.GetData() + beg; size = _size; };
      void NewSubVector (Vector &v, int beg, int _size)
         { data = v.GetData() + beg; size = _size; };
      SubVector & operator=(const Vector &v)
         { Vector::operator=(v); return (*this); };
      ~SubVector () { data = NULL; };
};

//--------------------------------------------------------------------


inline double & Vector::operator() (int i)
{
#ifdef DEBUG
  if (data == 0 || i < 0 || i >= size)
    {
      cerr << "Vector::operator()" << endl;
      return data[0];
    }
#endif
  
  return data[i];
}

inline const double & Vector::operator() (int i) const
{
#ifdef DEBUG
  if (data == 0 || i < 0 || i >= size)
    {
      cerr << "Vector::operator()" << endl;
      return data[0];
    }
#endif
  
  return data[i];
}

#endif
