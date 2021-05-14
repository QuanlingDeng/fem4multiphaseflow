#include <iostream>
#include <iomanip>
#include <math.h>
#include "vector.hpp"

using namespace std;

Vector::Vector (int s)
{
  size = s;
  if (s <= 0)
     data = 0;
  else
     data = new double[s];
}

void Vector::Load (istream ** in, int np, int * dim) {
  int i,j,s;

  s = 0;
  for (i = 0; i < np; i++)
    s += dim[i];

  SetSize (s);

  int p=0;
  for (i = 0; i < np; i++)
    for (j = 0; j < dim[i]; j++)
      *in[i] >> data[p++];
}


void Vector::Load(istream &in, int Size){
  SetSize (Size);

  for(int i=0; i<size; i++)
    in >> data[i];
}

void Vector::SetSize(int s)
{
   if (size == s)
      return;
   if (data != 0)
     delete [] data;
   size = s;
   if (s <= 0)
      data = 0;
   else
      data = new double[s];
}

double & Vector::Elem (int i)
{
  return operator()(i);
}

const double & Vector::Elem (int i) const
{
  return operator()(i);
}

//============================================================================

double Vector::operator*(const Vector & v) const{
#ifdef DEBUG
  if (v.size != size)
    cerr << "Vector::operator*\n";
#endif
      
  double sum = 0.;
  for(int i=0; i<size; i++)
    sum += v.data[i]*data[i];
  return sum;
}

//============================================================================

double Vector::InfinityNorm() const{
  double temp = 0.0, datAbs;
  for(int i=0; i<size; i++)
    {
      datAbs = fabs( data[i] );
      temp = ( temp>datAbs ) ? ( temp ) : datAbs;
    }
  return temp;
}

//============================================================================

Vector & Vector::operator=(const Vector &v){
  /*  if (size != (size = v.Size())){
    if (data != 0)
      delete [] data;
    data = new double[size];
  }
  for(int i=0; i<size; i++)
    data[i] = v.data[i];
  return *this;
  */
  this->SetSize(v.Size());
  for(int i=0; i<size; i++)
    data[i] = v.data[i];
  return *this;


}

//============================================================================

Vector & Vector::operator=(double value){
  for(int i =0; i<size; i++)
    data[i] = value;
  return *this;
}


//============================================================================

Vector & Vector::operator*=(double c)
{
   for(int i = 0; i < size; i++)
      data[i] *= c;
   return *this;
}

//============================================================================

Vector & Vector::operator-=(Vector &v)
{
   for(int i = 0; i < size; i++)
      data[i] -= v(i);
   return *this;
}

//============================================================================

Vector & Vector::operator+=(Vector &v)
{
   for(int i = 0; i < size; i++)
      data[i] += v(i);
   return *this;
}

//============================================================================

Vector & Vector::Add(double a, Vector &Va)
{
   for(int i = 0; i < size; i++)
      data[i] += a * Va(i);
   return *this;
}

//============================================================================

void Vector::Neg()
{
   for(int i = 0; i < size; i++)
      data[i] = -data[i];
}

//============================================================================

void subtract(const Vector &v1, const Vector &v2, Vector &v){
#ifdef DEBUG 
  if (v.size != v1.size || v.size != v2.size)
    cerr << "subtract\n";
#endif
  for(int i=0; i<v.size; i++)
    v.data[i] = v1.data[i] - v2.data[i];
}
  
//============================================================================

void add(const Vector &v1, const Vector &v2, Vector &v){
#ifdef DEBUG 
  if (v.size != v1.size || v.size != v2.size)
    cerr << "add(Vector &v1, Vector &v2, Vector &v)\n";
#endif
  for(int i=0; i<v.size; i++)
    v.data[i] = v1.data[i] + v2.data[i];
}

//============================================================================

void add(const Vector &v1,double alpha,const Vector &v2,Vector &v){
#ifdef DEBUG 
  if (v.size != v1.size || v.size != v2.size)
    cerr << "add(Vector &v1, double alpha, Vector &v2, Vector &v)\n";
#endif
  for(int i=0; i<v.size; i++)
    v.data[i] = v1.data[i] + alpha*v2.data[i];
}

//============================================================================

void Vector::GetSubVector(const Array<int> & dofs, Vector &elemvect) const
{
   int i, j, n = dofs.Size();

   elemvect.SetSize (n);

   for (i = 0; i < n; i++)
      if ((j=dofs[i]) >= 0)
         elemvect(i) = data[j];
      else
         elemvect(i) = -data[-1-j];
}

//============================================================================

void Vector::SetSubVector(const Array<int> & dofs, const Vector &elemvect)
{
   int i, j, n = dofs.Size();

   for (i = 0; i < n; i++)
      if ((j=dofs[i]) >= 0)
         data[j] = elemvect(i);
      else
         data[-1-j] = -elemvect(i);
}

//============================================================================

void Vector::AddElementVector(const Array<int> & dofs, const Vector & elemvect)
{
   int i, j, n = dofs.Size();

   for (i = 0; i < n; i++)
      if ((j=dofs[i]) >= 0)
         data[j] += elemvect(i);
      else
         data[-1-j] -= elemvect(i);
}

//============================================================================

void Vector::Print (ostream & out, int width)
{
  // output flags = scientific + show sign 
  out << setiosflags(ios::scientific | ios::showpos); 
  for (int i = 0; i < size; i++) 
    {
      out << data[i] << " ";
      if ( !((i+1) % width) )
	out << endl;
    }
  out << endl;
}

void Vector::Print_HYPRE (ostream &out)
{
   int i;
   ios::fmtflags old_fmt = out.setf(ios::scientific);
   int old_prec = out.precision(14);

   out << size << '\n';  // number of rows

   for (i = 0; i < size; i++)
      out << data[i] << '\n';

   out.precision(old_prec);
   out.setf(old_fmt);
}

Vector::~Vector() 
{
   if (data != 0)
      delete [] data;
}
