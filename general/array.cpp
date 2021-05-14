/**************************************************************************/
/* File:   array.cc                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/* 
   Abstract data type ARRAY
*/

using namespace std;

#include "array.hpp"
#include <string.h>

BaseArray :: BaseArray(int asize, int ainc, int elementsize)
{
 
  void * p;

  if (asize)
    {
      p = new char[asize * elementsize];

      data = p;
      size = allocsize = asize; inc = ainc;
      return;
    }

  data = 0;
  size = allocsize = 0;
  inc = ainc;
}

BaseArray :: ~BaseArray()
{
  if (data)
    delete (char*)data;
}

void BaseArray :: ReSize (int minsize, int elementsize)
{
  void * p;
  int nsize = (inc) ? allocsize + inc : 2 * allocsize;
  if (nsize < minsize) nsize = minsize;

  if (data)
    {
      p = new char [nsize * elementsize];
      memcpy (p, data, ((nsize < size) ? nsize : size) * elementsize);


      delete (char*)data;
      data = p;
    }
  else
    {
      p = new char[nsize * elementsize];
      data = p;
    }

  allocsize = nsize;
}

template <class T>
void Array<T>::Save(ostream & out)
{
   int i;

   out << size << '\n';

   for (i = 0; i < size; i++)
      out << operator[](i) << '\n';
}

template <class T>
T Array<T>::Max()
{
   int i;
   T max;

   if (size > 0)
      max = operator[](0);
   for (i = 1; i < size; i++)
      if (max < operator[](i))
         max = operator[](i);

   return max;
}

template <class T>
void swap (T * v,int i,int j)
{
  static T t;
  t = v[i];  v[i] = v[j];  v[j] = t;
}

template <class T>
void qsort (T * v,int l,int r)
{
  int i, last;
  if (l >= r)
    return;
  swap (v,l,(l+r)/2);
  last = l;
  for (i = l+1; i <= r; i++)
    if (v[i] < v[l])
      swap(v,++last,i);
  swap (v,l,last);
  qsort (v,l,last-1);
  qsort (v,last+1,r);
}

template <class T>
void Array<T>::Sort()
{
  qsort((T*)data,0,size-1);
}

template class Array<int>;
template class Array<double>;
