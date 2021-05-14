#ifndef FILE_ARRAY
#define FILE_ARRAY

/**************************************************************************/
/* File:   array.hh                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

using namespace std;

#include <iostream>

/// Base class for array container.
class BaseArray
{
protected:
  /// Pointer to data
  void * data;
  /// Size of the array
  int size;
  /// Size of the allocated memory
  int allocsize;
  /** Increment of allocated memory on overflow,
    allocsize = 0 doubles the array */
  int inc;

  /// Creates array of asize elements of size elementsize,
  BaseArray(int asize, int ainc, int elmentsize);
  ///
  ~BaseArray();
  /** Resizes allocated array to max(allocsize, minsize) 
      elements of size elementsize */
  void ReSize (int minsize, int elementsize);
};


/** 
   Abstract data type Array.
   
   Array<T> is an automatically increasing array containing elements of the
   generic type T. The allocated size may be larger then the logical size
   of the array. 
   The elements can be accessed by the [] operator, the range is 0 to size-1. 
*/
template <class T> 
class Array : public BaseArray
{
public:
  /// Creates array of asize elements 
  inline Array(int asize = 0, int ainc = 0)
    : BaseArray (asize, ainc, sizeof (T)) { ; }
  
  ///
  inline ~Array () 
  { ; }

  /// Return the data as 'T *'
  inline operator T *() { return (T *)data; }

  /// Returns the data
  inline void * GetData() { return data; }

  /// Logical size of the array
  inline int Size() const 
  { return size; };
  
  /// Change logical size of the array, keep existing entries
  inline void SetSize(int nsize);
  
  /// Change allocated size 
  inline void SetAllocSize (int nallocsize);

  /// access element
  inline T & operator[] (int i);

  /// access const element
  inline const T & operator[] (int i) const;

  /// append element to array, resize if necessary
  inline int Append (const T & el);

  /// append another array to this array, resize if necessary
  inline int Append (const Array<T> & els);

  /// return the last element in the array
  inline T & Last();

  /// append element when it is not yet in the array, return index
  inline int Union (const T & el);

  /// delete element i, move last one to position i
  inline void DeleteElement (int i);

  /// delete last element of array
  inline void DeleteLast ();

  /// delete whole array
  inline void DeleteAll ();

  /// Prints array to stream out.
  inline void Print(ostream & out, int width);

  void Save(ostream & out);

  /** Finds the maximal element in the array.
      (uses the comparison operator '<' for class T)  */
  T Max();

  /// Sorts the array.
  void Sort();
  
private:
  /// array copy is not supported
  Array<T> & operator= (Array<T> &);
  /// array copy is not supported
  Array (const Array<T> &);
};







template <class T>
inline void Array<T> :: SetSize(int nsize)
{
  if (nsize > allocsize)
    ReSize (nsize, sizeof(T));
  size = nsize;
}

template <class T>
inline void Array<T> :: SetAllocSize (int nallocsize)
{
  if (nallocsize > allocsize)
    ReSize (nallocsize, sizeof(T));
}




template <class T>
inline T &  Array<T> :: operator[] (int i)
{
#ifdef DEBUG    
  if (i < 0 || i >= size)
  {
     cerr << "Access element " << i << " of array, size = " << size << endl;
     return ((T*)(NULL))[0];
  }
#endif
  return ((T*)data)[i];
}

template <class T>
inline const T & Array<T> :: operator[] (int i) const
{
#ifdef DEBUG    
    if (i < 0 || i >= size)
    {
       cerr << "Access element " << i << " of array, size = " << size << endl;
       return ((T*)(NULL))[0];
    }
#endif
    return ((T*)data)[i];
}


template <class T>
inline int Array<T> :: Append (const T & el)
{
  if (size == allocsize) 
    ReSize (size+1, sizeof (T));
  ((T*)data)[size] = el;
  size++;
  return size;
}


template <class T>
inline int Array<T> :: Append (const Array<T> & els)
{
  int i, old_size;

  old_size = size;
  SetSize (size + els.Size());
  for (i = 0; i < els.Size(); i++)
     ((T*)data)[old_size+i] = els[i];
  return size;
}


template <class T>
inline T & Array<T> :: Last()
{
#ifdef DEBUG    
  if (size < 1)
  {
     cerr << "Array<T> :: Last()" << endl;
     return ((T*)data)[0];
  }
#endif
  return ((T*)data)[size-1];
}


template <class T>
inline int Array<T> :: Union (const T & el)
{
  int i=0;
  while ((i<size)&&(((T*)data)[i]!=el)) i++;
  if (i==size) 
    Append(el);
  return i;
}


template <class T>
inline void Array<T> :: DeleteElement (int i)
{
#ifdef DEBUG    
  if (i < 0 || i >= size)
  {
    cerr << "Delete element " << i << " of array, size = " << size << endl;
    return;
  }
#endif
  
  ((T*)data)[i] = ((T*)data)[size-1];
  size--;
}

  
template <class T>
inline void Array<T> :: DeleteLast ()
{
  size--;
}

template <class T>
inline void Array<T> :: DeleteAll ()
{
  if (data) delete (char*)data;
  data = NULL;
  size = allocsize = 0;
}

template <class T>
inline void Array<T> :: Print (ostream & out, int width)
{
  for (int i = 0; i < size; i++) 
    {
      out << ((T*)data)[i]  << " ";
      if ( !((i+1) % width) )
	out << endl;
    }
  out << endl;
}

#endif

