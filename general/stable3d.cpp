
#include <iostream>
#include "stable3d.hpp"

using namespace std;

STable3D::STable3D (int nr)
{
   int i;

   Size = nr;
   Rows = new STable3DNode *[nr];
   for (i = 0; i < nr; i++)
      Rows[i] = NULL;
   NElem = 0;
}

inline void Sort3 (int &r, int &c, int &f)
{
   int t;

   if (r > c)
      if (c > f)
      {
         t = r;  r = f;  f = t;  //  (r,c,f) --> (f,c,r)
      }
      else
         if (r > f)
         {
            t = r;  r = c;  c = f;  f = t;  //  (r,c,f) --> (c,f,r)
         }
         else
         {
            t = r;  r = c;  c = t;  //  (r,c,f) --> (c,r,f)
         }
   else
      if (c > f)
         if (r > f)
         {
            t = f;  f = c;  c = r;  r = t;  //  (r,c,f) --> (f,r,c)
         }
         else
         {
            t = c;  c = f;  f = t;  //  (r,c,f) --> (r,f,c)
         }
}

int STable3D::Push (int r, int c, int f)
{
   STable3DNode *node;

#ifdef DEBUG
   if (r == c || c == f || f == c)
   {
      cerr << "STable3D::Push : r = " << r << ", c = " << c << ", f = "
           << f << endl;
      return *((int *)NULL);
   }
#endif

   Sort3 (r, c, f);

   for (node = Rows[r]; node != NULL; node = node->Prev)
   {
      if (node->Column == c)
         if (node->Floor == f)
            return node->Number;
   }

   node = NodesMem.Alloc ();
   node->Column = c;
   node->Floor  = f;
   node->Number = NElem;
   node->Prev   = Rows[r];
   Rows[r] = node;

   NElem++;
   return (NElem-1);
}

int STable3D::operator() (int r, int c, int f) const
{
   STable3DNode *node;

   Sort3 (r, c, f);

   for (node = Rows[r]; node != NULL; node = node->Prev)
   {
      if (node->Column == c)
         if (node->Floor == f)
            return node->Number;
   }

   cerr << "STable3D::operator(): (r,c,f) = (" << r << "," << c << ","
        << f << ")" << endl;

   return *((int *)NULL);
}

int STable3D::Push4 (int r, int c, int f, int t)
{
#ifdef DEBUG
   if (r == c || r == f || r == t || c == f || c == t || f == t)
   {
      cerr << "STable3D::Push4 : r = " << r << ", c = " << c << ", f = "
           << f << ", t = " << t << endl;
      return *((int *)NULL);
   }
#endif

  int i = 0;
  int max = r;

  if (max < c) max = c, i = 1;
  if (max < f) max = f, i = 2;
  if (max < t) max = t, i = 3;

  switch(i) {
  case 0:
    return Push (c,f,t);
  case 1:
    return Push (r,f,t);
  case 2:
    return Push (r,c,t);
  case 3:
    return Push (r,c,f);
  }

  return -1;
}

int STable3D::operator() (int r, int c, int f, int t) const
{
  int i = 0;
  int max = r;

  if (max < c) max = c, i = 1;
  if (max < f) max = f, i = 2;
  if (max < t) max = t, i = 3;

  switch(i) {
  case 0:
    return (*this)(c,f,t);
  case 1:
    return (*this)(r,f,t);
  case 2:
    return (*this)(r,c,t);
  case 3:
    return (*this)(r,c,f);
  }

  return -1;
}

STable3D::~STable3D ()
{
   delete [] Rows;
}


ITable::ITable (int ND)
{
   int i;

   NDomains = ND;
   Domain = new ITableNode *[ND];
   NumNeighbors = new int[ND];
   for (i = 0; i < ND; i++)
   {
      Domain[i] = NULL;
      NumNeighbors[i] = 0;
   }
   TotWeight = NInterfaces = 0;
}

int ITable::ConnectDomains (int dom1, int dom2, int weight)
{
   int t;
   ITableNode *node;

   if (dom1 > dom2)
      t = dom1, dom1 = dom2, dom2 = t;

   TotWeight += weight;

   for (node = Domain[dom1]; node != NULL; node = node -> Prev)
      if (node -> OtherDomain == dom2)
      {
         node -> NFaces++;
         node -> IntWeight += weight;
         return node -> InterfaceNumber;
      }

   NumNeighbors[dom1]++;
   NumNeighbors[dom2]++;

   node = NodesMem.Alloc ();
   node -> OtherDomain = dom2;
   node -> InterfaceNumber = NInterfaces;
   node -> Prev = Domain[dom1];
   node -> NFaces = 1;
   node -> IntWeight = weight;
   Domain[dom1] = node;

   node = NodesMem.Alloc ();
   node -> OtherDomain = dom1;
   node -> InterfaceNumber = NInterfaces;
   node -> Prev = Domain[dom2];
   Domain[dom2] = node;

   return NInterfaces++;
}

int ITable::InterfaceNo (int dom1, int dom2)
{
   int t;
   ITableNode *node;

   if (dom1 > dom2)
      t = dom1, dom1 = dom2, dom2 = t;

   for (node = Domain[dom1]; node != NULL; node = node -> Prev)
      if (node -> OtherDomain == dom2)
         return node -> InterfaceNumber;

   cerr << "ITable::InterfaceNo (...)" << endl;

   return -1;
}

void ITable::GetInterfaceOffsetsAndNFaces (int *offsets, int *interfaces)
{
   int i, dom;
   ITableNode *node;

   offsets[0] = 0;
   for (dom = 0; dom < NDomains; dom++)
      for (node = Domain[dom]; node != NULL; node = node -> Prev)
         if (dom < node -> OtherDomain)
            offsets[node -> InterfaceNumber+1] = node -> IntWeight + 1; // !!!

   for (i = 1; i < NInterfaces; i++)
      offsets[i+1] += offsets[i];

   for (dom = 0; dom < NDomains; dom++)
      for (node = Domain[dom]; node != NULL; node = node -> Prev)
         if (dom < node -> OtherDomain)
            interfaces[offsets[node -> InterfaceNumber]] = node -> NFaces;

   for (i = 0; i < NInterfaces; i++)
      offsets[i]++;
}

ITable::~ITable()
{
   delete [] NumNeighbors;
   delete [] Domain;
}
