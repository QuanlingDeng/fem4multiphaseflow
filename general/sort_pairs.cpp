
#include <stdlib.h>
#include "sort_pairs.hpp"


template <class A, class B>
int ComparePairs (const void *_p, const void *_q)
{
   Pair<A, B> *p, *q;

   p = (Pair<A, B> *)_p;
   q = (Pair<A, B> *)_q;

   if (p -> one < q -> one)  return -1;
   if (q -> one < p -> one)  return +1;
   return 0;
}

template <class A, class B>
void SortPairs (Pair<A, B> *pairs, int size)
{
   qsort (pairs, size, sizeof(Pair<A, B>), ComparePairs<A, B>);
}


template int ComparePairs<int, int> (const void *, const void *);
template void SortPairs<int, int> (Pair<int, int> *, int );
