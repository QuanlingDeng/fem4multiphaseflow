
#ifdef DEBUG
#include <iostream.h>
#endif


int sorted_ok (int n, int *p)
{
   for (int i = 1; i < n; i++)
      if (p[i-1] > p[i])
         return 0;
   return 1;
}


void copy (int n, int *src, int *dst)
{
   for (int i = 0; i < n; i++)
      dst[i] = src[i];
}


void merge (int np, int nq, int *p, int *q, int *merged)
{
   int ip, iq, j;
   int pp, qq;

   if (np == 0)
   {
      copy (nq, q, merged);
      return;
   }
   if (nq == 0)
   {
      copy (np, p, merged);
      return;
   }

   j = ip = iq = 0;

   pp = p[0];
   qq = q[0];
   while (1)
      if (pp <= qq)
      {
         merged[j++] = pp;
         if (++ip == np)
         {
            merged[j++] = qq;
            for (iq++; iq < nq; iq++)
               merged[j++] = q[iq];
            return;
         }
         pp = p[ip];
      }
      else
      {
         merged[j++] = qq;
         if (++iq == nq)
         {
            merged[j++] = pp;
            for (ip++; ip < np; ip++)
               merged[j++] = p[ip];
            return;
         }
         qq = q[iq];
      }
}


void mergesort (int np, int *p, int *aux)
{
   int s, log2_np = 0;

   for (s = 1; s < np; s *= 2)
      log2_np++;

   int *_p, *_q, *_t;

   _p = p;
   _q = aux;

   if (log2_np % 2)
      for (int o = 1; o < np; o += 2)
      {
         int p1 = _p[o-1];
         int p2 = _p[o];

         if (p1 > p2)
            _p[o-1] = p2, _p[o] = p1;
      }
   else
   {
      for (int o = 1; o < np; o += 2)
      {
         int p1 = _p[o-1];
         int p2 = _p[o];

         if (p1 > p2)
            _q[o-1] = p2, _q[o] = p1;
         else
            _q[o-1] = p1, _q[o] = p2;
      }
      if (np % 2)
         _q[np-1] = _p[np-1];
      _p = aux; _q = p;
   }

   for (s = 2; s < np; s *= 2)
   {
      // 's' is the size of the chunks that are already sorted
      for (int o = 0; 1; o += 2*s)
         if (o+2*s < np)
            merge (s, s, _p+o, _p+o+s, _q+o);
         else
         {
            if (o+s < np)
               merge (s, np-s-o, _p+o, _p+o+s, _q+o);
            else
               copy (np-o, _p+o, _q+o);
            break;
         }
      _t = _p; _p = _q; _q = _t;
   }

#ifdef DEBUG
   if (_p != p && np != 1)
      cerr << "mergesort (...)" << endl;

   if (!sorted_ok (np, p))
      cerr << "mergesort (...) !sorted_ok(...)" << endl;
#endif
}


void mergecsrrows (int nr, int *row, int *ii, int *jj,
                   int *merged, int *aux, int *i_merged)
{
   int s, log2_nr, nnr, im;

   if (nr == 1)
   {
      int r1 = row[0];
      int b1 = ii[r1];
      int s1 = ii[r1+1]-b1;

      copy (s1, jj+b1, merged);
      i_merged[0] = 0;
      i_merged[1] = s1;
      return;
   }

   for (s = 1, log2_nr = 0; s < nr; s *= 2)
      log2_nr++;

   int *p, *q, *t;

   if (log2_nr % 2)
      p = aux, q = merged;
   else
      p = merged, q = aux;

   nnr = 0;
   i_merged[0] = im = 0;
   for (int o = 1; o < nr; o += 2)
   {
      int r1 = row[o-1];
      int r2 = row[o];
      int b1 = ii[r1];
      int b2 = ii[r2];
      int s1 = ii[r1+1]-b1;
      int s2 = ii[r2+1]-b2;

      merge (s1, s2, jj+b1, jj+b2, q+im);
      i_merged[++nnr] = (im += s1 + s2);
   }
   if (nr % 2)
   {
      int r1 = row[nr-1];
      int b1 = ii[r1];
      int s1 = ii[r1+1]-b1;

      copy (s1, jj+b1, q+im);
      i_merged[++nnr] = (im += s1);
   }
   t = p; p = q; q = t;

   while (nnr > 1)
   {
      im = 0;
      for (int o = 1; o < nnr; o += 2)
      {
         int b1 = im;
         int b2 = i_merged[o];
         im = i_merged[o+1];
         int s1 = b2-b1;
         int s2 = im-b2;

         merge (s1, s2, p+b1, p+b2, q+b1);

         i_merged[o/2+1] = im;
      }
      if (nnr % 2)
      {
         int b1 = im;
         im = i_merged[nnr];
         int s1 = im-b1;

         copy (s1, p+b1, q+b1);

         i_merged[(nnr+1)/2] = im;
      }
      t = p; p = q; q = t;
      nnr = (nnr+1)/2;
   }

#ifdef DEBUG
   if (!sorted_ok (im, merged))
      cerr << "mergecsrrows (...) !sorted_ok(...)" << endl;
#endif
}
