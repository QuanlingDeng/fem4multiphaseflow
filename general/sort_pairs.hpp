
template <class A, class B>
class Pair
{
   public:
      A one;
      B two;
};


template <class A, class B>
int ComparePairs (const void *_p, const void *_q);


template <class A, class B>
void SortPairs (Pair<A, B> *pairs, int size);
