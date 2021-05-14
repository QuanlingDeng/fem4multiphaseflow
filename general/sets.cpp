#include "array.hpp"
#include "table.hpp"
#include "sets.hpp"

using namespace std;

IntegerSet::IntegerSet (IntegerSet &s)
   : me (s.me.Size())
{
   for (int i = 0; i < me.Size(); i++)
      me[i] = s.me[i];
}

int IntegerSet::operator== (IntegerSet &s)
{
   if (me.Size() != s.me.Size())
      return 0;

   for (int i = 0; i < me.Size(); i++)
      if (me[i] != s.me[i])
         return 0;

   return 1;
}

void IntegerSet::Recreate (const int n, const int *p)
{
   int i, j;

   me.SetSize (n);

   for (i = 0; i < n; i++)
      me[i] = p[i];

   me.Sort();

   for (j = 0, i = 1; i < n; i++)
      if (me[i] != me[j])
         me[++j] = me[i];

   me.SetSize (j+1);
}


int ListOfIntegerSets::Insert (IntegerSet &s)
{
   for (int i = 0; i < TheList.Size(); i++)
      if (*TheList[i] == s)
         return i;

   TheList.Append (new IntegerSet (s));

   return TheList.Size()-1;
}

int ListOfIntegerSets::Lookup (IntegerSet &s)
{
   for (int i = 0; i < TheList.Size(); i++)
      if (*TheList[i] == s)
         return i;

   cerr << "ListOfIntegerSets::Lookup ()" << endl;
   return *((int *)NULL); // crash
}

void ListOfIntegerSets::AsTable (Table & t)
{
  int i;

  t.MakeI(Size());

  for (i = 0; i < Size(); i++)
    t.AddColumnsInRow (i, TheList[i] -> Size());

  t.MakeJ();

  for (i = 0; i < Size(); i++)
    t.SetRow (i, *TheList[i]);
}

ListOfIntegerSets::~ListOfIntegerSets()
{
   for (int i = 0; i < TheList.Size(); i++)
      delete TheList[i];
}
