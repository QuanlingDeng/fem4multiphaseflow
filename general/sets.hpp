
class IntegerSet
{
   private:
      Array<int> me;

   public:
      IntegerSet() { };

      IntegerSet (IntegerSet &s);

      IntegerSet (const int n, const int *p) { Recreate (n, p); };

      int Size() { return me.Size(); };

      operator Array<int>& () { return me; } ;

      int PickElement() { return me[0]; };

      int operator== (IntegerSet &s);

      void Recreate (const int n, const int *p);
};


class ListOfIntegerSets
{
   private:
      Array<IntegerSet *> TheList;

   public:

      int Size() { return TheList.Size(); };

      int PickElementInSet (int i) { return TheList[i] -> PickElement(); };

      int Insert (IntegerSet &s);

      int Lookup (IntegerSet &s);

      void AsTable (Table & t);

      ~ListOfIntegerSets();
};
