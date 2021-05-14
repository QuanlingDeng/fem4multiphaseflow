#ifndef FILE_MEM_ALLOC
#define FILE_MEM_ALLOC


template <class Elem, int Num>
class StackPart
{
   public:
      StackPart<Elem, Num> *Prev;
      Elem Elements[Num];
};

template <class Elem, int Num>
class Stack
{
   private:
      StackPart <Elem, Num> *TopPart, *TopFreePart;
      int UsedInTop, SSize;
   public:
      Stack() { TopPart = TopFreePart = NULL; UsedInTop = Num; SSize = 0; };
      int Size() { return SSize; };
      void Push (Elem E);
      Elem Pop();
      void Clear();
      ~Stack() { Clear(); };
};

template <class Elem, int Num>
void Stack <Elem, Num>::Push (Elem E)
{
   StackPart <Elem, Num> *aux;
   if (UsedInTop == Num)
   {
      if (TopFreePart == NULL)
         aux = new StackPart <Elem, Num>;
      else
         TopFreePart = (aux = TopFreePart)->Prev;
      aux->Prev = TopPart;
      TopPart = aux;
      UsedInTop = 0;
   }
   TopPart->Elements[UsedInTop++] = E;
   SSize++;
}

template <class Elem, int Num>
Elem Stack <Elem, Num>::Pop()
{
   StackPart <Elem, Num> *aux;
   if (UsedInTop == 0)
   {
      TopPart = (aux = TopPart)->Prev;
      aux->Prev = TopFreePart;
      TopFreePart = aux;
      UsedInTop = Num;
   }
   SSize--;
   return TopPart->Elements[--UsedInTop];
}

template <class Elem, int Num>
void Stack <Elem, Num>::Clear()
{
   StackPart <Elem, Num> *aux;
   while (TopPart != NULL)
   {
      TopPart = (aux = TopPart)->Prev;
      delete aux;
   }
   while (TopFreePart != NULL)
   {
      TopFreePart = (aux = TopFreePart)->Prev;
      delete aux;
   }
   UsedInTop = Num;
   SSize = 0;
}


template <class Elem, int Num>
class MemAllocNode
{
   public:
      MemAllocNode <Elem, Num> *Prev;
      Elem Elements[Num];
};

template <class Elem, int Num>
class MemAlloc
{
   private:
      MemAllocNode <Elem, Num> *Last;
      int AllocatedInLast;
      Stack <Elem *, Num> UsedMem;
   public:
      MemAlloc() { Last = NULL; AllocatedInLast = Num; };
      Elem *Alloc();
      void Free (Elem *);
      void Clear();
      ~MemAlloc() { Clear(); };
};

template <class Elem, int Num>
Elem *MemAlloc <Elem, Num>::Alloc()
{
   MemAllocNode <Elem, Num> *aux;
   if (UsedMem.Size() > 0)
      return UsedMem.Pop();
   if (AllocatedInLast == Num)
   {
      aux = Last;
      Last = new MemAllocNode <Elem, Num>;
      Last->Prev = aux;
      AllocatedInLast = 0;
   }
   return &(Last->Elements[AllocatedInLast++]);
}

template <class Elem, int Num>
void MemAlloc <Elem, Num>::Free (Elem *E)
{
   UsedMem.Push (E);
}

template <class Elem, int Num>
void MemAlloc <Elem, Num>::Clear()
{
   MemAllocNode <Elem, Num> *aux;
   while (Last != NULL)
   {
      aux = Last->Prev;
      delete Last;
      Last = aux;
   }
   AllocatedInLast = Num;
   UsedMem.Clear();
}


#endif
