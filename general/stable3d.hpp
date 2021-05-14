
#include "mem_alloc.hpp"

class STable3DNode
{
   public:
      STable3DNode *Prev;
      int Column, Floor, Number;
};


class STable3D
{
   private:
      int Size, NElem;
      STable3DNode **Rows;
      MemAlloc <STable3DNode, 1024> NodesMem;

   public:
      STable3D (int nr);

      int Push (int r, int c, int f);

      int operator() (int r, int c, int f) const;

      int Push4 (int r, int c, int f, int t);

      int operator() (int r, int c, int f, int t) const;

      int NumberOfElements() { return NElem; };

      ~STable3D ();
};


class ITableNode
{
   public:
      ITableNode *Prev;
      int OtherDomain, InterfaceNumber, NFaces, IntWeight;
};


class ITable
{
   private:
      int NDomains, NInterfaces, TotWeight, *NumNeighbors;
      ITableNode **Domain, *WalkNode;
      MemAlloc <ITableNode, 1024> NodesMem;

   public:
      ITable (int ND);

      int ConnectDomains (int dom1, int dom2, int weight);

      int InterfaceNo (int dom1, int dom2);

      void GetInterfaceOffsetsAndNFaces (int *offsets, int *interfaces);

      int NumberOfInterfaces() { return NInterfaces; };

      int NumberOfNeighbors (int d) { return NumNeighbors[d]; };

      int TotalWeight() { return TotWeight; };

      void NbWalkStart (int dom) { WalkNode = Domain[dom]; };

      int NbWalkValid() { return (WalkNode != NULL); };

      void NbWalkNext() { WalkNode = WalkNode -> Prev; };

      int NbWalkOtherDomain() { return WalkNode -> OtherDomain; };

      int NbWalkInterfaceNum() { return WalkNode -> InterfaceNumber; };

      int NbWalkNFaces() { return WalkNode -> NFaces; };

      ~ITable();
};
