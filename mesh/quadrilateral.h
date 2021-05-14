#ifndef FILE_QUADRILATERAL
#define FILE_QUADRILATERAL
#include "mesh_header.h"

class Quadrilateral : public Element
{
protected:
    int indices[4];
    
public:
    Quadrilateral() : Element() {};
    Quadrilateral(const int *ind, int bg);
    virtual void SetVertices(const int *ind);

    /// Get the element's geometry
    virtual int GetGeometryType() const { return Element::QUADRILATERAL; }

    virtual void GetVertices(Array<int> &v) const;

    virtual int GetNVertices() const {return 4; };

    virtual void GetEdges(Array<int> &e) const { ; }
    
    virtual ~Quadrilateral();
};

#endif
