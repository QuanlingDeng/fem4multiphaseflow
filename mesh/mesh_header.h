#ifndef FILE_MESH_HEADER
#define FILE_MESH_HEADER

#include <iostream>
using namespace std;
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include "../linalg/linalg_header.h"
#include "../general/table.hpp"
#include "../geometry/geometry.h"

#include "vertex.h"
#include "element.h"
#include "point.h"
#include "segment.h"
#include "triangle.h"
#include "quadrilateral.h"
#include "tetrahedral.h"
#include "mesh.h"
#include "controlvolume.h"
#include "dualmesh.h"


#include "TwoDNURBS.h"
#include "ThreeDNURBS.h"
#include "SurfaceNURBS.h"

#include "element_type.h"


#endif
