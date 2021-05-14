#ifndef FILE_FEM_HEADER
#define FILE_FEM_HEADER

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cstring>

#include "../linalg/linalg_header.h"
#include "../mesh/mesh_header.h"
#include "../general/general_header.h"
//#include <gsl/gsl_bspline.h>
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_cblas.h>
//#include <gsl/gsl_blas.h>

#include "Data.h"
#include "fem.h"
#include "postprocessing.h"
#include "basisfunctions.h"

#include "linearfem1d.h"
#include "quadraticfem1d.h"
#include "cubicfem1d.h"

#include "linearfem2d.h"
#include "quadraticfem2d.h"
#include "cubicfem2d.h"
#include "bilinearfem2d.h"
#include "linearfemelast2d.h"
#include "bilinearfemelast2d.h"

#include "linearfvem2d.h"
#include "quadraticfvem2d.h"
#include "cubicfvem2d.h"
//#include "quadraticbsplinefem2d.h"
#include "linearfemgeomech2d.h"

#include "linearfem3d.h"

//#include "stokes2dMINI.h"
//#include "stokes2dP1P0.h"
//#include "stokes2dCR.h"

#endif
