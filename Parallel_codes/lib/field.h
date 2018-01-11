#ifndef FIELD_H
#define FIELD_H

#include <blitz/array.h>
#include "mpidata.h"
#include "grid.h"

class field {

    public:
      blitz::Array<double, 3> F;
      bool xStag;
      bool yStag;
      bool zStag;
      int local_Nx;
      int local_Ny;
      int local_Nz;
      const grid &gridData;


      field(const grid &gridData, bool xStag_, bool yStag_, bool zStag_);

      mpidata *mpiHandle;

};

#endif
