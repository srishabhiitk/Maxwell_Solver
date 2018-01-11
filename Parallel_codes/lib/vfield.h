#ifndef VFIELD_H
#define VFIELD_H

#include <blitz/array.h>
#include "field.h"
#include "grid.h"

class vfield {

    public:
      field *Vx;
      field *Vy;
      field *Vz;

      const grid &gridData;

      bool isFaceCentered;

      vfield(const grid &gridData, bool isFaceCentered_);

      void curl_3d(vfield *curl);


};

#endif
