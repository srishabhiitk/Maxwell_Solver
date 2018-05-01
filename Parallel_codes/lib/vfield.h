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
      bool isPlanar;

      vfield(const grid &gridData, bool isFaceCentered_);

      void curl_3d(vfield *curl);
      void curl_3d_adv(vfield *curl1,vfield *curl2);
      void curl_planar(vfield *curl);
      void curl_planar_adv(vfield *curl1, vfield *curl2);


};

#endif
