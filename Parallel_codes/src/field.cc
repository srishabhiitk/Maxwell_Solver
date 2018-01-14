#include "field.h"


field::field(const grid &gridData, bool xStag_, bool yStag_, bool zStag_): gridData(gridData){
  local_Nx=gridData.local_colloq_x;
  local_Ny=gridData.local_colloq_y;
  local_Nz=gridData.local_colloq_z;
  xStag=xStag_;
  yStag=yStag_;
  zStag=zStag_;


  if (xStag_){
    local_Nx=local_Nx-1;
  }
  if (yStag_){
    local_Ny=local_Ny-1;
  }
  if (zStag_){
    local_Nz=local_Nz-1;
  }

  F.resize(local_Nx,local_Ny,local_Nz);
  F=0.0;

  mpiHandle = new mpidata(F, xStag, yStag, gridData.rankData, gridData);
  mpiHandle->createSubarrays(xStag, yStag, zStag);

}
