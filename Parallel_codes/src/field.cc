#include "field.h"


field::field(const grid &gridData, bool xStag_, bool yStag_, bool zStag_): gridData(gridData){
  local_Nx=gridData.local_colloq_x;
  local_Ny=gridData.local_colloq_y;
  local_Nz=gridData.local_colloq_z;
  xStag=xStag_;
  yStag=yStag_;
  zStag=zStag_;
  int local_Nx_=local_Nx;
  int local_Ny_=local_Ny;
  int local_Nz_=local_Nz;


  if (xStag_){
    local_Nx_=local_Nx-1;
  }
  if (yStag_){
    local_Ny_=local_Ny-1;
  }
  if (zStag_){
    local_Nz_=local_Nz-1;
  }

  F.resize(local_Nx_,local_Ny_,local_Nz_);
  F=0.0;

  mpiHandle = new mpidata(F, xStag, yStag, gridData.rankData, gridData);
  mpiHandle->createSubarrays(xStag, yStag, zStag);

}
