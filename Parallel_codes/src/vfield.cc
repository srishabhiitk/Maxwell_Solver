#include <iostream>
#include "vfield.h"

vfield::vfield(const grid &gridData, bool isFaceCentered_): gridData(gridData) {
  isFaceCentered=isFaceCentered_;

  if (isFaceCentered_){
    Vx = new field(gridData, false, true, true);
    Vy = new field(gridData, true, false, true);
    Vz = new field(gridData, true, true, false);
  }
  else{
    Vx = new field(gridData, true, false, false);
    Vy = new field(gridData, false, true, false);
    Vz = new field(gridData, false, false, true);
  }


}


void vfield::curl_3d(vfield *curl){
    bool isFaceCentered_c = curl->isFaceCentered;
    if (isFaceCentered!=isFaceCentered_c){
        double dx = gridData.inputData.dx;
        double dy = gridData.inputData.dy;
        double dz = gridData.inputData.dz;

        int Nx = gridData.local_colloq_x;
        int Ny = gridData.local_colloq_y;
        int Nz = gridData.local_colloq_z;

        int i,j,k;
        if (!isFaceCentered){
            for (i=0;i<Nx;i++){
                for (j=0;j<(Ny-1);j++){
                    for (k=0;k<(Nz-1);k++){
                      curl->Vx->F(i,j,k) = (Vz->F(i,j+1,k)-Vz->F(i,j,k))/dy - (Vy->F(i,j,k+1)-Vy->F(i,j,k))/dz;
                    }
                }
            }
            for (i=0;i<(Nx-1);i++){
                for (j=0;j<Ny;j++){
                    for (k=0;k<(Nz-1);k++){
                      curl->Vy->F(i,j,k) = (Vx->F(i,j,k+1)-Vx->F(i,j,k))/dz - (Vz->F(i+1,j,k)-Vz->F(i,j,k))/dx;
                    }
                }
            }
            for (i=0;i<(Nx-1);i++){
                for (j=0;j<(Ny-1);j++){
                    for (k=0;k<Nz;k++){
                      curl->Vz->F(i,j,k) = (Vy->F(i+1,j,k)-Vy->F(i,j,k))/dx - (Vx->F(i,j+1,k)-Vx->F(i,j,k))/dy;
                    }
                }
            }
        }
        else{
            for (i=0;i<(Nx-1);i++){
                for (j=1;j<(Ny-1);j++){
                    for (k=1;k<(Nz-1);k++){
                      curl->Vx->F(i,j,k) = (Vz->F(i,j,k)-Vz->F(i,j-1,k))/dy - (Vy->F(i,j,k)-Vy->F(i,j,k-1))/dz;
                    }
                }
            }
            for (i=1;i<(Nx-1);i++){
                for (j=0;j<(Ny-1);j++){
                    for (k=1;k<(Nz-1);k++){
                      curl->Vy->F(i,j,k) = (Vx->F(i,j,k)-Vx->F(i,j,k-1))/dz - (Vz->F(i,j,k)-Vz->F(i-1,j,k))/dx;
                    }
                }
            }
            for (i=1;i<(Nx-1);i++){
                for (j=1;j<(Ny-1);j++){
                    for (k=0;k<(Nz-1);k++){
                      curl->Vz->F(i,j,k) = (Vy->F(i,j,k)-Vy->F(i-1,j,k))/dx - (Vx->F(i,j,k)-Vx->F(i,j-1,k))/dy;
                    }
                }
            }
        }
    }
    else{
    std::cout<<"First two parameter should be staggered and unstaggered."<<std::endl;
    }

}
