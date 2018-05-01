#include <iostream>
#include "vfield.h"

vfield::vfield(const grid &gridData, bool isFaceCentered_): gridData(gridData) {
    isFaceCentered = isFaceCentered_;
    isPlanar = gridData.isPlanar;

    if (isPlanar){
        if (isFaceCentered_){
            Vx = new field(gridData, false, false, true);
            Vy = NULL;
            Vz = new field(gridData, true, false, false);
        }
        else{
            Vx = NULL;
            Vy = new field(gridData, false, false, false);
            Vz = NULL;
        }
    }
    else{
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

}


void vfield::curl_3d(vfield *curl){
    int n_threads = gridData.inputData.n_threads;
    bool isFaceCentered_c = curl->isFaceCentered;
    if (isFaceCentered!=isFaceCentered_c){

        double dx = gridData.inputData.dx;
        double dy = gridData.inputData.dy;
        double dz = gridData.inputData.dz;

        int Nx = gridData.local_colloq_x;
        int Ny = gridData.local_colloq_y;
        int Nz = gridData.local_colloq_z;



        if (!isFaceCentered){
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int i=0;i<Nx;i++){
                for (int j=0;j<(Ny-1);j++){
                    for (int k=0;k<(Nz-1);k++){
                      curl->Vx->F(i,j,k) = (Vz->F(i,j+1,k)-Vz->F(i,j,k))/dy - (Vy->F(i,j,k+1)-Vy->F(i,j,k))/dz;
                    }
                }
            }
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int j=0;j<Ny;j++){
                for (int i=0;i<(Nx-1);i++){
                    for (int k=0;k<(Nz-1);k++){
                      curl->Vy->F(i,j,k) = (Vx->F(i,j,k+1)-Vx->F(i,j,k))/dz - (Vz->F(i+1,j,k)-Vz->F(i,j,k))/dx;
                    }
                }
            }
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int k=0;k<Nz;k++){
                for (int j=0;j<(Ny-1);j++){
                    for (int i=0;i<(Nx-1);i++){
                      curl->Vz->F(i,j,k) = (Vy->F(i+1,j,k)-Vy->F(i,j,k))/dx - (Vx->F(i,j+1,k)-Vx->F(i,j,k))/dy;
                    }
                }
            }
        }
        else{
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int i=0;i<(Nx-1);i++){
                for (int j=1;j<(Ny-1);j++){
                    for (int k=1;k<(Nz-1);k++){
                      curl->Vx->F(i,j,k) = (Vz->F(i,j,k)-Vz->F(i,j-1,k))/dy - (Vy->F(i,j,k)-Vy->F(i,j,k-1))/dz;
                    }
                }
            }
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int j=0;j<(Ny-1);j++){
                for (int i=1;i<(Nx-1);i++){
                    for (int k=1;k<(Nz-1);k++){
                      curl->Vy->F(i,j,k) = (Vx->F(i,j,k)-Vx->F(i,j,k-1))/dz - (Vz->F(i,j,k)-Vz->F(i-1,j,k))/dx;
                    }
                }
            }
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int k=0;k<(Nz-1);k++){
                for (int j=1;j<(Ny-1);j++){
                    for (int i=1;i<(Nx-1);i++){
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



void vfield::curl_3d_adv(vfield *curl1, vfield *curl2){
    int n_threads = gridData.inputData.n_threads;
    bool isFaceCentered_c1 = curl1->isFaceCentered;
    bool isFaceCentered_c2 = curl2->isFaceCentered;
    if (isFaceCentered!=isFaceCentered_c1 && isFaceCentered_c1==isFaceCentered_c2){

        double dx = gridData.inputData.dx;
        double dy = gridData.inputData.dy;
        double dz = gridData.inputData.dz;

        int Nx = gridData.local_colloq_x;
        int Ny = gridData.local_colloq_y;
        int Nz = gridData.local_colloq_z;



        if (!isFaceCentered){
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int i=0;i<Nx;i++){
                for (int j=0;j<(Ny-1);j++){
                    for (int k=0;k<(Nz-1);k++){
                      curl1->Vx->F(i,j,k) = (Vz->F(i,j+1,k)-Vz->F(i,j,k))/dy;
                      curl2->Vx->F(i,j,k) = (Vy->F(i,j,k+1)-Vy->F(i,j,k))/dz;
                    }
                }
            }
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int j=0;j<Ny;j++){
                for (int i=0;i<(Nx-1);i++){
                    for (int k=0;k<(Nz-1);k++){
                      curl1->Vy->F(i,j,k) = (Vx->F(i,j,k+1)-Vx->F(i,j,k))/dz;
                      curl2->Vy->F(i,j,k) = (Vz->F(i+1,j,k)-Vz->F(i,j,k))/dx;
                    }
                }
            }
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int k=0;k<Nz;k++){
                for (int j=0;j<(Ny-1);j++){
                    for (int i=0;i<(Nx-1);i++){
                      curl1->Vz->F(i,j,k) = (Vy->F(i+1,j,k)-Vy->F(i,j,k))/dx;
                      curl2->Vz->F(i,j,k) = (Vx->F(i,j+1,k)-Vx->F(i,j,k))/dy;
                    }
                }
            }
        }
        else{
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int i=0;i<(Nx-1);i++){
                for (int j=1;j<(Ny-1);j++){
                    for (int k=1;k<(Nz-1);k++){
                      curl1->Vx->F(i,j,k) = (Vz->F(i,j,k)-Vz->F(i,j-1,k))/dy;
                      curl2->Vx->F(i,j,k) = (Vy->F(i,j,k)-Vy->F(i,j,k-1))/dz;
                    }
                }
            }
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int j=0;j<(Ny-1);j++){
                for (int i=1;i<(Nx-1);i++){
                    for (int k=1;k<(Nz-1);k++){
                      curl1->Vy->F(i,j,k) = (Vx->F(i,j,k)-Vx->F(i,j,k-1))/dz;
                      curl2->Vy->F(i,j,k) = (Vz->F(i,j,k)-Vz->F(i-1,j,k))/dx;
                    }
                }
            }
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int k=0;k<(Nz-1);k++){
                for (int j=1;j<(Ny-1);j++){
                    for (int i=1;i<(Nx-1);i++){
                      curl1->Vz->F(i,j,k) = (Vy->F(i,j,k)-Vy->F(i-1,j,k))/dx;
                      curl2->Vz->F(i,j,k) = (Vx->F(i,j,k)-Vx->F(i,j-1,k))/dy;
                    }
                }
            }
        }
    }
    else{
    std::cout<<"Error: Passed parameters should both be either face-centered or not face-centered. The calling object should have opposite respective property!"<<std::endl;
    }

}






void vfield::curl_planar(vfield *curl){
    int n_threads = gridData.inputData.n_threads;

    double dx = gridData.inputData.dx;
    double dz = gridData.inputData.dz;

    int Nx = gridData.local_colloq_x;
    int Nz = gridData.local_colloq_z;


    if (curl->Vy==NULL && Vx==NULL && Vz==NULL){
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0; i<(Nx); i++){
            for (int j=0; j<(Nz-1); j++){
                curl->Vx->F(i,0,j) = -(Vy->F(i,0,j+1)-Vy->F(i,0,j))/dz;
            }
        }
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0; i<(Nx-1); i++){
            for (int j=0; j<(Nz); j++){
                curl->Vz->F(i,0,j) = (Vy->F(i+1,0,j)-Vy->F(i,0,j))/dx;
            }
        }
    }
    else if (Vy==NULL && curl->Vx==NULL && curl->Vz==NULL){
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=1; i<(Nx-1); i++){
            for (int j=1; j<(Nz-1); j++){
                curl->Vy->F(i,0,j) = (Vx->F(i,0,j)-Vx->F(i,0,j-1))/dz - (Vz->F(i,0,j)-Vz->F(i-1,0,j))/dx;
            }
        }
    }
    else{
        std::cout<<"The parameters passed to curl_planar function are not correctly defined for planar transverse fields!"<<std::endl;
    }
}


void vfield::curl_planar_adv(vfield *curl1, vfield *curl2){
    int n_threads = gridData.inputData.n_threads;

    double dx = gridData.inputData.dx;
    double dz = gridData.inputData.dz;

    int Nx = gridData.local_colloq_x;
    int Nz = gridData.local_colloq_z;


    if (curl1->Vy==NULL && curl2->Vy==NULL && Vx==NULL && Vz==NULL){
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0; i<(Nx); i++){
            for (int j=0; j<(Nz-1); j++){
                curl2->Vx->F(i,0,j) = (Vy->F(i,0,j+1)-Vy->F(i,0,j))/dz;
            }
        }
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0; i<(Nx-1); i++){
            for (int j=0; j<(Nz); j++){
                curl1->Vz->F(i,0,j) = (Vy->F(i+1,0,j)-Vy->F(i,0,j))/dx;
            }
        }
    }
    else if (Vy==NULL && curl1->Vx==NULL && curl1->Vz==NULL && curl2->Vx==NULL && curl2->Vz==NULL){
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=1; i<(Nx-1); i++){
            for (int j=1; j<(Nz-1); j++){
                curl1->Vy->F(i,0,j) = (Vx->F(i,0,j)-Vx->F(i,0,j-1))/dz;
                curl2->Vy->F(i,0,j) = (Vz->F(i,0,j)-Vz->F(i-1,0,j))/dx;
            }
        }
    }
    else{
        std::cout<<"The parameters passed to curl_planar function are not correctly defined for planar transverse fields!"<<std::endl;
    }
}
