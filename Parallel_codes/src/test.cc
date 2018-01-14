#include <iostream>
#include <math.h>
#include <blitz/array.h>
#include "parallel.h"
#include "reader.h"
#include "grid.h"
#include "vfield.h"
#include "writer.h"
#include <sstream>



bool check_global_limits(blitz::TinyVector<int,3> v, grid gridData, bool xstag, bool ystag){
    int Nx_e = gridData.colloq_end_index_x;
    int Ny_e = gridData.colloq_end_index_y;
    if (xstag){
        Nx_e = Nx_e - 1;
    }
    if (ystag){
        Ny_e = Ny_e - 1;
    }
    if ((v(0)-gridData.colloq_start_index_x)>=0 && (v(1)-gridData.colloq_start_index_y)>=0){
        if ((v(0)-Nx_e)<=0 && (v(1)-Ny_e)<=0){
            return true;
        }
    }
    return false;
}

blitz::TinyVector<int,3> global_to_local(blitz::TinyVector<int,3> glob, grid gridData){
    blitz::TinyVector<int,3> loc = glob;
    loc(0) = glob(0) - gridData.colloq_start_index_x;
    loc(1) = glob(1) - gridData.colloq_start_index_y;
    return loc;
}



int main() {
    // INITIALIZE MPI
    MPI_Init(NULL, NULL);

    // ALL PROCESSES READ THE INPUT PARAMETERS
    reader inputData;

    // INITIALIZE PARALLELIZATION DATA
    parallel mpi(inputData);

    grid gridData(inputData, mpi);


    int num_timesteps = inputData.num_timesteps;
    double epsilon0 = 8.85*pow(10,-12);
    double mew0 = M_PI*4*pow(10,-7);
    double S = 1/pow(3,0.5);//inputData.S;
    double c = 1/pow(epsilon0*mew0,0.5);
    double dt=S*inputData.dx/c;
    double Ca = 1.0, Da = 1.0;
    double Cb = dt/epsilon0;
    double Db = dt/mew0;
    double sigma=pow(10,0);
    double init=20.0;
    int t = 0; //timestep
    int i,j,k; //loop vars
    int Nx = gridData.local_colloq_x;
    int Ny = gridData.local_colloq_y;
    int Nz = gridData.local_colloq_z;
    double local_max, max = 0;


    vfield *E = new vfield(gridData,false);
    vfield *H = new vfield(gridData,true);

    vfield *Ecurl = new vfield(gridData,true);
    vfield *Hcurl = new vfield(gridData,false);

    vfield *Etotal = new vfield(gridData,false);

    blitz::TinyVector<int,3> source;
    source = int(inputData.Nx/2), int(inputData.Ny/2), int(inputData.Nz/2);

    writer Ez_writer("Ez_data.h", gridData, E->Vz);
    // Ez_writer.closeWriter();


    for (t=0; t<num_timesteps; t++){

        if (t==0){
            if (check_global_limits(source, gridData, false, false)){
                E->Vz->F(global_to_local(source, gridData))=exp(-pow(-init+t,2)/sigma);
            }
        }


        MPI_Reduce(&local_max, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        // std::cout<<"Processor: "<<mpi.rank<<", max: "<<local_max<<std::endl;
        // MPI_Barrier( MPI_COMM_WORLD );
        if (mpi.rank==0){
            std::cout<<"t: "<<t<<", max: "<<max<<std::endl;
            // std::cout<<std::endl;
        }


        E->curl_3d(Ecurl);

        for (i=0;i<Nx;i++){
          for (j=0;j<(Ny-1);j++){
            for (k=0;k<(Nz-1);k++){
              H->Vx->F(i,j,k)=Da*H->Vx->F(i,j,k)-Db*Ecurl->Vx->F(i,j,k);
            }
          }
        }

        for (i=0;i<(Nx-1);i++){
          for (j=0;j<Ny;j++){
            for (k=0;k<(Nz-1);k++){
              H->Vy->F(i,j,k)=Da*H->Vy->F(i,j,k)-Db*Ecurl->Vy->F(i,j,k);
            }
          }
        }

        for (i=0;i<(Nx-1);i++){
          for (j=0;j<(Ny-1);j++){
            for (k=0;k<Nz;k++){
              H->Vz->F(i,j,k)=Da*H->Vz->F(i,j,k)-Db*Ecurl->Vz->F(i,j,k);
            }
          }
        }

        H->curl_3d(Hcurl);

        for (i=0;i<(Nx-1);i++){
          for (j=0;j<Ny;j++){
            for (k=0;k<Nz;k++){
              E->Vx->F(i,j,k)=Ca*E->Vx->F(i,j,k)+Cb*Hcurl->Vx->F(i,j,k);
            }
          }
        }

        for (i=0;i<Nx;i++){
          for (j=0;j<(Ny-1);j++){
            for (k=0;k<Nz;k++){
              E->Vy->F(i,j,k)=Ca*E->Vy->F(i,j,k)+Cb*Hcurl->Vy->F(i,j,k);
            }
          }
        }

        for (i=0;i<Nx;i++){
          for (j=0;j<Ny;j++){
            for (k=0;k<(Nz-1);k++){
              E->Vz->F(i,j,k)=Ca*E->Vz->F(i,j,k)+Cb*Hcurl->Vz->F(i,j,k);
            }
          }
        }




        for (i=0;i<(Nx-1);i++){
          for (j=0;j<(Ny-1);j++){
            for (k=0;k<(Nz-1);k++){
              Etotal->Vz->F(i,j,k)=E->Vz->F(i,j,k)*E->Vz->F(i,j,k)+E->Vy->F(i,j,k)*E->Vy->F(i,j,k)+E->Vx->F(i,j,k)*E->Vx->F(i,j,k);
            }
          }
        }

        E->Vz->mpiHandle->syncData();
        E->Vy->mpiHandle->syncData();
        E->Vx->mpiHandle->syncData();


        if (check_global_limits(source, gridData, false, false)){
            E->Vz->F(global_to_local(source, gridData))=exp(-pow(-init+t,2)/sigma);
        }

        local_max = blitz::max(Etotal->Vz->F);


        if (t>=20&&t<=200&&(t%10==0)){
            Ez_writer.writeHDF5(t);
        }


    }

    Ez_writer.closeWriter();

    MPI_Finalize();

    return 0;
}
