#include "maxwell.h"



maxwell::maxwell(reader &_inputData, parallel &_mpi, grid &_gridData): inputData(_inputData), mpi(_mpi), gridData(_gridData) {
    n_threads = inputData.n_threads;
    num_timesteps = inputData.num_timesteps;
    epsilon0 = 8.85*pow(10,-12);
    mew0 = M_PI*4*pow(10,-7);
    S = 1/pow(3,0.5);//inputData.S;
    c = 1/pow(epsilon0*mew0,0.5);
    dt=S*inputData.dx/c;
    dx=inputData.dx;
    dy=inputData.dy;
    dz=inputData.dz;
    sigma=pow(10,-4);
    init=60;
    lambda=10*dx;
    usePW = inputData.usePW;
    if (usePW){
        // double theta = M_PI/2.0;
        // H_x_dir = sin(theta);
        // H_z_dir = -cos(theta);
        // E_y_dir = 1;
        // kx = cos(theta);
        // kz = sin(theta);
        // ky = 0;

        H_x_dir = -1.0/sqrt(6);
        H_y_dir = 2.0/sqrt(6);
        H_z_dir = -1.0/sqrt(6);
        E_x_dir = 1.0/sqrt(2);
        E_y_dir = 0;
        E_z_dir = -1.0/sqrt(2);
        kx = 1.0/sqrt(3);
        kz = 1.0/sqrt(3);
        ky = 1.0/sqrt(3);

        // H_x_dir = 0.0;
        // H_y_dir = 0.0;
        // H_z_dir = -1.0;
        // E_x_dir = 1.0;
        // E_y_dir = 0.0;
        // E_z_dir = 0.0;
        // kx = 0.0;
        // ky = 1.0;
        // kz = 0.0;

        p1 = 20;
        q1 = 20;
        r1 = 20;
        p2 = 80;
        q2 = 80;
        r2 = 80;
        // pw_start_cord = 300,0,300;
        // pw_end_cord = 500,0,500;
        plane_wave_initialise();
    }
}



void maxwell::solve(){
    if (gridData.isPlanar){
        solve_planar();
    }
    else{
        solve3d_pml();
    }
}



void maxwell::solve3d(){
    double Ca = 1.0, Da = 1.0;
    double Cb = dt/epsilon0;
    double Db = dt/mew0;
    // double sigma=pow(10,0);
    // double init=20.0;
    int t = 0; //timestep
    //int i,j,k; //loop vars
    int Nx = gridData.local_colloq_x;
    int Ny = gridData.local_colloq_y;
    int Nz = gridData.local_colloq_z;
    double local_max, max = 0;

    struct timeval begin,end;


    vfield *E = new vfield(gridData,false);
    vfield *H = new vfield(gridData,true);

    vfield *Ecurl = new vfield(gridData,true);
    vfield *Hcurl = new vfield(gridData,false);

    vfield *Etotal = new vfield(gridData,false);

    blitz::TinyVector<int,3> source;
    source = int(inputData.Nx/2), int(inputData.Ny/2), int(inputData.Nz/2);

    writer Ez_writer("Ez_data2.h5", false, gridData, Etotal->Vz);
    // Ez_writer.closeWriter();

    gettimeofday(&begin,NULL);
    for (t=0; t<num_timesteps; t++){

        // if (t==0){
        //     if (check_global_limits(source, false, false)){
        //         E->Vz->F(global_to_local(source))=exp(-pow(-init+t,2)/sigma);
        //     }
        // }


        MPI_Reduce(&local_max, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (mpi.rank==0){
            std::cout<<"t: "<<t<<", max: "<<max<<std::endl;
            // std::cout<<std::endl;
        }


        E->curl_3d(Ecurl);
        if (usePW){
            plane_wave_execute(Ecurl,t);
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<Nx;i++){
          for (int j=0;j<(Ny-1);j++){
            for (int k=0;k<(Nz-1);k++){
              H->Vx->F(i,j,k)=Da*H->Vx->F(i,j,k)-Db*Ecurl->Vx->F(i,j,k);
            }
          }
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<Ny;j++){
            for (int k=0;k<(Nz-1);k++){
              H->Vy->F(i,j,k)=Da*H->Vy->F(i,j,k)-Db*Ecurl->Vy->F(i,j,k);
            }
          }
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<(Ny-1);j++){
            for (int k=0;k<Nz;k++){
              H->Vz->F(i,j,k)=Da*H->Vz->F(i,j,k)-Db*Ecurl->Vz->F(i,j,k);
            }
          }
        }

        H->curl_3d(Hcurl);
        if (usePW){
            plane_wave_execute(Hcurl,t);
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<Ny;j++){
            for (int k=0;k<Nz;k++){
              E->Vx->F(i,j,k)=Ca*E->Vx->F(i,j,k)+Cb*Hcurl->Vx->F(i,j,k);
            }
          }
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<Nx;i++){
          for (int j=0;j<(Ny-1);j++){
            for (int k=0;k<Nz;k++){
              E->Vy->F(i,j,k)=Ca*E->Vy->F(i,j,k)+Cb*Hcurl->Vy->F(i,j,k);
            }
          }
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<Nx;i++){
          for (int j=0;j<Ny;j++){
            for (int k=0;k<(Nz-1);k++){
              E->Vz->F(i,j,k)=Ca*E->Vz->F(i,j,k)+Cb*Hcurl->Vz->F(i,j,k);
            }
          }
        }




        E->Vz->mpiHandle->syncData();
        E->Vy->mpiHandle->syncData();
        E->Vx->mpiHandle->syncData();



        //Has to be placed after syncing because sometimes in plane-wave simulation
        //curl gets updated by plane-wave where it is not updated (is at boundary)
        //it does not affect anything because things get reset after syncing but
        //curl specifically at that point keeps on getting edited! (NOT RESOLVED YET!)
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<(Ny-1);j++){
            for (int k=0;k<(Nz-1);k++){
              Etotal->Vz->F(i,j,k)=E->Vz->F(i,j,k)*E->Vz->F(i,j,k)+E->Vy->F(i,j,k)*E->Vy->F(i,j,k)+E->Vx->F(i,j,k)*E->Vx->F(i,j,k);
            }
          }
        }


        // if (check_global_limits(source, false, false)){
        //     E->Vz->F(global_to_local(source))=exp(-pow(-init+t,2)/sigma);
        // }

        local_max = blitz::max(Etotal->Vz->F);
        // std::cout<<"rank: "<<mpi.rank<<", max:"<<local_max<<std::endl;


        if (t%10==0){
            Ez_writer.writeHDF5(t);
        }

        // if (mpi.rank==0){
        //     std::cout<<"t="<<t<<std::endl;
        // }


    }
    gettimeofday(&end,NULL);
    std::cout<<"Time taken: "<<((end.tv_sec-begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.0e6<<" sec"<<std::endl;
    Ez_writer.closeWriter();
}



void maxwell::solve3d_pml(){
    double Ca = 1.0, Da = 1.0;
    // double sigma=pow(10,0);
    // double init=20.0;
    int t = 0;
    int Nx = gridData.local_colloq_x;
    int Ny = gridData.local_colloq_y;
    int Nz = gridData.local_colloq_z;
    double local_max, max = 0;
    double dx = gridData.inputData.dx;


    // For PML BC
    double a_max=0.0;
    double k_max=15;
    double m_pml=3;
    double ma_pml=1;
    double sigma_optimal=0.8*(m_pml+1)/dx/pow(mew0/epsilon0,0.5);
    double sigma_max=sigma_optimal*0.65;
    double d_pml=10*dx;

    // double a_max=0;
    // double k_max=1;
    // double m_pml=3;
    // double ma_pml=1;
    // double sigma_optimal=0.8*(m_pml+1)/dx/pow(mew0/epsilon0,0.5);
    // double sigma_max=sigma_optimal*0;
    // double d_pml=10*dx;


    blitz::Array<double,1> k_ex(Nx),k_ey(Ny),k_ez(Nz);
    blitz::Array<double,1> k_hx(Nx-1),k_hy(Ny-1),k_hz(Nz-1);

    k_fn(0,true,k_ex,d_pml,m_pml,k_max);
    k_fn(1,true,k_ey,d_pml,m_pml,k_max);
    k_fn(2,true,k_ez,d_pml,m_pml,k_max);
    k_fn(0,false,k_hx,d_pml,m_pml,k_max);
    k_fn(1,false,k_hy,d_pml,m_pml,k_max);
    k_fn(2,false,k_hz,d_pml,m_pml,k_max);



    blitz::Array<double,1> b_ex(Nx),b_ey(Ny),b_ez(Nz),c_ex(Nx),c_ey(Ny),c_ez(Nz);
    blitz::Array<double,1> b_hx(Nx-1),b_hy(Ny-1),b_hz(Nz-1),c_hx(Nx-1),c_hy(Ny-1),c_hz(Nz-1);

    b_fn(0,true,b_ex,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    b_fn(1,true,b_ey,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    b_fn(2,true,b_ez,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    c_fn(0,true,c_ex,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    c_fn(1,true,c_ey,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    c_fn(2,true,c_ez,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    b_fn(0,false,b_hx,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    b_fn(1,false,b_hy,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    b_fn(2,false,b_hz,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    c_fn(0,false,c_hx,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    c_fn(1,false,c_hy,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    c_fn(2,false,c_hz,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);

    // std::cout<<"rank:"<<mpi.rank<<std::endl<<c_ez<<std::endl;


    vfield *E = new vfield(gridData,false);
    vfield *Etotal = new vfield(gridData,false);
    vfield *H = new vfield(gridData,true);

    vfield *Cb = new vfield(gridData,false);
    vfield *Db = new vfield(gridData,true);

    vfield *Ecurl1 = new vfield(gridData,true);
    vfield *Ecurl2 = new vfield(gridData,true);
    vfield *Hcurl1 = new vfield(gridData,false);
    vfield *Hcurl2 = new vfield(gridData,false);

    vfield *psi_E1 = new vfield(gridData,false);
    vfield *psi_E2 = new vfield(gridData,false);
    vfield *psi_H1 = new vfield(gridData,true);
    vfield *psi_H2 = new vfield(gridData,true);

    Cb->Vx->F=dt/epsilon0;
    Cb->Vy->F=dt/epsilon0;
    Cb->Vz->F=dt/epsilon0;
    Db->Vx->F=dt/mew0;
    Db->Vy->F=dt/mew0;
    Db->Vz->F=dt/mew0;


    // //Solid block at center
    // blitz::TinyVector<int,3> loBound;
    // loBound=45,45,45;
    // loBound = global_to_local(loBound);
    // blitz::TinyVector<int,3> upBound;
    // upBound=55,55,55;
    // upBound = global_to_local(upBound);
    // blitz::RectDomain<3> rd;
    // rd=blitz::RectDomain<3>(loBound,upBound);
    // Cb->Vx->F(rd) = dt/(epsilon0*6);
    // Cb->Vy->F(rd) = dt/(epsilon0*6);
    // Cb->Vz->F(rd) = dt/(epsilon0*6);




    blitz::TinyVector<int,3> source;
    source = int(inputData.Nx/2), int(inputData.Ny/2), int(inputData.Nz/2);


    writer Ez_writer("Ez_data.h5", false, gridData, Etotal->Vz);



    for (t=0; t<num_timesteps; t++){
        // if (t==0){
        //     if (check_global_limits(source, false, true)){
        //         // E->Vz->F(global_to_local(source))=exp(-pow(-init+t,2)/sigma);
        //         E->Vz->F(global_to_local(source))=sin(t/10.0*M_PI);
        //     }
        // }

        MPI_Reduce(&local_max, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (mpi.rank==0){
            std::cout<<"t: "<<t<<", max: "<<max<<std::endl;
            // std::cout<<std::endl;
        }

        E->curl_3d_adv(Ecurl1,Ecurl2);
        if (usePW){
            plane_wave_execute(Ecurl1,Ecurl2,t);
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<Nx;i++){
          for (int j=0;j<(Ny-1);j++){
            for (int k=0;k<(Nz-1);k++){
              psi_H1->Vx->F(i,j,k)=b_hy(j)*psi_H1->Vx->F(i,j,k)+c_hy(j)*Ecurl1->Vx->F(i,j,k);
              psi_H2->Vx->F(i,j,k)=b_hz(k)*psi_H2->Vx->F(i,j,k)+c_hz(k)*Ecurl2->Vx->F(i,j,k);
              H->Vx->F(i,j,k)=Da*H->Vx->F(i,j,k)-Db->Vx->F(i,j,k)*(Ecurl1->Vx->F(i,j,k)/k_hy(j)-Ecurl2->Vx->F(i,j,k)/k_hz(k))-Db->Vx->F(i,j,k)*(psi_H1->Vx->F(i,j,k)-psi_H2->Vx->F(i,j,k));
            }
          }
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<Ny;j++){
            for (int k=0;k<(Nz-1);k++){
              psi_H1->Vy->F(i,j,k)=b_hz(k)*psi_H1->Vy->F(i,j,k)+c_hz(k)*Ecurl1->Vy->F(i,j,k);
              psi_H2->Vy->F(i,j,k)=b_hx(i)*psi_H2->Vy->F(i,j,k)+c_hx(i)*Ecurl2->Vy->F(i,j,k);
              H->Vy->F(i,j,k)=Da*H->Vy->F(i,j,k)-Db->Vy->F(i,j,k)*(Ecurl1->Vy->F(i,j,k)/k_hz(k)-Ecurl2->Vy->F(i,j,k)/k_hx(i))-Db->Vy->F(i,j,k)*(psi_H1->Vy->F(i,j,k)-psi_H2->Vy->F(i,j,k));
            }
          }
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<(Ny-1);j++){
            for (int k=0;k<Nz;k++){
              psi_H1->Vz->F(i,j,k)=b_hx(i)*psi_H1->Vz->F(i,j,k)+c_hx(i)*Ecurl1->Vz->F(i,j,k);
              psi_H2->Vz->F(i,j,k)=b_hy(j)*psi_H2->Vz->F(i,j,k)+c_hy(j)*Ecurl2->Vz->F(i,j,k);
              H->Vz->F(i,j,k)=Da*H->Vz->F(i,j,k)-Db->Vz->F(i,j,k)*(Ecurl1->Vz->F(i,j,k)/k_hx(i)-Ecurl2->Vz->F(i,j,k)/k_hy(j))-Db->Vz->F(i,j,k)*(psi_H1->Vz->F(i,j,k)-psi_H2->Vz->F(i,j,k));
            }
          }
        }




        H->curl_3d_adv(Hcurl1,Hcurl2);
        if (usePW){
            plane_wave_execute(Hcurl1,Hcurl2,t);
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<Ny;j++){
            for (int k=0;k<Nz;k++){
              psi_E1->Vx->F(i,j,k)=b_ey(j)*psi_E1->Vx->F(i,j,k)+c_ey(j)*Hcurl1->Vx->F(i,j,k);
              psi_E2->Vx->F(i,j,k)=b_ez(k)*psi_E2->Vx->F(i,j,k)+c_ez(k)*Hcurl2->Vx->F(i,j,k);
              E->Vx->F(i,j,k)=Ca*E->Vx->F(i,j,k)+Cb->Vx->F(i,j,k)*(Hcurl1->Vx->F(i,j,k)/k_ey(j)-Hcurl2->Vx->F(i,j,k)/k_ez(k))+Cb->Vx->F(i,j,k)*(psi_E1->Vx->F(i,j,k)-psi_E2->Vx->F(i,j,k));
            }
          }
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<Nx;i++){
          for (int j=0;j<(Ny-1);j++){
            for (int k=0;k<Nz;k++){
              psi_E1->Vy->F(i,j,k)=b_ez(k)*psi_E1->Vy->F(i,j,k)+c_ez(k)*Hcurl1->Vy->F(i,j,k);
              psi_E2->Vy->F(i,j,k)=b_ex(i)*psi_E2->Vy->F(i,j,k)+c_ex(i)*Hcurl2->Vy->F(i,j,k);
              E->Vy->F(i,j,k)=Ca*E->Vy->F(i,j,k)+Cb->Vy->F(i,j,k)*(Hcurl1->Vy->F(i,j,k)/k_ez(k)-Hcurl2->Vy->F(i,j,k)/k_ex(i))+Cb->Vy->F(i,j,k)*(psi_E1->Vy->F(i,j,k)-psi_E2->Vy->F(i,j,k));
            }
          }
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<Nx;i++){
          for (int j=0;j<Ny;j++){
            for (int k=0;k<(Nz-1);k++){
              psi_E1->Vz->F(i,j,k)=b_ex(i)*psi_E1->Vz->F(i,j,k)+c_ex(i)*Hcurl1->Vz->F(i,j,k);
              psi_E2->Vz->F(i,j,k)=b_ey(j)*psi_E2->Vz->F(i,j,k)+c_ey(j)*Hcurl2->Vz->F(i,j,k);
              E->Vz->F(i,j,k)=Ca*E->Vz->F(i,j,k)+Cb->Vz->F(i,j,k)*(Hcurl1->Vz->F(i,j,k)/k_ex(i)-Hcurl2->Vz->F(i,j,k)/k_ey(j))+Cb->Vz->F(i,j,k)*(psi_E1->Vz->F(i,j,k)-psi_E2->Vz->F(i,j,k));
            }
          }
        }




        E->Vz->mpiHandle->syncData();
        E->Vy->mpiHandle->syncData();
        E->Vx->mpiHandle->syncData();



        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<(Ny-1);j++){
            for (int k=0;k<(Nz-1);k++){
              Etotal->Vz->F(i,j,k)=E->Vz->F(i,j,k)*E->Vz->F(i,j,k)+E->Vy->F(i,j,k)*E->Vy->F(i,j,k)+E->Vx->F(i,j,k)*E->Vx->F(i,j,k);
            }
          }
        }




        // std::cout<<E->Vy->F(blitz::Range(49,54),0,blitz::Range(49,54))<<std::endl;


        // if (check_global_limits(source, false, true)){
        //   // E->Vy->F(global_to_local(source))=exp(-pow(-init+t,2)/sigma);
        //   E->Vz->F(global_to_local(source))=sin(t/10.0*M_PI);
        // }

        local_max = blitz::max(Etotal->Vz->F);


        if ((t%10==0)){
            Ez_writer.writeHDF5(t);
        }

    }//end of for loop

    Ez_writer.closeWriter();

}



void maxwell::solve_planar(){
    double Ca = 1.0, Da = 1.0;
    // double Cb = dt/epsilon0;
    double Db = dt/mew0;
    // double sigma=pow(10,0);
    // double init=20.0;
    int t = 0; //timestep
    //int i,j,k; //loop vars
    int Nx = gridData.local_colloq_x;
    int Nz = gridData.local_colloq_z;
    double local_max, max = 0;


    struct timeval begin,end;


    vfield *E = new vfield(gridData,false);
    vfield *Cb = new vfield(gridData,false);
    vfield *H = new vfield(gridData,true);

    vfield *Ecurl = new vfield(gridData,true);
    vfield *Hcurl = new vfield(gridData,false);

    blitz::TinyVector<int,3> source;
    // source = int(inputData.Nx/2), 0, int(inputData.Nz/2);
    source = 688, 0, 120;

    writer Ey_writer("Ey_data_2d.h5", false, gridData, E->Vy);

    writer Cb_reader("lens_distort.h5", true, gridData, Cb->Vy);


    Cb_reader.readHDF5("data");

    // Cb->Vy->F=(Cb->Vy->F*1.5+1)*epsilon0;
    Cb->Vy->F=dt/epsilon0;

    // std::cout<<"Db:"<<Db<<std::endl;
    // Ez_writer.closeWriter();

    gettimeofday(&begin,NULL);
    for (t=0; t<num_timesteps; t++){
        //
        // if (t==0){
        //     if (check_global_limits(source, false, true)){
        //         E->Vy->F(global_to_local(source))=sin(t/20.0*M_PI);//exp(-pow(-init+t,2)/sigma);
        //     }
        // }

        // std::cout<<"Before max: "<<blitz::max(Ecurl->Vx->F*Ecurl->Vx->F)<<std::endl;
        E->curl_planar(Ecurl);
        // std::cout<<"After max: "<<blitz::max(Ecurl->Vx->F*Ecurl->Vx->F)<<std::endl;
        if (usePW){
            plane_wave_execute(Ecurl,t);
        }
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<Nx;i++){
          for (int j=0;j<(Nz-1);j++){
             H->Vx->F(i,0,j)=Da*H->Vx->F(i,0,j)-Db*Ecurl->Vx->F(i,0,j);
          }
        }
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<Nz;j++){
             H->Vz->F(i,0,j)=Da*H->Vz->F(i,0,j)-Db*Ecurl->Vz->F(i,0,j);
          }
        }

        H->curl_planar(Hcurl);
        if (usePW){
            plane_wave_execute(Hcurl,t);
        }
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<(Nz-1);j++){
             E->Vy->F(i,0,j)=Ca*E->Vy->F(i,0,j)+Cb->Vy->F(i,0,j)*Hcurl->Vy->F(i,0,j);
          }
        }





        E->Vy->mpiHandle->syncData();


        // if (check_global_limits(source, false, true)){
        //     E->Vy->F(global_to_local(source))=sin(t/20.0*M_PI);
        // }

        local_max = blitz::max(E->Vy->F*E->Vy->F);

        MPI_Reduce(&local_max, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (mpi.rank==0){
            std::cout<<"t: "<<t<<", max: "<<max<<std::endl;
            // std::cout<<std::endl;
        }


        if (t%5==0 && t>=100){
            Ey_writer.writeHDF5(t);
        }

    }



    gettimeofday(&end,NULL);
    if (mpi.rank==0){
        std::cout<<"Time taken: "<<((end.tv_sec-begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.0e6<<" sec"<<std::endl;
    }
    Ey_writer.closeWriter();
    Cb_reader.closeWriter();
}



void maxwell::solve_planar_pml(){
    double Ca = 1.0, Da = 1.0;
    double sigma=pow(10,0);
    double init=20.0;
    int t = 0;
    int Nx = gridData.local_colloq_x;
    int Nz = gridData.local_colloq_z;
    double local_max, max = 0;
    double dx = gridData.inputData.dx;


    // For PML BC
    // double a_max=0.0;
    // double k_max=15;
    // double m_pml=3;
    // double ma_pml=1;
    // double sigma_optimal=0.8*(m_pml+1)/dx/pow(mew0/epsilon0,0.5);
    // double sigma_max=sigma_optimal*0.65;
    // double d_pml=10*dx;

    double a_max=0;
    double k_max=1;
    double m_pml=3;
    double ma_pml=1;
    double sigma_optimal=0.8*(m_pml+1)/dx/pow(mew0/epsilon0,0.5);
    double sigma_max=sigma_optimal*0;
    double d_pml=10*dx;



    blitz::Array<double,1> k_ex(Nx),k_ez(Nz);
    blitz::Array<double,1> k_hx(Nx-1),k_hz(Nz-1);

    k_fn(0,true,k_ex,d_pml,m_pml,k_max);
    k_fn(2,true,k_ez,d_pml,m_pml,k_max);
    k_fn(0,false,k_hx,d_pml,m_pml,k_max);
    k_fn(2,false,k_hz,d_pml,m_pml,k_max);


    blitz::Array<double,1> b_ex(Nx),b_ez(Nz),c_ex(Nx),c_ez(Nz);
    blitz::Array<double,1> b_hx(Nx-1),b_hz(Nz-1),c_hx(Nx-1),c_hz(Nz-1);

    b_fn(0,true,b_ex,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    b_fn(2,true,b_ez,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    c_fn(0,true,c_ex,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    c_fn(2,true,c_ez,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    b_fn(0,false,b_hx,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    b_fn(2,false,b_hz,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    c_fn(0,false,c_hx,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
    c_fn(2,false,c_hz,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);


    vfield *E = new vfield(gridData,false);
    vfield *H = new vfield(gridData,true);

    vfield *Cb = new vfield(gridData,false);
    vfield *Db = new vfield(gridData,true);

    vfield *Ecurl1 = new vfield(gridData,true);
    vfield *Ecurl2 = new vfield(gridData,true);
    vfield *Hcurl1 = new vfield(gridData,false);
    vfield *Hcurl2 = new vfield(gridData,false);
    Ecurl1->Vx->F.resize(0,0,0);
    Ecurl2->Vz->F.resize(0,0,0);

    vfield *psi_E1 = new vfield(gridData,false);
    vfield *psi_E2 = new vfield(gridData,false);
    vfield *psi_H1 = new vfield(gridData,true);
    vfield *psi_H2 = new vfield(gridData,true);
    psi_H1->Vx->F.resize(0,0,0);
    psi_H2->Vz->F.resize(0,0,0);



    // writer Cb_reader("wg2.h5", true, gridData, Cb->Vy);
    //
    // Cb_reader.readHDF5("data");

    // Cb->Vy->F=(Cb->Vy->F*11 + 1)*epsilon0;
    // Cb->Vy->F=dt/Cb->Vy->F;
    Cb->Vy->F=dt/epsilon0;
    Db->Vx->F=dt/mew0;
    Db->Vz->F=dt/mew0;

    blitz::TinyVector<int,3> source;
    // source = int(inputData.Nx/2), 0, int(inputData.Nz/2);
    source = 688, 0, 120;//int(inputData.Nx/2), 0, 100;

    writer Ey_writer("Ey_data_2d.h5", false, gridData, E->Vy);


    // std::cout<<"Db:"<<Db->Vx->F(0,0,0)<<std::endl;
    // std::cout<<"Cb:"<<Cb->Vy->F(0,0,0)<<std::endl;


    for (t=0; t<num_timesteps; t++){
        // if (t==0){
        //     if (check_global_limits(source, false, true)){
        //         // E->Vy->F(global_to_local(source))=1;//exp(-pow(-init+t,2)/sigma);
        //         E->Vy->F(global_to_local(source))=sin(t/30.0*M_PI);
        //     }
        // }

        MPI_Reduce(&local_max, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (mpi.rank==0){
            std::cout<<"t: "<<t<<", max: "<<max<<std::endl;
            // std::cout<<std::endl;
        }

        E->curl_planar_adv(Ecurl1, Ecurl2);
        if (usePW){
            plane_wave_execute(Ecurl1, Ecurl2, t);
        }
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<Nx;i++){
            for (int j=0;j<(Nz-1);j++){
                psi_H2->Vx->F(i,0,j)=b_hz(j)*psi_H2->Vx->F(i,0,j)+c_hz(j)*Ecurl2->Vx->F(i,0,j);
                H->Vx->F(i,0,j)=Da*H->Vx->F(i,0,j)-Db->Vx->F(i,0,j)*(-Ecurl2->Vx->F(i,0,j)/k_hz(j)-psi_H2->Vx->F(i,0,j));
            }
        }
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
            for (int j=0;j<Nz;j++){
                psi_H1->Vz->F(i,0,j)=b_hx(i)*psi_H1->Vz->F(i,0,j)+c_hx(i)*Ecurl1->Vz->F(i,0,j);
                H->Vz->F(i,0,j)=Da*H->Vz->F(i,0,j)-Db->Vz->F(i,0,j)*(Ecurl1->Vz->F(i,0,j)/k_hx(i)+psi_H1->Vz->F(i,0,j));
            }
        }

        H->curl_planar_adv(Hcurl1, Hcurl2);
        if (usePW){
            plane_wave_execute(Hcurl1, Hcurl2, t);
        }
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
            for (int j=0;j<(Nz-1);j++){
                psi_E1->Vy->F(i,0,j)=b_ez(j)*psi_E1->Vy->F(i,0,j)+c_ez(j)*Hcurl1->Vy->F(i,0,j);
                psi_E2->Vy->F(i,0,j)=b_ex(i)*psi_E2->Vy->F(i,0,j)+c_ex(i)*Hcurl2->Vy->F(i,0,j);
                E->Vy->F(i,0,j)=Ca*E->Vy->F(i,0,j)+Cb->Vy->F(i,0,j)*(Hcurl1->Vy->F(i,0,j)/k_ez(j)-Hcurl2->Vy->F(i,0,j)/k_ex(i)+psi_E1->Vy->F(i,0,j)-psi_E2->Vy->F(i,0,j));
            }
        }




        E->Vy->mpiHandle->syncData();

        // std::cout<<E->Vy->F(blitz::Range(49,54),0,blitz::Range(49,54))<<std::endl;


        // if (check_global_limits(source, false, true)){
        //   // E->Vy->F(global_to_local(source))=1;//exp(-pow(-init+t,2)/sigma);
        //   E->Vy->F(global_to_local(source))=sin(t/30.0*M_PI);
        // }

        local_max = blitz::max(E->Vy->F*E->Vy->F);


        if (t%10==0){
          Ey_writer.writeHDF5(t);
        }

    }//end of for loop

    Ey_writer.closeWriter();
    // Cb_reader.closeWriter();
}



void maxwell::plane_wave_execute(vfield *curl, int timestep){
    blitz::TinyVector<int,3> temp;


    if (gridData.isPlanar){

        if (curl->isFaceCentered){
            // Changes on YZ plane
            if ((p1-1)>=gridData.colloq_start_index_x && (p1-1)<=gridData.colloq_end_index_x && start_z!=-1){
                for (int k=start_z; k<=end_z; k++){
                    temp = (p1-1),0,k;
                    temp = global_to_local(temp);
                    curl->Vz->F(temp) = curl->Vz->F(temp) + exp(-pow(((p1*dx*kx + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dx;
                }
            }
            if (p2>=gridData.colloq_start_index_x && p2<=gridData.colloq_end_index_x && start_z!=-1){
                for (int k=start_z; k<=end_z; k++){
                    temp = p2,0,k;
                    temp = global_to_local(temp);
                    curl->Vz->F(temp) = curl->Vz->F(temp) - exp(-pow(((p2*dx*kx + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dx;
                }
            }



            // Changes on XY plane
            if ((r1-1)>=gridData.colloq_start_index_z && (r1-1)<=gridData.colloq_end_index_z && start_x!=-1){
                for (int i=start_x; i<=end_x; i++){
                    temp = i,0,(r1-1);
                    temp = global_to_local(temp);
                    curl->Vx->F(temp) = curl->Vx->F(temp) - exp(-pow(((i*dx*kx + r1*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dz;
                }
            }
            if (r2>=gridData.colloq_start_index_z && r2<=gridData.colloq_end_index_z && start_x!=-1){
                for (int i=start_x; i<=end_x; i++){
                    temp = i,0,r2;
                    temp = global_to_local(temp);
                    curl->Vx->F(temp) = curl->Vx->F(temp) + exp(-pow(((i*dx*kx + r2*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dz;
                }
            }
        }
        else{
            // *************  if curl is not face-centered  *****************

            // Changes on YZ plane
            if (p1>=gridData.colloq_start_index_x && p1<=gridData.colloq_end_index_x && start_z!=-1){
                for (int k=start_z; k<=end_z; k++){
                    temp = p1,0,k;
                    temp = global_to_local(temp);
                    curl->Vy->F(temp) = curl->Vy->F(temp) + exp(-pow((((p1-0.5)*dx*kx + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dx/c/mew0;
                }
            }
            if (p2>=gridData.colloq_start_index_x && p2<=gridData.colloq_end_index_x && start_z!=-1){
                for (int k=start_z; k<=end_z; k++){
                    temp = p2,0,k;
                    temp = global_to_local(temp);
                    curl->Vy->F(temp) = curl->Vy->F(temp) - exp(-pow((((p2+0.5)*dx*kx + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dx/c/mew0;
                }
            }



            // Changes on XY plane
            if (r1>=gridData.colloq_start_index_z && r1<=gridData.colloq_end_index_z && start_x!=-1){
                for (int i=start_x; i<=end_x; i++){
                    temp = i,0,r1;
                    temp = global_to_local(temp);
                    curl->Vy->F(temp) = curl->Vy->F(temp) - exp(-pow(((i*dx*kx + (r1-0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dz/c/mew0;
                }
            }
            if (r2>=gridData.colloq_start_index_z && r2<=gridData.colloq_end_index_z && start_x!=-1){
                for (int i=start_x; i<=end_x; i++){
                    temp = i,0,r2;
                    temp = global_to_local(temp);
                    curl->Vy->F(temp) = curl->Vy->F(temp) + exp(-pow(((i*dx*kx + (r2+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dz/c/mew0;
                }
            }
        }



    }


    else{
        // **********************************************************************
        // ************************  For 3-D plane wave  ************************
        // **********************************************************************
        if (curl->isFaceCentered){
            // Changes on YZ plane
            if ((p1-1)>=gridData.colloq_start_index_x && (p1-1)<=gridData.colloq_end_index_x && start_y!=-1 && start_z!=-1){
                for (int j=start_y; j<=end_y-1; j++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = (p1-1),j,k;
                        temp = global_to_local(temp);
                        curl->Vz->F(temp) = curl->Vz->F(temp) - exp(-pow(((p1*dx*kx + (j+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dx;
                    }
                }
                for (int j=start_y; j<=end_y; j++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = (p1-1),j,k;
                        temp = global_to_local(temp);
                        curl->Vy->F(temp) = curl->Vy->F(temp) + exp(-pow(((p1*dx*kx + j*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_z_dir/dx;
                    }
                }
            }
            if (p2>=gridData.colloq_start_index_x && p2<=gridData.colloq_end_index_x && start_y!=-1 && start_z!=-1){
                for (int j=start_y; j<=end_y-1; j++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = p2,j,k;
                        temp = global_to_local(temp);
                        curl->Vz->F(temp) = curl->Vz->F(temp) + exp(-pow(((p2*dx*kx + (j+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dx;
                    }
                }
                for (int j=start_y; j<=end_y; j++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = p2,j,k;
                        temp = global_to_local(temp);
                        curl->Vy->F(temp) = curl->Vy->F(temp) - exp(-pow(((p2*dx*kx + j*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_z_dir/dx;
                    }
                }
            }


            // Changes on XZ plane
            if ((q1-1)>=gridData.colloq_start_index_y && (q1-1)<=gridData.colloq_end_index_y && start_x!=-1 && start_z!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = i,(q1-1),k;
                        temp = global_to_local(temp);
                        curl->Vz->F(temp) = curl->Vz->F(temp) + exp(-pow((((i+0.5)*dx*kx + q1*dy*ky + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_x_dir/dy;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = i,(q1-1),k;
                        temp = global_to_local(temp);
                        curl->Vx->F(temp) = curl->Vx->F(temp) - exp(-pow(((i*dx*kx + q1*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_z_dir/dy;
                    }
                }
            }
            if (q2>=gridData.colloq_start_index_y && q2<=gridData.colloq_end_index_y && start_x!=-1 && start_z!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = i,q2,k;
                        temp = global_to_local(temp);
                        curl->Vz->F(temp) = curl->Vz->F(temp) - exp(-pow((((i+0.5)*dx*kx + q2*dy*ky + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_x_dir/dy;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = i,q2,k;
                        temp = global_to_local(temp);
                        curl->Vx->F(temp) = curl->Vx->F(temp) + exp(-pow(((i*dx*kx + q2*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_z_dir/dy;
                    }
                }
            }



            // Changes on XY plane
            if ((r1-1)>=gridData.colloq_start_index_z && (r1-1)<=gridData.colloq_end_index_z && start_x!=-1 && start_y!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int j=start_y; j<=end_y; j++){
                        temp = i,j,(r1-1);
                        temp = global_to_local(temp);
                        curl->Vy->F(temp) = curl->Vy->F(temp) - exp(-pow((((i+0.5)*dx*kx + j*dy*ky + r1*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_x_dir/dz;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int j=start_y; j<=end_y-1; j++){
                        temp = i,j,(r1-1);
                        temp = global_to_local(temp);
                        curl->Vx->F(temp) = curl->Vx->F(temp) + exp(-pow(((i*dx*kx + (j+0.5)*dy*ky + r1*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dz;
                    }
                }
            }
            if (r2>=gridData.colloq_start_index_z && r2<=gridData.colloq_end_index_z && start_x!=-1 && start_y!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int j=start_y; j<=end_y; j++){
                        temp = i,j,r2;
                        temp = global_to_local(temp);
                        curl->Vy->F(temp) = curl->Vy->F(temp) + exp(-pow((((i+0.5)*dx*kx + j*dy*ky + r2*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_x_dir/dz;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int j=start_y; j<=end_y-1; j++){
                        temp = i,j,r2;
                        temp = global_to_local(temp);
                        curl->Vx->F(temp) = curl->Vx->F(temp) - exp(-pow(((i*dx*kx + (j+0.5)*dy*ky + r2*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dz;
                    }
                }
            }
        }
        else{
            // *************  if curl is not face-centered  *****************

            // Changes on YZ plane
            if (p1>=gridData.colloq_start_index_x && p1<=gridData.colloq_end_index_x && start_y!=-1 && start_z!=-1){
                for (int j=start_y; j<=end_y-1; j++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = p1,j,k;
                        temp = global_to_local(temp);
                        curl->Vy->F(temp) = curl->Vy->F(temp) + exp(-pow((((p1-0.5)*dx*kx + (j+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dx/c/mew0;
                    }
                }
                for (int j=start_y; j<=end_y; j++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = p1,j,k;
                        temp = global_to_local(temp);
                        curl->Vz->F(temp) = curl->Vz->F(temp) - exp(-pow((((p1-0.5)*dx*kx + j*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_y_dir/dx/c/mew0;
                    }
                }
            }
            if (p2>=gridData.colloq_start_index_x && p2<=gridData.colloq_end_index_x && start_y!=-1 && start_z!=-1){
                for (int j=start_y; j<=end_y-1; j++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = p2,j,k;
                        temp = global_to_local(temp);
                        curl->Vy->F(temp) = curl->Vy->F(temp) - exp(-pow((((p2+0.5)*dx*kx + (j+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dx/c/mew0;
                    }
                }
                for (int j=start_y; j<=end_y; j++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = p2,j,k;
                        temp = global_to_local(temp);
                        curl->Vz->F(temp) = curl->Vz->F(temp) + exp(-pow((((p2+0.5)*dx*kx + j*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_y_dir/dx/c/mew0;
                    }
                }
            }


            // Changes on XZ plane
            if (q1>=gridData.colloq_start_index_y && q1<=gridData.colloq_end_index_y && start_x!=-1 && start_z!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = i,q1,k;
                        temp = global_to_local(temp);
                        curl->Vx->F(temp) = curl->Vx->F(temp) - exp(-pow((((i+0.5)*dx*kx + (q1-0.5)*dy*ky + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dy/c/mew0;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = i,q1,k;
                        temp = global_to_local(temp);
                        curl->Vz->F(temp) = curl->Vz->F(temp) + exp(-pow(((i*dx*kx + (q1-0.5)*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dy/c/mew0;
                    }
                }
            }
            if (q2>=gridData.colloq_start_index_y && q2<=gridData.colloq_end_index_y && start_x!=-1 && start_z!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = i,q2,k;
                        temp = global_to_local(temp);
                        curl->Vx->F(temp) = curl->Vx->F(temp) + exp(-pow((((i+0.5)*dx*kx + (q2+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dy/c/mew0;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = i,q2,k;
                        temp = global_to_local(temp);
                        curl->Vz->F(temp) = curl->Vz->F(temp) - exp(-pow(((i*dx*kx + (q2+0.5)*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dy/c/mew0;
                    }
                }
            }



            // Changes on XY plane
            if (r1>=gridData.colloq_start_index_z && r1<=gridData.colloq_end_index_z && start_x!=-1 && start_y!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int j=start_y; j<=end_y; j++){
                        temp = i,j,r1;
                        temp = global_to_local(temp);
                        curl->Vx->F(temp) = curl->Vx->F(temp) + exp(-pow((((i+0.5)*dx*kx + j*dy*ky + (r1-0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_y_dir/dz/c/mew0;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int j=start_y; j<=end_y-1; j++){
                        temp = i,j,r1;
                        temp = global_to_local(temp);
                        curl->Vy->F(temp) = curl->Vy->F(temp) - exp(-pow(((i*dx*kx + (j+0.5)*dy*ky + (r1-0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dz/c/mew0;
                    }
                }
            }
            if (r2>=gridData.colloq_start_index_z && r2<=gridData.colloq_end_index_z && start_x!=-1 && start_y!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int j=start_y; j<=end_y; j++){
                        temp = i,j,r2;
                        temp = global_to_local(temp);
                        curl->Vx->F(temp) = curl->Vx->F(temp) - exp(-pow((((i+0.5)*dx*kx + j*dy*ky + (r2+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_y_dir/dz/c/mew0;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int j=start_y; j<=end_y-1; j++){
                        temp = i,j,r2;
                        temp = global_to_local(temp);
                        curl->Vy->F(temp) = curl->Vy->F(temp) + exp(-pow(((i*dx*kx + (j+0.5)*dy*ky + (r2+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dz/c/mew0;
                    }
                }
            }
        }
    }  //End of else
}



void maxwell::plane_wave_execute(vfield *curl1, vfield *curl2, int timestep){
    blitz::TinyVector<int,3> temp;


    if (gridData.isPlanar){ //For 2-D plane wave
        if (curl1->isFaceCentered){
            // Changes on YZ plane
            if ((p1-1)>=gridData.colloq_start_index_x && (p1-1)<=gridData.colloq_end_index_x && start_z!=-1){
                for (int k=start_z; k<=end_z; k++){
                    temp = (p1-1),0,k;
                    temp = global_to_local(temp);
                    curl1->Vz->F(temp) = curl1->Vz->F(temp) + exp(-pow(((p1*dx*kx + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dx;
                }

            }
            if (p2>=gridData.colloq_start_index_x && p2<=gridData.colloq_end_index_x && start_z!=-1){
                for (int k=start_z; k<=end_z; k++){
                    temp = p2,0,k;
                    temp = global_to_local(temp);
                    curl1->Vz->F(temp) = curl1->Vz->F(temp) - exp(-pow(((p2*dx*kx + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dx;
                }

            }



            // Changes on XY plane
            if ((r1-1)>=gridData.colloq_start_index_z && (r1-1)<=gridData.colloq_end_index_z && start_x!=-1){
                for (int i=start_x; i<=end_x; i++){
                    temp = i,0,(r1-1);
                    temp = global_to_local(temp);
                    curl2->Vx->F(temp) = curl2->Vx->F(temp) + exp(-pow(((i*dx*kx + r1*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dz;

                }
            }
            if (r2>=gridData.colloq_start_index_z && r2<=gridData.colloq_end_index_z && start_x!=-1){
                for (int i=start_x; i<=end_x; i++){
                    temp = i,0,r2;
                    temp = global_to_local(temp);
                    curl2->Vx->F(temp) = curl2->Vx->F(temp) - exp(-pow(((i*dx*kx + r2*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dz;
                }
            }
        }
        else{
            // *************  if curl is not face-centered  *****************

            // Changes on YZ plane
            if (p1>=gridData.colloq_start_index_x && p1<=gridData.colloq_end_index_x && start_z!=-1){
                for (int k=start_z; k<=end_z; k++){
                    temp = p1,0,k;
                    temp = global_to_local(temp);
                    curl2->Vy->F(temp) = curl2->Vy->F(temp) - exp(-pow((((p1-0.5)*dx*kx + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dx/c/mew0;
                }
            }
            if (p2>=gridData.colloq_start_index_x && p2<=gridData.colloq_end_index_x && start_z!=-1){
                for (int k=start_z; k<=end_z; k++){
                    temp = p2,0,k;
                    temp = global_to_local(temp);
                    curl2->Vy->F(temp) = curl2->Vy->F(temp) + exp(-pow((((p2+0.5)*dx*kx + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dx/c/mew0;
                }
            }



            // Changes on XY plane
            if (r1>=gridData.colloq_start_index_z && r1<=gridData.colloq_end_index_z && start_x!=-1){
                for (int i=start_x; i<=end_x; i++){
                    temp = i,0,r1;
                    temp = global_to_local(temp);
                    curl1->Vy->F(temp) = curl1->Vy->F(temp) - exp(-pow(((i*dx*kx + (r1-0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dz/c/mew0;
                }
            }
            if (r2>=gridData.colloq_start_index_z && r2<=gridData.colloq_end_index_z && start_x!=-1){
                for (int i=start_x; i<=end_x; i++){
                    temp = i,0,r2;
                    temp = global_to_local(temp);
                    curl1->Vy->F(temp) = curl1->Vy->F(temp) + exp(-pow(((i*dx*kx + (r2+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dz/c/mew0;
                }
            }
        }
    }
    else{
        // **********************************************************************
        // ************************  For 3-D plane wave  ************************
        // **********************************************************************
        if (curl1->isFaceCentered){
            // Changes on YZ plane
            if ((p1-1)>=gridData.colloq_start_index_x && (p1-1)<=gridData.colloq_end_index_x && start_y!=-1 && start_z!=-1){
                for (int j=start_y; j<=end_y-1; j++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = (p1-1),j,k;
                        temp = global_to_local(temp);
                        curl1->Vz->F(temp) = curl1->Vz->F(temp) - exp(-pow(((p1*dx*kx + (j+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dx;
                    }
                }
                for (int j=start_y; j<=end_y; j++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = (p1-1),j,k;
                        temp = global_to_local(temp);
                        curl2->Vy->F(temp) = curl2->Vy->F(temp) - exp(-pow(((p1*dx*kx + j*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_z_dir/dx;
                    }
                }
            }
            if (p2>=gridData.colloq_start_index_x && p2<=gridData.colloq_end_index_x && start_y!=-1 && start_z!=-1){
                for (int j=start_y; j<=end_y-1; j++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = p2,j,k;
                        temp = global_to_local(temp);
                        curl1->Vz->F(temp) = curl1->Vz->F(temp) + exp(-pow(((p2*dx*kx + (j+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dx;
                    }
                }
                for (int j=start_y; j<=end_y; j++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = p2,j,k;
                        temp = global_to_local(temp);
                        curl2->Vy->F(temp) = curl2->Vy->F(temp) + exp(-pow(((p2*dx*kx + j*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_z_dir/dx;
                    }
                }
            }


            // Changes on XZ plane
            if ((q1-1)>=gridData.colloq_start_index_y && (q1-1)<=gridData.colloq_end_index_y && start_x!=-1 && start_z!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = i,(q1-1),k;
                        temp = global_to_local(temp);
                        curl2->Vz->F(temp) = curl2->Vz->F(temp) - exp(-pow((((i+0.5)*dx*kx + q1*dy*ky + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_x_dir/dy;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = i,(q1-1),k;
                        temp = global_to_local(temp);
                        curl1->Vx->F(temp) = curl1->Vx->F(temp) - exp(-pow(((i*dx*kx + q1*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_z_dir/dy;
                    }
                }
            }
            if (q2>=gridData.colloq_start_index_y && q2<=gridData.colloq_end_index_y && start_x!=-1 && start_z!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = i,q2,k;
                        temp = global_to_local(temp);
                        curl2->Vz->F(temp) = curl2->Vz->F(temp) + exp(-pow((((i+0.5)*dx*kx + q2*dy*ky + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_x_dir/dy;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = i,q2,k;
                        temp = global_to_local(temp);
                        curl1->Vx->F(temp) = curl1->Vx->F(temp) + exp(-pow(((i*dx*kx + q2*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_z_dir/dy;
                    }
                }
            }



            // Changes on XY plane
            if ((r1-1)>=gridData.colloq_start_index_z && (r1-1)<=gridData.colloq_end_index_z && start_x!=-1 && start_y!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int j=start_y; j<=end_y; j++){
                        temp = i,j,(r1-1);
                        temp = global_to_local(temp);
                        curl1->Vy->F(temp) = curl1->Vy->F(temp) - exp(-pow((((i+0.5)*dx*kx + j*dy*ky + r1*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_x_dir/dz;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int j=start_y; j<=end_y-1; j++){
                        temp = i,j,(r1-1);
                        temp = global_to_local(temp);
                        curl2->Vx->F(temp) = curl2->Vx->F(temp) - exp(-pow(((i*dx*kx + (j+0.5)*dy*ky + r1*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dz;
                    }
                }
            }
            if (r2>=gridData.colloq_start_index_z && r2<=gridData.colloq_end_index_z && start_x!=-1 && start_y!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int j=start_y; j<=end_y; j++){
                        temp = i,j,r2;
                        temp = global_to_local(temp);
                        curl1->Vy->F(temp) = curl1->Vy->F(temp) + exp(-pow((((i+0.5)*dx*kx + j*dy*ky + r2*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_x_dir/dz;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int j=start_y; j<=end_y-1; j++){
                        temp = i,j,r2;
                        temp = global_to_local(temp);
                        curl2->Vx->F(temp) = curl2->Vx->F(temp) + exp(-pow(((i*dx*kx + (j+0.5)*dy*ky + r2*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dz;
                    }
                }
            }
        }
        else{
            // *************  if curl is not face-centered  *****************

            // Changes on YZ plane
            if (p1>=gridData.colloq_start_index_x && p1<=gridData.colloq_end_index_x && start_y!=-1 && start_z!=-1){
                for (int j=start_y; j<=end_y-1; j++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = p1,j,k;
                        temp = global_to_local(temp);
                        curl2->Vy->F(temp) = curl2->Vy->F(temp) - exp(-pow((((p1-0.5)*dx*kx + (j+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dx/c/mew0;
                    }
                }
                for (int j=start_y; j<=end_y; j++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = p1,j,k;
                        temp = global_to_local(temp);
                        curl1->Vz->F(temp) = curl1->Vz->F(temp) - exp(-pow((((p1-0.5)*dx*kx + j*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_y_dir/dx/c/mew0;
                    }
                }
            }
            if (p2>=gridData.colloq_start_index_x && p2<=gridData.colloq_end_index_x && start_y!=-1 && start_z!=-1){
                for (int j=start_y; j<=end_y-1; j++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = p2,j,k;
                        temp = global_to_local(temp);
                        curl2->Vy->F(temp) = curl2->Vy->F(temp) + exp(-pow((((p2+0.5)*dx*kx + (j+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dx/c/mew0;
                    }
                }
                for (int j=start_y; j<=end_y; j++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = p2,j,k;
                        temp = global_to_local(temp);
                        curl1->Vz->F(temp) = curl1->Vz->F(temp) + exp(-pow((((p2+0.5)*dx*kx + j*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_y_dir/dx/c/mew0;
                    }
                }
            }


            // Changes on XZ plane
            if (q1>=gridData.colloq_start_index_y && q1<=gridData.colloq_end_index_y && start_x!=-1 && start_z!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = i,q1,k;
                        temp = global_to_local(temp);
                        curl1->Vx->F(temp) = curl1->Vx->F(temp) - exp(-pow((((i+0.5)*dx*kx + (q1-0.5)*dy*ky + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dy/c/mew0;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = i,q1,k;
                        temp = global_to_local(temp);
                        curl2->Vz->F(temp) = curl2->Vz->F(temp) - exp(-pow(((i*dx*kx + (q1-0.5)*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dy/c/mew0;
                    }
                }
            }
            if (q2>=gridData.colloq_start_index_y && q2<=gridData.colloq_end_index_y && start_x!=-1 && start_z!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = i,q2,k;
                        temp = global_to_local(temp);
                        curl1->Vx->F(temp) = curl1->Vx->F(temp) + exp(-pow((((i+0.5)*dx*kx + (q2+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dy/c/mew0;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = i,q2,k;
                        temp = global_to_local(temp);
                        curl2->Vz->F(temp) = curl2->Vz->F(temp) + exp(-pow(((i*dx*kx + (q2+0.5)*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dy/c/mew0;
                    }
                }
            }



            // Changes on XY plane
            if (r1>=gridData.colloq_start_index_z && r1<=gridData.colloq_end_index_z && start_x!=-1 && start_y!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int j=start_y; j<=end_y; j++){
                        temp = i,j,r1;
                        temp = global_to_local(temp);
                        curl2->Vx->F(temp) = curl2->Vx->F(temp) - exp(-pow((((i+0.5)*dx*kx + j*dy*ky + (r1-0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_y_dir/dz/c/mew0;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int j=start_y; j<=end_y-1; j++){
                        temp = i,j,r1;
                        temp = global_to_local(temp);
                        curl1->Vy->F(temp) = curl1->Vy->F(temp) - exp(-pow(((i*dx*kx + (j+0.5)*dy*ky + (r1-0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dz/c/mew0;
                    }
                }
            }
            if (r2>=gridData.colloq_start_index_z && r2<=gridData.colloq_end_index_z && start_x!=-1 && start_y!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int j=start_y; j<=end_y; j++){
                        temp = i,j,r2;
                        temp = global_to_local(temp);
                        curl2->Vx->F(temp) = curl2->Vx->F(temp) + exp(-pow((((i+0.5)*dx*kx + j*dy*ky + (r2+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_y_dir/dz/c/mew0;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int j=start_y; j<=end_y-1; j++){
                        temp = i,j,r2;
                        temp = global_to_local(temp);
                        curl1->Vy->F(temp) = curl1->Vy->F(temp) + exp(-pow(((i*dx*kx + (j+0.5)*dy*ky + (r2+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dz/c/mew0;
                    }
                }
            }
        }
    }
}



void maxwell::plane_wave_initialise(){

    start_x=-1, end_x=-1, start_y=-1, end_y=-1, start_z=-1, end_z=-1;

    if (p1>gridData.colloq_end_index_x){
        start_x = -1;
        end_x = -1;
    }
    else if(p1>=gridData.colloq_start_index_x){
        start_x = p1;
        if (p2>gridData.colloq_end_index_x){
            end_x = gridData.colloq_end_index_x;
        }
        else{
            end_x = p2;
        }
    }
    else{
        if(p2>gridData.colloq_end_index_x){
            start_x = gridData.colloq_start_index_x;
            end_x = gridData.colloq_end_index_x;
        }
        else if (p2>=gridData.colloq_start_index_x){
            start_x = gridData.colloq_start_index_x;
            end_x = p2;
        }
        else{
            start_x = -1;
            end_x = -1;
        }
    }



    if (q1>gridData.colloq_end_index_y){
        start_y = -1;
        end_y = -1;
    }
    else if(q1>=gridData.colloq_start_index_y){
        start_y = q1;
        if (q2>gridData.colloq_end_index_y){
            end_y = gridData.colloq_end_index_y;
        }
        else{
            end_y = q2;
        }
    }
    else{
        if(q2>gridData.colloq_end_index_y){
            start_y = gridData.colloq_start_index_y;
            end_y = gridData.colloq_end_index_y;
        }
        else if (q2>=gridData.colloq_start_index_y){
            start_y = gridData.colloq_start_index_y;
            end_y = q2;
        }
        else{
            start_y = -1;
            end_y = -1;
        }
    }

    if (r1>gridData.colloq_end_index_z){
        start_z = -1;
        end_z = -1;
    }
    else if(r1>=gridData.colloq_start_index_z){
        start_z = r1;
        if (r2>gridData.colloq_end_index_z){
            end_z = gridData.colloq_end_index_z;
        }
        else{
            end_z = r2;
        }
    }
    else{
        if(r2>gridData.colloq_end_index_z){
            start_z = gridData.colloq_start_index_z;
            end_z = gridData.colloq_end_index_z;
        }
        else if (r2>=gridData.colloq_start_index_z){
            start_z = gridData.colloq_start_index_z;
            end_z = r2;
        }
        else{
            start_z = -1;
            end_z = -1;
        }
    }


    // std::cout<<"Rank: "<<mpi.rank<<", X-start: "<<start_x<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", X-end: "<<end_x<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Y-start: "<<start_y<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Y-end: "<<end_y<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Z-start: "<<start_z<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Z-end: "<<end_z<<std::endl;
    //
    // std::cout<<"Rank: "<<mpi.rank<<", Colloq X-start: "<<gridData.colloq_start_index_x<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Colloq X-end: "<<gridData.colloq_end_index_x<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Colloq Y-start: "<<gridData.colloq_start_index_y<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Colloq Y-end: "<<gridData.colloq_end_index_y<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Colloq Z-start: "<<gridData.colloq_start_index_z<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Colloq Z-end: "<<gridData.colloq_end_index_z<<std::endl;
}



bool maxwell::check_global_limits(blitz::TinyVector<int,3> v, bool xstag, bool ystag){
    int Nx_e = gridData.colloq_end_index_x;
    int Ny_e = gridData.colloq_end_index_y;
    if (xstag){
        Nx_e = Nx_e - 1;
    }
    if (ystag && (!gridData.isPlanar)){
        Ny_e = Ny_e - 1;
    }
    if ((v(0)-gridData.colloq_start_index_x)>=0 && (v(1)-gridData.colloq_start_index_y)>=0){
        if ((v(0)-Nx_e)<=0 && (v(1)-Ny_e)<=0){
            return true;
        }
    }
    return false;
}



blitz::TinyVector<int,3> maxwell::global_to_local(blitz::TinyVector<int,3> glob){
    blitz::TinyVector<int,3> loc = glob;
    loc(0) = glob(0) - gridData.colloq_start_index_x;
    loc(1) = glob(1) - gridData.colloq_start_index_y;
    return loc;
}



int maxwell::sigma_fn(int dim, bool is_E, blitz::Array<double,1> w, double d_pml, double m_pml,double sigma_max){
  double delta=0;
  int N=0;
  int N_local=0;
  int offset=0;
  if (dim==0){
    delta=gridData.inputData.dx;
    N_local=gridData.local_colloq_x;
    N=gridData.inputData.Nx;
    offset=gridData.colloq_start_index_x;
  }
  else if (dim==1){
    delta=gridData.inputData.dy;
    N_local=gridData.local_colloq_y;
    N=gridData.inputData.Ny;
    offset=gridData.colloq_start_index_y;
  }
  else if (dim==2){
    delta=gridData.inputData.dz;
    N_local=gridData.local_colloq_z;
    N=gridData.inputData.Nz;
    offset=0;
  }
  else{
    std::cout<<"Warning: Dim should be 1, 2 or 3. Returning Default!"<<std::endl;
    w=0;
    return 0;
  }

  blitz::firstIndex ii;
  blitz::Array<double,1> temp;
  if (is_E){
    temp.resize(N_local);
    temp=ii+offset;
  }
  else{
    temp.resize(N_local-1);
    temp=ii+0.5+offset;
  }

  double d=d_pml/delta;
  w=0.0;
  w=(temp<=(d-1))*pow((((d-1)-temp)/(d-1)),m_pml)*sigma_max+(temp>(d-1))*w;
  w=(temp>=(N-d))*pow((temp-N+d)/(d-1),m_pml)*sigma_max+(temp<(N-d))*w;
  return 0;
}



int maxwell::k_fn(int dim, bool is_E, blitz::Array<double,1> w, double d_pml, double m_pml,double k_max){
    double delta=0;
    int N=0;
    int N_local=0;
    int offset=0;
    if (dim==0){
      delta=gridData.inputData.dx;
      N_local=gridData.local_colloq_x;
      N=gridData.inputData.Nx;
      offset=gridData.colloq_start_index_x;
    }
    else if (dim==1){
      delta=gridData.inputData.dy;
      N_local=gridData.local_colloq_y;
      N=gridData.inputData.Ny;
      offset=gridData.colloq_start_index_y;
    }
    else if (dim==2){
      delta=gridData.inputData.dz;
      N_local=gridData.local_colloq_z;
      N=gridData.inputData.Nz;
      offset=0;
    }
    else{
      std::cout<<"Warning: Dim should be 1, 2 or 3. Returning Default!"<<std::endl;
      w=1.0;
      return 0;
    }

    blitz::firstIndex ii;
    blitz::Array<double,1> temp;
    if (is_E){
      temp.resize(N_local);
      temp=ii+offset;
    }
    else{
      temp.resize(N_local-1);
      temp=ii+0.5+offset;
    }


    double d=d_pml/delta;
    w=1;
    w=(temp<=(d-1))*(pow((((d-1)-temp)/(d-1)),m_pml)*(k_max-1)+1)+(temp>(d-1))*w;
    w=(temp>=(N-d))*(pow((temp-N+d)/(d-1),m_pml)*(k_max-1)+1)+(temp<(N-d))*w;
    return 0;

}



int maxwell::a_fn(int dim, bool is_E, blitz::Array<double,1> w, double d_pml, double ma_pml,double a_max){
    double delta=0;
    int N=0;
    int N_local=0;
    int offset=0;
    if (dim==0){
      delta=gridData.inputData.dx;
      N_local=gridData.local_colloq_x;
      N=gridData.inputData.Nx;
      offset=gridData.colloq_start_index_x;
    }
    else if (dim==1){
      delta=gridData.inputData.dy;
      N_local=gridData.local_colloq_y;
      N=gridData.inputData.Ny;
      offset=gridData.colloq_start_index_y;
    }
    else if (dim==2){
      delta=gridData.inputData.dz;
      N_local=gridData.local_colloq_z;
      N=gridData.inputData.Nz;
      offset=0;
    }
    else{
      std::cout<<"Warning: Dim should be 1, 2 or 3. Returning Default!"<<std::endl;
      w=0;
      return 0;
    }

    blitz::firstIndex ii;
    blitz::Array<double,1> temp;
    if (is_E){
      temp.resize(N_local);
      temp=ii+offset;
    }
    else{
      temp.resize(N_local-1);
      temp=ii+0.5+offset;
    }


    double d=d_pml/delta;
    w=0.0;
    w=(temp<=(d-1))*pow((temp/(d-1)),ma_pml)*a_max+(temp>(d-1))*w;
    w=(temp>=(N-d))*pow((N-1-temp)/(d-1),ma_pml)*a_max+(temp<(N-d))*w;
    return 0;
}



int maxwell::b_fn(int dim, bool is_E, blitz::Array<double,1> w, double d_pml, double ma_pml, double m_pml, double k_max, double sigma_max, double a_max){
  blitz::Array<double,1> temp_sigma, temp_a, temp_k;
  int N=0;
  if (dim==0){
    N=gridData.local_colloq_x;
  }
  else if (dim==1){
    N=gridData.local_colloq_y;
  }
  else if (dim==2){
    N=gridData.local_colloq_z;
  }
  else{
    std::cout<<"Warning: Dim should be 1, 2 or 3. Returning Default!"<<std::endl;
  }
  if (is_E){
    temp_sigma.resize(N);
    sigma_fn(dim,true,temp_sigma,d_pml,m_pml,sigma_max);
    temp_a.resize(N);
    a_fn(dim,true,temp_a,d_pml,ma_pml,a_max);
    temp_k.resize(N);
    k_fn(dim,true,temp_k,d_pml,m_pml,k_max);
  }
  else{
    temp_sigma.resize(N-1);
    sigma_fn(dim,false,temp_sigma,d_pml,m_pml,sigma_max);
    temp_a.resize(N-1);
    a_fn(dim,false,temp_a,d_pml,ma_pml,a_max);
    temp_k.resize(N-1);
    k_fn(dim,false,temp_k,d_pml,m_pml,k_max);
  }

  w=exp(-((temp_sigma/temp_k/epsilon0)-(temp_a/epsilon0))*dt);

  return 0;
}



int maxwell::c_fn(int dim, bool is_E, blitz::Array<double,1> w, double d_pml, double ma_pml, double m_pml, double k_max, double sigma_max, double a_max){
  blitz::Array<double,1> temp_sigma, temp_a, temp_k;
  int N=0;
  if (dim==0){
    N=gridData.local_colloq_x;
  }
  else if (dim==1){
    N=gridData.local_colloq_y;
  }
  else if (dim==2){
    N=gridData.local_colloq_z;
  }
  else{
    std::cout<<"Warning: Dim should be 1, 2 or 3. Returning Default!"<<std::endl;
  }
  if (is_E){
    temp_sigma.resize(N);
    sigma_fn(dim,true,temp_sigma,d_pml,m_pml,sigma_max);
    temp_a.resize(N);
    a_fn(dim,true,temp_a,d_pml,ma_pml,a_max);
    temp_k.resize(N);
    k_fn(dim,true,temp_k,d_pml,m_pml,k_max);
  }
  else{
    temp_sigma.resize(N-1);
    sigma_fn(dim,false,temp_sigma,d_pml,m_pml,sigma_max);
    temp_a.resize(N-1);
    a_fn(dim,false,temp_a,d_pml,ma_pml,a_max);
    temp_k.resize(N-1);
    k_fn(dim,false,temp_k,d_pml,m_pml,k_max);
  }

  b_fn(dim,is_E,w,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max);
  if (!is_E){
    N=N-1;
  }
  int i;
  for (i=0;i<N;i++){
    if (temp_sigma(i)==0){
      w(i)=0.0;
    }
    else{
      w(i)=temp_sigma(i)/(temp_sigma(i)*temp_k(i)+pow(temp_k(i),2)*temp_a(i))*(w(i)-1);
    }
  }
  return 0;
}
