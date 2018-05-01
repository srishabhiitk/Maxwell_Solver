#ifndef MAXWELL_H
#define MAXWELL_H

#include <iostream>
#include <math.h>
#include <blitz/array.h>
#include "parallel.h"
#include "reader.h"
#include "grid.h"
#include "vfield.h"
#include "writer.h"
#include <sstream>
#include <sys/time.h>

class maxwell {
    public:
        int n_threads;
        int num_timesteps;
        double epsilon0;
        double mew0;
        double S;
        double c;
        double dt, dx, dy, dz;
        double sigma, init;
        double lambda;
        reader &inputData;
        parallel &mpi;
        grid &gridData;
        blitz::Array<blitz::TinyVector<int,3>,1> H_x_pw_index, H_z_pw_index, E_y_pw_index1, E_y_pw_index2;
        blitz::Array<double,1> H_x_pw_coeff, H_z_pw_coeff, E_y_pw_coeff1, E_y_pw_coeff2;
        blitz::Array<int,1> H_x_pw_sign, H_z_pw_sign, E_y_pw_sign1, E_y_pw_sign2;
        double H_x_dir, H_z_dir, E_y_dir;
        double kx, ky, kz;
        blitz::TinyVector<int,3> pw_start_cord, pw_end_cord;
        int p1,p2,r1,r2;

        maxwell(reader &_inputData, parallel &_mpi, grid &_gridData);
        void solve();
        void solve3d();
        void solve3d_pml();
        void solve_planar();
        void solve_planar_pml();
        void plane_wave_execute(vfield *curl, int timestep);
        void plane_wave_initialise();
        int sigma_fn(int dim, bool is_E, blitz::Array<double,1> w, double d_pml, double m_pml,double sigma_max);
        int k_fn(int dim, bool is_E, blitz::Array<double,1> w, double d_pml, double m_pml,double k_max);
        int a_fn(int dim, bool is_E, blitz::Array<double,1> w, double d_pml, double ma_pml,double a_max);
        int b_fn(int dim, bool is_E, blitz::Array<double,1> w, double d_pml, double ma_pml, double m_pml, double k_max, double sigma_max, double a_max);
        int c_fn(int dim, bool is_E, blitz::Array<double,1> w, double d_pml, double ma_pml, double m_pml, double k_max, double sigma_max, double a_max);
        bool check_global_limits(blitz::TinyVector<int,3> v, bool xstag, bool ystag);
        blitz::TinyVector<int,3> global_to_local(blitz::TinyVector<int,3> glob);
};



#endif
