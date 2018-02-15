#include <iostream>
#include <blitz/array.h>
#include <math.h>
#include <sys/time.h>
#include <sstream>
// #include "H5Cpp.h"
// using namespace H5;

int num_processors=8;


class field{
public:
  int size;
  bool is2D;
  blitz::Array<double,3> x,y,z;
  bool staggered;

  field(int N, bool staggered_, bool is2D_){
    size=N;
    is2D=is2D_;
    staggered=staggered_;
    if (is2D){
      if (staggered_){
        x.resize(N,N-1,1);
        y.resize(N-1,N,1);
        z.resize(0,0,0);
      }
      else{
        x.resize(0,0,0);
        y.resize(0,0,0);
        z.resize(N,N,1);
      }
    }
    else{
      if (staggered_){
        x.resize(N-1,N,N);
        y.resize(N,N-1,N);
        z.resize(N,N,N-1);
      }
      else{
        x.resize(N,N-1,N-1);
        y.resize(N-1,N,N-1);
        z.resize(N-1,N-1,N);
      }
    }
  }

  void intialize_to_zero(){
    x=0;
    y=0;
    z=0;
  }
};



blitz::RectDomain<3> shift(blitz::RectDomain<3> rd, int dim, int steps){
  rd.lbound()(dim)+=steps;
  rd.ubound()(dim)+=steps;
  return rd;
}





void curl_3d(field *F, field *Fcurl, double dx, double dy, double dz){
  //Checking the requirements
  if (!(F->is2D) && (F->staggered^Fcurl->staggered)){
    int N=F->size;
    blitz::RectDomain<3> rd;
    blitz::TinyVector<int,3> loBound,upBound;
    Fcurl->intialize_to_zero();
    if (F->staggered){
      upBound=N-2,N-2,N-2;

      loBound=1,0,0;
      rd=blitz::RectDomain<3>(loBound,upBound);
      Fcurl->x(rd)=(F->z(shift(rd,1,1))-F->z(rd))/dy-(F->y(shift(rd,2,1))-F->y(rd))/dz;

      loBound=0,1,0;
      rd=blitz::RectDomain<3>(loBound,upBound);
      Fcurl->y(rd)=(F->x(shift(rd,2,1))-F->x(rd))/dz-(F->z(shift(rd,0,1))-F->z(rd))/dx;

      loBound=0,0,1;
      rd=blitz::RectDomain<3>(loBound,upBound);
      Fcurl->z(rd)=(F->y(shift(rd,0,1))-F->y(rd))/dx-(F->x(shift(rd,1,1))-F->x(rd))/dy;
    }
    else{
      upBound=N-2,N-2,N-2;

      loBound=0,1,1;
      rd=blitz::RectDomain<3>(loBound,upBound);
      Fcurl->x(rd)=(F->z(rd)-F->z(shift(rd,1,-1)))/dy-(F->y(rd)-F->y(shift(rd,2,-1)))/dz;

      loBound=1,0,1;
      rd=blitz::RectDomain<3>(loBound,upBound);
      Fcurl->y(rd)=(F->x(rd)-F->x(shift(rd,2,-1)))/dz-(F->z(rd)-F->z(shift(rd,0,-1)))/dx;

      loBound=1,1,0;
      rd=blitz::RectDomain<3>(loBound,upBound);
      Fcurl->z(rd)=(F->y(rd)-F->y(shift(rd,0,-1)))/dx-(F->x(rd)-F->x(shift(rd,1,-1)))/dy;
    }
  }
  //Error if requirements of curl_3d are not met!
  else{
    if ((F->is2D)){
      std::cout<<"Error! curl_3d is defined for 3D fields only."<<std::endl;
    }
    else{
      std::cout<<"First two parameter should be staggered and unstaggered."<<std::endl;
    }

  }
}


void curl_3d_adv(field *F, field *Fcurl1, field *Fcurl2, double dx, double dy, double dz){
  //Checking the requirements
  if (!(F->is2D) && (F->staggered^Fcurl1->staggered)&& (F->staggered^Fcurl2->staggered)){
    int N=F->size;
    // Fcurl1->intialize_to_zero();
    // Fcurl2->intialize_to_zero();
    int i,j,k;
    if (F->staggered){
      #pragma omp parallel for num_threads(num_processors) schedule(dynamic)
      for (int i=0;i<N;i++){
        for (int j=0;j<(N-1);j++){
          for (int k=0;k<(N-1);k++){
            Fcurl1->x(i,j,k)=(F->z(i,j+1,k)-F->z(i,j,k))/dy;
            Fcurl2->x(i,j,k)=(F->y(i,j,k+1)-F->y(i,j,k))/dz;
          }
        }
      }
      #pragma omp parallel for num_threads(num_processors) schedule(dynamic)
      for (int i=0;i<(N-1);i++){
        for (int j=0;j<N;j++){
          for (int k=0;k<(N-1);k++){
            Fcurl1->y(i,j,k)=(F->x(i,j,k+1)-F->x(i,j,k))/dz;
            Fcurl2->y(i,j,k)=(F->z(i+1,j,k)-F->z(i,j,k))/dx;
          }
        }
      }
      #pragma omp parallel for num_threads(num_processors) schedule(dynamic)
      for (int i=0;i<(N-1);i++){
        for (int j=0;j<(N-1);j++){
          for (int k=0;k<N;k++){
            Fcurl1->z(i,j,k)=(F->y(i+1,j,k)-F->y(i,j,k))/dx;
            Fcurl2->z(i,j,k)=(F->x(i,j+1,k)-F->x(i,j,k))/dy;
          }
        }
      }
    }
    else{
      #pragma omp parallel for num_threads(num_processors) schedule(dynamic)
      for (int i=0;i<(N-1);i++){
        for (int j=1;j<(N-1);j++){
          for (int k=1;k<(N-1);k++){
            Fcurl1->x(i,j,k)=(F->z(i,j,k)-F->z(i,j-1,k))/dy;
            Fcurl2->x(i,j,k)=(F->y(i,j,k)-F->y(i,j,k-1))/dz;
          }
        }
      }
      #pragma omp parallel for num_threads(num_processors) schedule(dynamic)
      for (int i=1;i<(N-1);i++){
        for (int j=0;j<(N-1);j++){
          for (int k=1;k<(N-1);k++){
            Fcurl1->y(i,j,k)=(F->x(i,j,k)-F->x(i,j,k-1))/dz;
            Fcurl2->y(i,j,k)=(F->z(i,j,k)-F->z(i-1,j,k))/dx;
          }
        }
      }
      #pragma omp parallel for num_threads(num_processors) schedule(dynamic)
      for (int i=1;i<(N-1);i++){
        for (int j=1;j<(N-1);j++){
          for (int k=0;k<(N-1);k++){
            Fcurl1->z(i,j,k)=(F->y(i,j,k)-F->y(i-1,j,k))/dx;
            Fcurl2->z(i,j,k)=(F->x(i,j,k)-F->x(i,j-1,k))/dy;
          }
        }
      }
    }
  }
  //Error if requirements of curl_3d are not met!
  else{
    if ((F->is2D)){
      std::cout<<"Error! curl_3d is defined for 3D fields only."<<std::endl;
    }
    else{
      std::cout<<"First two parameter should be staggered and unstaggered."<<std::endl;
    }

  }
}




using blitz::Range;


void plane_wave_as_BC(field *E,int n,int timestep,double dx,double dy,double dz,double dt,double c,double sigma,double init){
  int i,j,k;
  for (j=0;j<(n-1);j++){
    for (k=0;k<(n-1);k++){
      i=0;
      E->x(i,j,k)=exp(-pow(((i*dx+(j+0.5)*dy+(k+0.5)*dz)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(6,0.5);
      E->y(j,i,k)=exp(-pow(((i*dy+(j+0.5)*dx+(k+0.5)*dz)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(1.5,0.5);
      E->z(k,j,i)=exp(-pow(((i*dz+(j+0.5)*dy+(k+0.5)*dx)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(6,0.5);
      i=n-1;
      E->x(i,j,k)=exp(-pow(((i*dx+(j+0.5)*dy+(k+0.5)*dz)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(6,0.5);
      E->y(j,i,k)=exp(-pow(((i*dy+(j+0.5)*dx+(k+0.5)*dz)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(1.5,0.5);
      E->z(k,j,i)=exp(-pow(((i*dz+(j+0.5)*dy+(k+0.5)*dx)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(6,0.5);
    }
  }
  for (k=0;k<(n-1);k++){
    for (i=1;i<(n-1);i++){
      j=0;
      E->x(i,j,k)=exp(-pow(((i*dx+(j+0.5)*dy+(k+0.5)*dz)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(6,0.5);
      E->y(j,i,k)=-exp(-pow(((i*dy+(j+0.5)*dx+(k+0.5)*dz)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(1.5,0.5);
      E->z(k,j,i)=exp(-pow(((i*dz+(j+0.5)*dy+(k+0.5)*dx)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(6,0.5);
      j=n-2;
      E->x(i,j,k)=exp(-pow(((i*dx+(j+0.5)*dy+(k+0.5)*dz)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(6,0.5);
      E->y(j,i,k)=-exp(-pow(((i*dy+(j+0.5)*dx+(k+0.5)*dz)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(1.5,0.5);
      E->z(k,j,i)=exp(-pow(((i*dz+(j+0.5)*dy+(k+0.5)*dx)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(6,0.5);
    }
  }
  for (j=1;j<(n-2);j++){
    for (i=1;i<(n-1);i++){
      k=0;
      E->x(i,j,k)=exp(-pow(((i*dx+(j+0.5)*dy+(k+0.5)*dz)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(6,0.5);
      E->y(j,i,k)=-exp(-pow(((i*dy+(j+0.5)*dx+(k+0.5)*dz)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(1.5,0.5);
      E->z(k,j,i)=exp(-pow(((i*dz+(j+0.5)*dy+(k+0.5)*dx)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(6,0.5);
      k=n-2;
      E->x(i,j,k)=exp(-pow(((i*dx+(j+0.5)*dy+(k+0.5)*dz)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(6,0.5);
      E->y(j,i,k)=-exp(-pow(((i*dy+(j+0.5)*dx+(k+0.5)*dz)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(1.5,0.5);
      E->z(k,j,i)=exp(-pow(((i*dz+(j+0.5)*dy+(k+0.5)*dx)/pow(3,0.5)+init*dx-timestep*c*dt),2)/sigma)/pow(6,0.5);
    }
  }
}


int sigma_fn(int dim, bool is_E, blitz::Array<double,1> w, double d_pml, double m_pml,double sigma_max,double dx, double dy, double dz, int N){
  double delta=0;
  if (dim==0){
    delta=dx;
  }
  else if (dim==1){
    delta=dy;
  }
  else if (dim==2){
    delta=dz;
  }
  else{
    std::cout<<"Dim should be 1, 2 or 3. Returning Default!"<<std::endl;
    w=0;
    return 0;
  }

  blitz::firstIndex ii;
  blitz::Array<double,1> temp;
  if (is_E){
    temp.resize(N);
    temp=ii;
  }
  else{
    temp.resize(N-1);
    temp=ii+0.5;
  }

  double d=d_pml/delta;
  w=0.0;
  w=(temp<=(d-1))*pow((((d-1)-temp)/(d-1)),m_pml)*sigma_max+(temp>(d-1))*w;
  w=(temp>=(N-d))*pow((temp-N+d)/(d-1),m_pml)*sigma_max+(temp<(N-d))*w;
  return 0;
}


int k_fn(int dim, bool is_E, blitz::Array<double,1> w, double d_pml, double m_pml,double k_max,double dx, double dy, double dz, int N){
  double delta=0;
  if (dim==0){
    delta=dx;
  }
  else if (dim==1){
    delta=dy;
  }
  else if (dim==2){
    delta=dz;
  }
  else{
    std::cout<<"Dim should be 1, 2 or 3. Returning Default!"<<std::endl;
    w=1.0;
    return 0;
  }

  blitz::firstIndex ii;
  blitz::Array<double,1> temp;
  if (is_E){
    temp.resize(N);
    temp=ii;
  }
  else{
    temp.resize(N-1);
    temp=ii+0.5;
  }

  double d=d_pml/delta;
  w=1;
  w=(temp<=(d-1))*(pow((((d-1)-temp)/(d-1)),m_pml)*(k_max-1)+1)+(temp>(d-1))*w;
  w=(temp>=(N-d))*(pow((temp-N+d)/(d-1),m_pml)*(k_max-1)+1)+(temp<(N-d))*w;
  return 0;

}


int a_fn(int dim, bool is_E, blitz::Array<double,1> w, double d_pml, double ma_pml,double a_max,double dx, double dy, double dz, int N){
  double delta=0;
  if (dim==0){
    delta=dx;
  }
  else if (dim==1){
    delta=dy;
  }
  else if (dim==2){
    delta=dz;
  }
  else{
    std::cout<<"Dim should be 1, 2 or 3. Returning Default!"<<std::endl;
    w=0.0;
    return 0;
  }

  blitz::firstIndex ii;
  blitz::Array<double,1> temp;
  if (is_E){
    temp.resize(N);
    temp=ii;
  }
  else{
    temp.resize(N-1);
    temp=ii+0.5;
  }

  double d=d_pml/delta;
  w=0.0;
  w=(temp<=(d-1))*pow((temp/(d-1)),ma_pml)*a_max+(temp>(d-1))*w;
  w=(temp>=(N-d))*pow((N-1-temp)/(d-1),ma_pml)*a_max+(temp<(N-d))*w;
  return 0;
}


int b_fn(int dim, bool is_E, blitz::Array<double,1> w, double d_pml, double ma_pml, double m_pml, double k_max, double sigma_max, double a_max, double epsilon0, double dt, double dx, double dy, double dz, int N){
  double delta=0;
  if (dim==0){
    delta=dx;
  }
  else if (dim==1){
    delta=dy;
  }
  else if (dim==2){
    delta=dz;
  }
  else{
    std::cout<<"Dim should be 1, 2 or 3. Returning Default!"<<std::endl;
    w=0.0;
    return 0;
  }

  blitz::firstIndex ii;
  blitz::Array<double,1> temp_sigma, temp_a, temp_k;
  if (is_E){
    temp_sigma.resize(N);
    sigma_fn(dim,true,temp_sigma,d_pml,m_pml,sigma_max,dx,dy,dz,N);
    temp_a.resize(N);
    a_fn(dim,true,temp_a,d_pml,ma_pml,a_max,dx,dy,dz,N);
    temp_k.resize(N);
    k_fn(dim,true,temp_k,d_pml,m_pml,k_max,dx,dy,dz,N);
  }
  else{
    temp_sigma.resize(N-1);
    sigma_fn(dim,false,temp_sigma,d_pml,m_pml,sigma_max,dx,dy,dz,N);
    temp_a.resize(N-1);
    a_fn(dim,false,temp_a,d_pml,ma_pml,a_max,dx,dy,dz,N);
    temp_k.resize(N-1);
    k_fn(dim,false,temp_k,d_pml,m_pml,k_max,dx,dy,dz,N);
  }



  w=exp(-((temp_sigma/temp_k/epsilon0)-(temp_a/epsilon0))*dt);
  return 0;
}


int c_fn(int dim, bool is_E, blitz::Array<double,1> w, double d_pml, double ma_pml, double m_pml, double k_max, double sigma_max, double a_max, double epsilon0, double dt, double dx, double dy, double dz, int N){
  double delta=0;
  if (dim==0){
    delta=dx;
  }
  else if (dim==1){
    delta=dy;
  }
  else if (dim==2){
    delta=dz;
  }
  else{
    std::cout<<"Dim should be 1, 2 or 3. Returning Default!"<<std::endl;
    w=0;
    return 0;
  }

  blitz::firstIndex ii;
  blitz::Array<double,1> temp_sigma, temp_a, temp_k;
  if (is_E){
    temp_sigma.resize(N);
    sigma_fn(dim,true,temp_sigma,d_pml,m_pml,sigma_max,dx,dy,dz,N);
    temp_a.resize(N);
    a_fn(dim,true,temp_a,d_pml,ma_pml,a_max,dx,dy,dz,N);
    temp_k.resize(N);
    k_fn(dim,true,temp_k,d_pml,m_pml,k_max,dx,dy,dz,N);
  }
  else{
    temp_sigma.resize(N-1);
    sigma_fn(dim,false,temp_sigma,d_pml,m_pml,sigma_max,dx,dy,dz,N);
    temp_a.resize(N-1);
    a_fn(dim,false,temp_a,d_pml,ma_pml,a_max,dx,dy,dz,N);
    temp_k.resize(N-1);
    k_fn(dim,false,temp_k,d_pml,m_pml,k_max,dx,dy,dz,N);
  }

  // std::cout << "k_ex:" << std::endl;
  // std::cout << temp_k << std::endl;
  // std::cout << "k_ey:" << std::endl;
  // std::cout << temp_k << std::endl;
  // std::cout << "k_ez:" << std::endl;
  // std::cout << temp_k << std::endl;
  // std::cout << "k_hx:" << std::endl;
  // std::cout << temp_k << std::endl;
  // std::cout << "k_hy:" << std::endl;
  // std::cout << temp_k << std::endl;
  // std::cout << "k_hz:" << std::endl;
  // std::cout << temp_k << std::endl;

  b_fn(dim,is_E,w,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max,epsilon0,dt,dx,dy,dz,N);
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




int main(){
  int num_timesteps=50;
  int N=52;
  double epsilon0=8.85*pow(10,-12);
  double mew0=M_PI*4*pow(10,-7);

  double c=1/pow(epsilon0*mew0,0.5); //Speed of light
  double S=1/pow(3,0.5); //Courant factor
  double dx=pow(10,-3);
  double dy=dx, dz=dx;
  double dt=S*dx/c;
  double Ca=1.0,Da=1.0;
  // double Cb=dt/epsilon0;
  // double Db=dt/mew0;
  double sigma=pow(10,0);
  double init=20.0;



  // For PML BC
  // double a_max=0.1;
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

  struct timeval begin,end;

  blitz::Array<double,1> k_ex(N),k_ey(N),k_ez(N);
  blitz::Array<double,1> k_hx(N-1),k_hy(N-1),k_hz(N-1);

  k_fn(0,true,k_ex,d_pml,m_pml,k_max,dx,dy,dz,N);
  k_fn(1,true,k_ey,d_pml,m_pml,k_max,dx,dy,dz,N);
  k_fn(2,true,k_ez,d_pml,m_pml,k_max,dx,dy,dz,N);
  k_fn(0,false,k_hx,d_pml,m_pml,k_max,dx,dy,dz,N);
  k_fn(1,false,k_hy,d_pml,m_pml,k_max,dx,dy,dz,N);
  k_fn(2,false,k_hz,d_pml,m_pml,k_max,dx,dy,dz,N);





  blitz::Array<double,1> b_ex(N),b_ey(N),b_ez(N),c_ex(N),c_ey(N),c_ez(N);
  blitz::Array<double,1> b_hx(N-1),b_hy(N-1),b_hz(N-1),c_hx(N-1),c_hy(N-1),c_hz(N-1);

  b_fn(0,true,b_ex,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max,epsilon0,dt,dx,dy,dz,N);
  b_fn(1,true,b_ey,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max,epsilon0,dt,dx,dy,dz,N);
  b_fn(2,true,b_ez,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max,epsilon0,dt,dx,dy,dz,N);
  c_fn(0,true,c_ex,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max,epsilon0,dt,dx,dy,dz,N);
  c_fn(1,true,c_ey,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max,epsilon0,dt,dx,dy,dz,N);
  c_fn(2,true,c_ez,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max,epsilon0,dt,dx,dy,dz,N);
  b_fn(0,false,b_hx,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max,epsilon0,dt,dx,dy,dz,N);
  b_fn(1,false,b_hy,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max,epsilon0,dt,dx,dy,dz,N);
  b_fn(2,false,b_hz,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max,epsilon0,dt,dx,dy,dz,N);
  c_fn(0,false,c_hx,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max,epsilon0,dt,dx,dy,dz,N);
  c_fn(1,false,c_hy,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max,epsilon0,dt,dx,dy,dz,N);
  c_fn(2,false,c_hz,d_pml,ma_pml,m_pml,k_max,sigma_max,a_max,epsilon0,dt,dx,dy,dz,N);



  field *E=new field(N,true,false);
  E->intialize_to_zero();
  field *Ecurl1=new field(N,false,false);
  Ecurl1->intialize_to_zero();
  field *Ecurl2=new field(N,false,false);
  Ecurl2->intialize_to_zero();
  field *Etotal=new field(N,true,false);
  Etotal->intialize_to_zero();
  field *Cb=new field(N,true,false);
  Cb->intialize_to_zero();


  field *H=new field(N,false,false);
  H->intialize_to_zero();
  field *Hcurl1=new field(N,true,false);
  Hcurl1->intialize_to_zero();
  field *Hcurl2=new field(N,true,false);
  Hcurl2->intialize_to_zero();
  field *Db=new field(N,false,false);
  Db->intialize_to_zero();


  field *psi_E1=new field(N,true,false);
  psi_E1->intialize_to_zero();
  field *psi_E2=new field(N,true,false);
  psi_E2->intialize_to_zero();
  field *psi_H1=new field(N,false,false);
  psi_H1->intialize_to_zero();
  field *psi_H2=new field(N,false,false);
  psi_H2->intialize_to_zero();



  Cb->x=dt/epsilon0;
  Cb->y=dt/epsilon0;
  Cb->z=dt/epsilon0;
  Db->x=dt/mew0;
  Db->y=dt/mew0;
  Db->z=dt/mew0;


  // blitz::RectDomain<3> rd;
  // blitz::TinyVector<int,3> loBound,upBound;

  //Solid block at center
  // loBound=45,45,45;
  // upBound=55,55,55;
  // rd=blitz::RectDomain<3>(loBound,upBound);
  // Cb->z(shift(rd,2,-1))=Cb->z(shift(rd,2,-1))/pow(2,2);
  // Cb->y(shift(rd,1,-1))=Cb->y(shift(rd,1,-1))/pow(2,2);
  // Cb->x(shift(rd,0,-1))=Cb->x(shift(rd,0,-1))/pow(2,2);


  // //HDF5 file declaration
  // H5File file("Ezdata2.h5", H5F_ACC_TRUNC);
	// hsize_t dims[3];               // dataset dimensions
  // std::stringstream set_name;



  gettimeofday(&begin,NULL);

  int t,i,j,k;
  for (t=0; t<num_timesteps; t++){
    if (t==0){
      // E->x(N/2,N/2,N/2)=1.0/pow(3,0.5);
      // E->y(N/2,N/2,N/2)=1.0/pow(3,0.5);
      // E->z(N/2,N/2,N/2)=1.0/pow(3,0.5);
      // E->z(N/2,N/2,N/2)=1.0;
      E->z(N/2,N/2,N/2)=exp(-pow(-init+t,2)/sigma);
      // plane_wave_as_BC(E,N,t,dx,dy,dz,dt,c,sigma,init);
    }


    // //Solid block at center
    // loBound=45,45,45;
    // upBound=55,55,55;
    // rd=blitz::RectDomain<3>(loBound,upBound);
    // E->z(shift(rd,2,-1))=0;
    // E->y(shift(rd,1,-1))=0;
    // E->x(shift(rd,0,-1))=0;
    // Eadv->z(shift(rd,2,-1))=0;
    // Eadv->y(shift(rd,1,-1))=0;
    // Eadv->x(shift(rd,0,-1))=0;



    curl_3d_adv(E,Ecurl1,Ecurl2,dx,dy,dz);

    #pragma omp parallel for num_threads(num_processors) schedule(dynamic)
    for (int i=0;i<N;i++){
      for (int j=0;j<(N-1);j++){
        for (int k=0;k<(N-1);k++){
          psi_H1->x(i,j,k)=b_hy(j)*psi_H1->x(i,j,k)+c_hy(j)*Ecurl1->x(i,j,k);
          psi_H2->x(i,j,k)=b_hz(k)*psi_H2->x(i,j,k)+c_hz(k)*Ecurl2->x(i,j,k);
          H->x(i,j,k)=Da*H->x(i,j,k)-Db->x(i,j,k)*(Ecurl1->x(i,j,k)/k_hy(j)-Ecurl2->x(i,j,k)/k_hz(k))-Db->x(i,j,k)*(psi_H1->x(i,j,k)-psi_H2->x(i,j,k));
        }
      }
    }

    #pragma omp parallel for num_threads(num_processors) schedule(dynamic)
    for (int i=0;i<(N-1);i++){
      for (int j=0;j<N;j++){
        for (int k=0;k<(N-1);k++){
          psi_H1->y(i,j,k)=b_hz(k)*psi_H1->y(i,j,k)+c_hz(k)*Ecurl1->y(i,j,k);
          psi_H2->y(i,j,k)=b_hx(i)*psi_H2->y(i,j,k)+c_hx(i)*Ecurl2->y(i,j,k);
          H->y(i,j,k)=Da*H->y(i,j,k)-Db->y(i,j,k)*(Ecurl1->y(i,j,k)/k_hz(k)-Ecurl2->y(i,j,k)/k_hx(i))-Db->y(i,j,k)*(psi_H1->y(i,j,k)-psi_H2->y(i,j,k));
        }
      }
    }

    #pragma omp parallel for num_threads(num_processors) schedule(dynamic)
    for (int i=0;i<(N-1);i++){
      for (int j=0;j<(N-1);j++){
        for (int k=0;k<N;k++){
          psi_H1->z(i,j,k)=b_hx(i)*psi_H1->z(i,j,k)+c_hx(i)*Ecurl1->z(i,j,k);
          psi_H2->z(i,j,k)=b_hy(j)*psi_H2->z(i,j,k)+c_hy(j)*Ecurl2->z(i,j,k);
          H->z(i,j,k)=Da*H->z(i,j,k)-Db->z(i,j,k)*(Ecurl1->z(i,j,k)/k_hx(i)-Ecurl2->z(i,j,k)/k_hy(j))-Db->z(i,j,k)*(psi_H1->z(i,j,k)-psi_H2->z(i,j,k));
        }
      }
    }




    curl_3d_adv(H,Hcurl1,Hcurl2,dx,dy,dz);

    #pragma omp parallel for num_threads(num_processors) schedule(dynamic)
    for (int i=0;i<(N-1);i++){
      for (int j=0;j<N;j++){
        for (int k=0;k<N;k++){
          psi_E1->x(i,j,k)=b_ey(j)*psi_E1->x(i,j,k)+c_ey(j)*Hcurl1->x(i,j,k);
          psi_E2->x(i,j,k)=b_ez(k)*psi_E2->x(i,j,k)+c_ez(k)*Hcurl2->x(i,j,k);
          E->x(i,j,k)=Ca*E->x(i,j,k)+Cb->x(i,j,k)*(Hcurl1->x(i,j,k)/k_ey(j)-Hcurl2->x(i,j,k)/k_ez(k))+Cb->x(i,j,k)*(psi_E1->x(i,j,k)-psi_E2->x(i,j,k));
        }
      }
    }

    #pragma omp parallel for num_threads(num_processors) schedule(dynamic)
    for (int i=0;i<N;i++){
      for (int j=0;j<(N-1);j++){
        for (int k=0;k<N;k++){
          psi_E1->y(i,j,k)=b_ez(k)*psi_E1->y(i,j,k)+c_ez(k)*Hcurl1->y(i,j,k);
          psi_E2->y(i,j,k)=b_ex(i)*psi_E2->y(i,j,k)+c_ex(i)*Hcurl2->y(i,j,k);
          E->y(i,j,k)=Ca*E->y(i,j,k)+Cb->y(i,j,k)*(Hcurl1->y(i,j,k)/k_ez(k)-Hcurl2->y(i,j,k)/k_ex(i))+Cb->y(i,j,k)*(psi_E1->y(i,j,k)-psi_E2->y(i,j,k));
        }
      }
    }

    #pragma omp parallel for num_threads(num_processors) schedule(dynamic)
    for (int i=0;i<N;i++){
      for (int j=0;j<N;j++){
        for (int k=0;k<(N-1);k++){
          psi_E1->z(i,j,k)=b_ex(i)*psi_E1->z(i,j,k)+c_ex(i)*Hcurl1->z(i,j,k);
          psi_E2->z(i,j,k)=b_ey(j)*psi_E2->z(i,j,k)+c_ey(j)*Hcurl2->z(i,j,k);
          E->z(i,j,k)=Ca*E->z(i,j,k)+Cb->z(i,j,k)*(Hcurl1->z(i,j,k)/k_ex(i)-Hcurl2->z(i,j,k)/k_ey(j))+Cb->z(i,j,k)*(psi_E1->z(i,j,k)-psi_E2->z(i,j,k));
        }
      }
    }




    #pragma omp parallel for num_threads(num_processors) schedule(dynamic)
    for (int i=0;i<(N-1);i++){
      for (int j=0;j<(N-1);j++){
        for (int k=0;k<(N-1);k++){
          Etotal->z(i,j,k)=E->z(i,j,k)*E->z(i,j,k)+E->y(i,j,k)*E->y(i,j,k)+E->x(i,j,k)*E->x(i,j,k);
        }
      }
    }





    // E->z(N/2,N/2,N/2)=1.0;
    // Eadv->z(N/2,N/2,N/2)=1.0;
    E->z(N/2,N/2,N/2)=exp(-pow(-init+t,2)/sigma);
    // Eadv->z(N/2,N/2,N/2)=exp(-pow(-init+t,2)/sigma);
    // plane_wave_as_BC(E,N,t,dx,dy,dz,dt,c,sigma,init);
    // plane_wave_as_BC(Eadv,N,t,dx,dy,dz,dt,c,sigma,init);




    // //Solid block at center
    // loBound=45,45,45;
    // upBound=55,55,55;
    // rd=blitz::RectDomain<3>(loBound,upBound);
    // E->z(shift(rd,2,-1))=0;
    // E->y(shift(rd,1,-1))=0;
    // E->x(shift(rd,0,-1))=0;
    // Eadv->z(shift(rd,2,-1))=0;
    // Eadv->y(shift(rd,1,-1))=0;
    // Eadv->x(shift(rd,0,-1))=0;




    // //Writing to file
    // if (t%3==0&&t>=10&&t<=300){
    //   dims[0] = E->z.shape()[0];
    // 	dims[1] = E->z.shape()[1];
    //   dims[2] = E->z.shape()[2];
    // 	DataSpace dataspace(3, dims);
    //   set_name<<"tz"<<t;
    // 	DataSet dataset = file.createDataSet(set_name.str(), PredType::NATIVE_DOUBLE, dataspace);
    //   dataset.write(Etotal->z.data(), PredType::NATIVE_DOUBLE);
    //   set_name.str("");
    // }


    // std::cout<<E->z(blitz::Range(45,55),blitz::Range(45,55),N/2)<<std::endl;

    std::cout<<t<<", max="<<blitz::max(Etotal->z)<<std::endl;
    // std::cout<<E->z(rd)<<std::endl;
  }
  std::cout<<"max="<<blitz::max(E->z)<<std::endl;

  gettimeofday(&end,NULL);
  std::cout<<"Time taken: "<<((end.tv_sec-begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.0e6<<" sec"<<std::endl;


  return 0;
}
