#include <iostream>
#include <blitz/array.h>
#include <math.h>
#include <sstream>
#include "H5Cpp.h"
using namespace H5;




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
      loBound=1,1,1;

      upBound=N-2,N-3,N-3;
      rd=blitz::RectDomain<3>(loBound,upBound);
      Fcurl->x(rd)=(F->z(shift(rd,1,1))-F->z(rd))/dy-(F->y(shift(rd,2,1))-F->y(rd))/dz;

      upBound=N-3,N-2,N-3;
      rd=blitz::RectDomain<3>(loBound,upBound);
      Fcurl->y(rd)=(F->x(shift(rd,2,1))-F->x(rd))/dz-(F->z(shift(rd,0,1))-F->z(rd))/dx;

      upBound=N-3,N-3,N-2;
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





int main(){
  int num_timesteps=600;
  int N=100;
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
  double sigma=pow(10,-4);
  double init=100.0;









  field *E=new field(N,false,false);
  E->intialize_to_zero();
  field *Eadv=new field(N,false,false);
  Eadv->intialize_to_zero();
  field *Ecurl=new field(N,true,false);
  Ecurl->intialize_to_zero();
  field *Etotal=new field(N,false,false);
  Etotal->intialize_to_zero();
  field *Cb=new field(N,false,false);
  Cb->intialize_to_zero();


  field *H=new field(N,true,false);
  H->intialize_to_zero();
  field *Hadv=new field(N,true,false);
  Hadv->intialize_to_zero();
  field *Hcurl=new field(N,false,false);
  Hcurl->intialize_to_zero();
  field *Db=new field(N,true,false);
  Db->intialize_to_zero();


  field *temp;//for swapping


  Cb->x=dt/epsilon0;
  Cb->y=dt/epsilon0;
  Cb->z=dt/epsilon0;
  Db->x=dt/mew0;
  Db->y=dt/mew0;
  Db->z=dt/mew0;


  blitz::RectDomain<3> rd;
  blitz::TinyVector<int,3> loBound,upBound;

  //Solid block at center
  loBound=45,45,45;
  upBound=55,55,55;
  rd=blitz::RectDomain<3>(loBound,upBound);
  Cb->z(shift(rd,2,-1))=Cb->z(shift(rd,2,-1))/pow(2,2);
  Cb->y(shift(rd,1,-1))=Cb->y(shift(rd,1,-1))/pow(2,2);
  Cb->x(shift(rd,0,-1))=Cb->x(shift(rd,0,-1))/pow(2,2);



  H5File file("Ezdata2.h5", H5F_ACC_TRUNC);
	hsize_t dims[3];               // dataset dimensions
  std::stringstream set_name;





  int i;
  for (i=0; i<num_timesteps; i++){
    if (i==0){
      // E->x(N/2,N/2,N/2)=1.0/pow(3,0.5);
      // E->y(N/2,N/2,N/2)=1.0/pow(3,0.5);
      // E->z(N/2,N/2,N/2)=1.0/pow(3,0.5);
      // E->z(N/2,N/2,N/2)=1.0;
      // E->z(N/2,N/2,N/2)=exp(-pow(-init+i,2)/sigma);
      plane_wave_as_BC(E,N,i,dx,dy,dz,dt,c,sigma,init);
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



    loBound=0,0,0;
    curl_3d(E,Ecurl,dx,dy,dz);

    upBound=N-2,N-1,N-1;
    rd=blitz::RectDomain<3>(loBound,upBound);
    Hadv->x(rd)=Da*H->x(rd)-Db->x(rd)*Ecurl->x(rd);


    upBound=N-1,N-2,N-1;
    rd=blitz::RectDomain<3>(loBound,upBound);
    Hadv->y(rd)=Da*H->y(rd)-Db->y(rd)*Ecurl->y(rd);


    upBound=N-1,N-1,N-2;
    rd=blitz::RectDomain<3>(loBound,upBound);
    Hadv->z(rd)=Da*H->z(rd)-Db->z(rd)*Ecurl->z(rd);



    // Swap H field
    temp=Hadv;
    Hadv=H;
    H=temp;



    curl_3d(H,Hcurl,dx,dy,dz);

    upBound=N-1,N-2,N-2;
    rd=blitz::RectDomain<3>(loBound,upBound);
    Eadv->x(rd)=Ca*E->x(rd)+Cb->x(rd)*Hcurl->x(rd);


    upBound=N-2,N-1,N-2;
    rd=blitz::RectDomain<3>(loBound,upBound);
    Eadv->y(rd)=Ca*E->y(rd)+Cb->y(rd)*Hcurl->y(rd);


    upBound=N-2,N-2,N-1;
    rd=blitz::RectDomain<3>(loBound,upBound);
    Eadv->z(rd)=Ca*E->z(rd)+Cb->z(rd)*Hcurl->z(rd);



    // Swap E field
    temp=Eadv;
    Eadv=E;
    E=temp;



    upBound=N-2,N-2,N-2;
    rd=blitz::RectDomain<3>(loBound,upBound);
    Etotal->z(rd)=E->z(rd)*E->z(rd)+E->y(rd)*E->y(rd)+E->x(rd)*E->x(rd);



    // E->z(N/2,N/2,N/2)=1.0;
    // Eadv->z(N/2,N/2,N/2)=1.0;
    // E->z(N/2,N/2,N/2)=exp(-pow(-init+i,2)/sigma);
    // Eadv->z(N/2,N/2,N/2)=exp(-pow(-init+i,2)/sigma);
    plane_wave_as_BC(E,N,i,dx,dy,dz,dt,c,sigma,init);
    plane_wave_as_BC(Eadv,N,i,dx,dy,dz,dt,c,sigma,init);




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




    //Writing to file
    if (i%3==0&&i>=150&&i<=480){
      dims[0] = E->z.shape()[0];
    	dims[1] = E->z.shape()[1];
      dims[2] = E->z.shape()[2];
    	DataSpace dataspace(3, dims);
      set_name<<"tz"<<i;
    	DataSet dataset = file.createDataSet(set_name.str(), PredType::NATIVE_DOUBLE, dataspace);
      dataset.write(Etotal->z.data(), PredType::NATIVE_DOUBLE);
      set_name.str("");
    }


    // std::cout<<E->z(blitz::Range(45,55),blitz::Range(45,55),N/2)<<std::endl;

    std::cout<<i<<", max="<<blitz::max(Etotal->z)<<std::endl;
    // std::cout<<E->z(rd)<<std::endl;
  }
  std::cout<<"max="<<blitz::max(E->z)<<std::endl;



  return 0;
}
