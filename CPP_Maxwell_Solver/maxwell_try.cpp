#include <iostream>
#include <blitz/array.h>
#include <math.h>
#include <fstream>



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
      Fcurl->y(rd)=(F->x(shift(rd,2,1))-F->x(rd))/dy-(F->z(shift(rd,0,1))-F->y(rd))/dz;

      upBound=N-3,N-3,N-2;
      rd=blitz::RectDomain<3>(loBound,upBound);
      Fcurl->z(rd)=(F->y(shift(rd,0,1))-F->y(rd))/dy-(F->x(shift(rd,1,1))-F->y(rd))/dz;
    }
    else{
      upBound=N-2,N-2,N-2;

      loBound=0,1,1;
      rd=blitz::RectDomain<3>(loBound,upBound);
      Fcurl->x(rd)=(F->z(rd)-F->z(shift(rd,1,-1)))/dy-(F->y(rd)-F->y(shift(rd,2,-1)))/dz;

      loBound=1,0,1;
      rd=blitz::RectDomain<3>(loBound,upBound);
      Fcurl->y(rd)=(F->x(rd)-F->x(shift(rd,2,-1)))/dy-(F->y(rd)-F->z(shift(rd,0,-1)))/dz;

      loBound=1,1,0;
      rd=blitz::RectDomain<3>(loBound,upBound);
      Fcurl->z(rd)=(F->y(rd)-F->y(shift(rd,0,-1)))/dy-(F->y(rd)-F->x(shift(rd,1,-1)))/dz;
    }
  }
  //Error if requirements of curl_3d are not met!
  else{
    if (!(F->is2D)){
      std::cout<<"Error! curl_3d is defined for 3D fields only."<<std::endl;
    }
    else{
      std::cout<<"First two parameter should be staggered and unstaggered."<<std::endl;
    }

  }
}



using blitz::Range;

int main(){
  int num_timesteps=100;
  int N=100;
  double epsilon0=8.85*pow(10,-12);
  double mew0=M_PI*4*pow(10,-7);

  double c=1/pow(epsilon0*mew0,0.5); //Speed of light
  double S=1/pow(2,0.5); //Courant factor
  double dx=pow(10,-3);
  double dy=dx, dz=dx;
  double dt=S*dx/c;
  double Ca=1.0,Da=1.0;
  double Cb=dt/epsilon0;
  double Db=dt/mew0;





  field *E=new field(N,false,false);
  E->intialize_to_zero();
  field *Eadv=new field(N,false,false);
  Eadv->intialize_to_zero();
  field *Ecurl=new field(N,true,false);
  Ecurl->intialize_to_zero();


  field *H=new field(N,true,false);
  H->intialize_to_zero();
  field *Hadv=new field(N,true,false);
  Hadv->intialize_to_zero();
  field *Hcurl=new field(N,false,false);
  Hcurl->intialize_to_zero();


  field *temp;//for swapping



  blitz::RectDomain<3> rd;
  blitz::TinyVector<int,3> loBound,upBound;

  int i;
  for (i=0; i<num_timesteps; i++){
    if (i==0){
      E->z(N/2,N/2,N/2)=1;
    }


    loBound=0,0,0;
    curl_3d(E,Ecurl,dx,dy,dz);

    upBound=N-2,N-1,N-1;
    rd=blitz::RectDomain<3>(loBound,upBound);
    Hadv->x(rd)=Da*H->x(rd)-Db*Ecurl->x(rd);

    upBound=N-1,N-2,N-1;
    rd=blitz::RectDomain<3>(loBound,upBound);
    Hadv->y(rd)=Da*H->y(rd)-Db*Ecurl->y(rd);

    upBound=N-1,N-1,N-2;
    rd=blitz::RectDomain<3>(loBound,upBound);
    Hadv->z(rd)=Da*H->z(rd)-Db*Ecurl->z(rd);


    // Swap H field
    temp=Hadv;
    Hadv=H;
    H=temp;



    curl_3d(H,Hcurl,dx,dy,dz);

    upBound=N-1,N-2,N-2;
    rd=blitz::RectDomain<3>(loBound,upBound);
    Eadv->x(rd)=Ca*E->x(rd)+Cb*Hcurl->x(rd);

    upBound=N-2,N-1,N-2;
    rd=blitz::RectDomain<3>(loBound,upBound);
    Eadv->y(rd)=Ca*E->y(rd)+Cb*Hcurl->y(rd);

    upBound=N-2,N-2,N-1;
    rd=blitz::RectDomain<3>(loBound,upBound);
    Eadv->z(rd)=Ca*E->z(rd)+Cb*Hcurl->z(rd);


    // Swap E field
    temp=Eadv;
    Eadv=E;
    E=temp;




    E->z(N/2,N/2,N/2)=1;
    Eadv->z(N/2,N/2,N/2)=1;

    std::cout<<i<<std::endl;
  }



  // // Print result to file
  // std::ofstream myfile;
  // myfile.open("Ez.data");
  // int j;
  //
  //
  // for (i=0;i<N;i++){
  //   for (j=0;j<N;j++)
  //     myfile<<E->z(i,j,0)<<" ";
  //   myfile<<"\n";
  // }
  std::cout<<"max="<<blitz::max(E->z)<<std::endl;




  return 0;
}
