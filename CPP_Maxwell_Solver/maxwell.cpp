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
        x.resize(N,N,N);
        y.resize(N,N,N);
        z.resize(N,N,N);
      }
    }
  }

  void intialize_to_zero(){
    x=0;
    y=0;
    z=0;
  }
};






using blitz::Range;






void plane_wave_as_BC(field *E,int n,int timestep,double dx,double dt,double c){
  int i,j;
  for (j=0;j<n;j++){
    i=0;
    E->z(i,j,0)=exp(-pow(((i*dx+j*dx)/pow(2,0.5)+300*dx-timestep*c*dt),2)/pow(10,-2));
    i=n-1;
    E->z(i,j,0)=exp(-pow(((i*dx+j*dx)/pow(2,0.5)+300*dx-timestep*c*dt),2)/pow(10,-2));
  }
  for (i=1;i<(n-1);i++){
    j=0;
    E->z(i,j,0)=exp(-pow(((i*dx+j*dx)/pow(2,0.5)+300*dx-timestep*c*dt),2)/pow(10,-2));
    j=n-1;
    E->z(i,j,0)=exp(-pow(((i*dx+j*dx)/pow(2,0.5)+300*dx-timestep*c*dt),2)/pow(10,-2));
  }
}





int main(){
  int num_timesteps=2800;
  int N=1000;
  double epsilon0=8.85*pow(10,-12);
  double mew0=M_PI*4*pow(10,-7);

  double c=1/pow(epsilon0*mew0,0.5); //Speed of light
  double S=1/pow(2,0.5); //Courant factor
  double dx=pow(10,-3);
  double dt=S*dx/c;
  double Ca=1.0,Da=1.0;
  double Cb=dt/epsilon0/dx;
  double Db=dt/mew0/dx;





  field *E=new field(N,false,true);
  E->intialize_to_zero();
  field *E1=new field(N,false,true);
  E1->intialize_to_zero();


  field *H=new field(N,true,true);
  H->intialize_to_zero();
  field *H1=new field(N,true,true);
  H1->intialize_to_zero();


  field *temp;


  /***************************************************************************/
  // Timestep loop
  /***************************************************************************/
  int i;
  for (i=0; i<num_timesteps; i++){
    // if (i==0){
    //   E->z(N/2,N/2,0)=1;
    // }

    plane_wave_as_BC(E,N,i,dx,dt,c);



    H1->y(Range(0,N-2),Range(0,N-1),0)=Da*H->y(Range(0,N-2),Range(0,N-1),0)+Db*(E->z(Range(1,N-1),Range(0,N-1),0)-E->z(Range(0,N-2),Range(0,N-1),0));
    H1->x(Range(0,N-1),Range(0,N-2),0)=Da*H->x(Range(0,N-1),Range(0,N-2),0)-Db*(E->z(Range(0,N-1),Range(1,N-1),0)-E->z(Range(0,N-1),Range(0,N-2),0));

    // Swap H field
    temp=H1;
    H1=H;
    H=temp;


    E1->z(Range(1,N-2),Range(1,N-2),0)=Ca*E->z(Range(1,N-2),Range(1,N-2),0)+Cb*(H->y(Range(1,N-2),Range(1,N-2),0)-H->y(Range(0,N-3),Range(1,N-2),0)+H->x(Range(1,N-2),Range(0,N-3),0)-H->x(Range(1,N-2),Range(1,N-2),0));

    // Swap E field
    temp=E1;
    E1=E;
    E=temp;



    //Extra conditions
    E1->z(Range(N/2-50,N/2+50),Range(N/2-50,N/2+50),0)=0;
    E->z(Range(N/2-50,N/2+50),Range(N/2-50,N/2+50),0)=0;
    // E->z(N/2,N/2,0)=1;
    // E1->z(N/2,N/2,0)=1;

    std::cout<<i<<std::endl;
  }



  // Print result to file
  std::ofstream myfile;
  myfile.open("Ez.data");
  int j;


  for (i=0;i<N;i++){
    for (j=0;j<N;j++)
      myfile<<E->z(i,j,0)<<" ";
    myfile<<"\n";
  }
  std::cout<<blitz::max(E->z)<<std::endl;




  return 0;
}
