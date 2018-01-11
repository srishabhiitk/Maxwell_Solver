#include <iostream>
#include <blitz/array.h>
#include <math.h>
// #include "H5Cpp.h"
// using namespace H5;


// using namespace std;


class test{
public:
  double data;

  test(double data_){
    data=data_;
  }

  void set_to_zero(){
    data=0;
  }
};


void test_fn(blitz::Array<double,1> h, blitz::firstIndex x){
  h=x/2.0;
}



int main()
{
    // cout << "Thank you for running my first C++ program!" << endl;
    // blitz::Array<int,2> a(3,3),b(3,3),sum(3,3),F(3,3);
    bool first=true,second=true;
    blitz::Array<int,2> a,b,sum;
    blitz::Array<double,3> F;
    blitz::Array<int,2> c(3,3);
    blitz::TinyVector<int, 3> tiny1;
    blitz::TinyVector<int, 3> tiny2;

    a.resize(3,3);
    b.resize(3,3);
    c.resize(0,0);
    // c(3,3);
    sum.resize(3,3);
    F.resize(3,3,3);
    F=0;
    a=0;
    a=1,2,3,4,5,6,7,8,9;
    b=2,3,4,5,6,7,1,2,3;
    std::cout << "a:" << std::endl;
    std::cout << a(blitz::Range(0,2),blitz::Range(0,2)) << std::endl;
    std::cout << "b:" << std::endl;
    std::cout << b(blitz::Range(0,2),blitz::Range(0,2)) << std::endl;
    b(blitz::Range(0,2),blitz::Range(0,2))=a(blitz::Range(0,1),blitz::Range(0,1));
    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::thirdIndex k;
    sum=b/a;
    F(blitz::Range(0,2),blitz::Range(0,2),blitz::Range(0,2))=i+j+k;
    double d=M_PI*4*pow(10,-7);
    test var(d);
    std::cout<< "Before setting to zero: "<< var.data <<std::endl;
    var.set_to_zero();
    std::cout<< "After setting to zero: "<< var.data <<std::endl;
    // int a, b, sum;

    blitz::Array<double,1> h(3);
    // b.resize(3,1);
    h=0;
    test_fn(h,i);
    std::cout << "F_old:" << std::endl;
    std::cout << F(blitz::Range(0,2),blitz::Range(0,2),blitz::Range(0,2)) << std::endl;
    F=F(i,j,k)*h(i);
    F=(F<3)*F;
    std::cout << "F_new:" << std::endl;
    std::cout << F(blitz::Range(0,2),blitz::Range(0,2),blitz::Range(0,2)) << std::endl;
    // cin >> a;
    // cin >> b;
    // a=2;
    // b=3;
    // sum = a + b;
    blitz::TinyVector<int,2> lb,ub;
    lb=0,0;
    ub=2,2;
    blitz::RectDomain<2> rd;
    rd= blitz::RectDomain<2>(lb,ub);
    // rd.ubound()(0)+=1;

    std::cout << "a:" << std::endl;
    std::cout << a(rd) << std::endl;
    std::cout << "h:" << std::endl;
    std::cout << h << std::endl;
    std::cout << "sum:" << std::endl;
    std::cout << sum(rd) << std::endl;
    std::cout << "F:" << std::endl;
    std::cout << F(blitz::Range(0,2),0,0) << std::endl;

    tiny1(0) = 0;
    tiny1(1) = 1;
    tiny1(2) = 2;
    tiny2 = 22;
    tiny2(0) = tiny2(0) -1;
    int foo[3] = {tiny1(0), tiny1(1), tiny1(2)};
    std::cout<<"tiny1: "<<tiny1<<std::endl;
    std::cout<<"tiny2: "<<tiny2<<std::endl;
    std::cout<<"foo: "<<foo[2]<<std::endl;

    lb=1,0;
    ub=ub-lb;
    rd= blitz::RectDomain<2>(lb,ub);
    lb=-1,3;
    ub=0,0;
    bool chk=lb(0)<0;

    std::cout<<chk<<std::endl;



    // F=F*0.33;
    // H5File file("try.h5", H5F_ACC_TRUNC);
  	// hsize_t dims[3];               // dataset dimensions
  	// dims[0] = 3;
  	// dims[1] = 3;
    // dims[2] = 1;
  	// DataSpace dataspace(3, dims);

  	// Create the dataset.
  	// DataSet dataset = file.createDataSet("t0", PredType::NATIVE_DOUBLE, dataspace);
    // dataset.write(F.data(), PredType::NATIVE_DOUBLE);



    return 0;
}
