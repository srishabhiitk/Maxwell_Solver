#include <omp.h>
#include <sys/time.h>
#include <iostream>
static long num_steps = 100000000; double step;
#define NUM_THREADS 6
int main(){
	int i; double pi,sum=0.0;
	struct timeval begin,end;
	step = 1.0/(double) num_steps;
	omp_set_num_threads(NUM_THREADS);
	gettimeofday(&begin,NULL);
	#pragma omp parallel
	{
		double x;
		#pragma omp for  reduction (+:sum)
		for (i=0; i<num_steps; i++){
			x = (i+0.5)*step;
			sum += 4.0/(1.0+x*x);
		}
	}
	pi += sum*step;
	//for (i=0, pi=0.0; i<nthreads; i++) pi += sum[i]*step;
	gettimeofday(&end,NULL);
	std::cout<<"Time taken: "<<((end.tv_sec-begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.0e6<<" sec"<<std::endl;
	return 0;
}
