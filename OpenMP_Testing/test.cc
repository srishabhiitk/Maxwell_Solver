#include <cstdlib> 
#include <iostream>
#include <sys/time.h>
using namespace std;
int main() 
{
	int N=5000;
	int n_threads=4;
	int i,j;
	int **arr1 = (int**)malloc(N * sizeof(int *));
	for(int i = 0; i < N; i++)
		arr1[i] = (int*)malloc(N * sizeof(int));
	int **arr2 = (int**)malloc(N * sizeof(int *));
	for(int i = 0; i < N; i++)
		arr2[i] = (int*)malloc(N * sizeof(int));
	int **arr3 = (int**)malloc(N * sizeof(int *));
	for(int i = 0; i < N; i++)
		arr3[i] = (int*)malloc(N * sizeof(int));
	struct timeval begin,end;
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){ 
			arr1[i][j] = rand(); 
			arr2[i][j] = rand();
		}
	}
	gettimeofday(&begin,NULL);
	#pragma omp parallel for num_threads(n_threads) default(none) shared(arr1,arr2,arr3,N)
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){ 
			arr3[i][j] = arr1[i][j] + arr2[i][j];
			arr3[i][j] = arr1[i][j] - arr2[i][j];
			arr3[i][j] = arr1[i][j] * arr2[i][j];
		}
	}
	gettimeofday(&end,NULL);
	std::cout<<"Time taken: "<<((end.tv_sec-begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.0e6<<" sec"<<std::endl;	
}
