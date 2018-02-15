#include <omp.h>

int main() {
	int a = omp_get_num_threads();
	return 0;
}
