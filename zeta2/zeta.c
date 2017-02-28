#include <omp.h>

double zeta(int n) {
	double sum = 0;
	#pragma omp parallel for schedule(static)
	for (int i = 1; i <= n; i++) {
		sum += 1.0/(i*i);
	}

	return sum;
}