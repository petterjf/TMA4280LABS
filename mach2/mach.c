#include <math.h>
#include <omp.h>

double mach(double x, int n) {
	double sum = 0;
	#pragma omp parallel for schedule(static)
	for (int i = 1; i <= n; i++) {
		sum += powl(-1, i-1)*pow(x, 2*i-1)/(2*i-1);
	}
	return sum;
}