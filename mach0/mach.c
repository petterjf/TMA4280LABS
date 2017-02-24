#include <math.h>

double mach(double x, unsigned long int n) {
	double sum = 0;
	for (unsigned long int i = 1; i <= n; i++) {
		sum += pow(-1, i-1)*pow(x, 2*i-1)/(2*i-1); 
	}
	return sum;
}