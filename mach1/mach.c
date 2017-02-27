#include <math.h>

double mach(double x, int l, int u) {
	double sum = 0;
	for (int i = l; i <= u; i++) {
		sum += pow(-1, i-1)*pow(x, 2*i-1)/(2*i-1); 
	}
	return sum;
}
