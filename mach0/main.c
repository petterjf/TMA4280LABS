#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double num_arctan(double, int);

int main(int argc, char **argv) {
	if (argc != 2) return 1;

	int n = atoi(argv[1]);
	double num_pi = 4*(4*num_arctan(1.0/5.0, n) - num_arctan(1.0/239.0, n));
	printf("%f\n", num_pi);

	return 0;
}

double num_arctan(double x, int n) {
	double sum = 0;
	for (int i = 1; i <= n; i++) {
		sum += pow(-1, i-1)*pow(x, 2*i-1)/(2*i-1); 
	}
	return sum;
}