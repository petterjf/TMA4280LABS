#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mach.h"

int main(int argc, char **argv) {
	if (argc != 2) return 1;

	unsigned long int n = atoi(argv[1]);
	double num_pi = 4*(4*mach(1.0/5.0, n) - mach(1.0/239.0, n));
	printf("Accuracy: %16g\n", fabs(M_PI - num_pi));
	return 0;
}