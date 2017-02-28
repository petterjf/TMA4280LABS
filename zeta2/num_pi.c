#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "zeta.h"

int main(int argc, char **argv) {
	if (argc != 2) return 1;

	int n = atoi(argv[1]);

	double sum = zeta(n);

	double num_pi = sqrt(6*sum);
	printf("Accuracy: %f\n", fabs(M_PI - num_pi));

	return 0;
}
