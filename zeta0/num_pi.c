#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "zeta.h"

int main(int argc, char **argv) {
	if (argc != 2) return 1;

	unsigned long int n = atoi(argv[1]);
	double sum = zeta(n);

	double num_pi = sqrt(6*sum);
	printf("%f\n", num_pi);

	return 0;
}
