#include <stdio.h>
#include <stdlib.h>
#include "mach.h"

int main(int argc, char **argv) {
	if (argc != 1) return 1;

	double sum = mach(1.0/5.0, 3);
	double e_sum = 9253.0/46875.0; // expected mach sum when x = 1/5 and n = 3
	printf("%f\n", sum - e_sum);

	return 0;
}