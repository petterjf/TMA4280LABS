#include <stdio.h>
#include <math.h>
#include "zeta.h"

int main(int argc, char **argv) {
	if (argc != 1) return 1;

	double sum = zeta(3); 
	double e_sum = 49.0/36.0; // expected sum

	printf("%.16g\n", sum - e_sum);

	return 0;
}
