#include <stdio.h>
#include <math.h>
#include "zeta.h"

int main(int argc, char **argv) {
	double pi = M_PI;
	double num_pi;

	unsigned long int n;
	int k_max = 24;
	double err_vec[k_max];

	for (int k = 0; k < k_max; k++) {
		n = pow(2, k+1);
		num_pi = sqrt(6*zeta(n));
		err_vec[k] = fabs(num_pi - pi);
	}

	// save content to file
	FILE *fp = fopen("vtest.dat", "w");
	for (int k = 0; k < k_max-1; k++) {
		fprintf(fp, "%.16g ", err_vec[k]);
	}
	fprintf(fp, "%.16g\n", err_vec[k_max-1]);
	fclose(fp);

	return 0;
}