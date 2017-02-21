#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv) {
	if (argc != 2) return 1;

	int n = atoi(argv[1]);
	double sum = 0;
	for (int i = 1; i <= n; i++) {
		sum += 1.0/(i*i);
	}

	double num_pi = sqrt(6*sum);
	printf("%f\n", num_pi);

	return 0;
}
