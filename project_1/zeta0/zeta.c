double zeta(unsigned long int n) {
	double sum = 0;
	for (unsigned long int i = 1; i <= n; i++) {
		sum += 1.0/(i*i);
	}

	return sum;
}