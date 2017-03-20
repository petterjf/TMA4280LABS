double zeta(unsigned long int l, unsigned long int u) {
	double sum = 0;
	for (unsigned long int i = l; i <= u; i++) {
		sum += 1.0/(i*i);
	}

	return sum;
}