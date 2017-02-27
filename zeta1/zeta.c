double zeta(int l, int u) {
	double sum = 0;
	for (int i = l; i <= u; i++) {
		sum += 1.0/(i*i);
	}

	return sum;
}