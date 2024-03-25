double inner_product(int len, double *x, double *y) {
	double res = 0;
	for (int i = 0; i < len; i++) res += x[i] * y[i];
	return res;
}

double sqr(double t){
	return t * t;
}

double max_quadratic(double A, double B, double C, double xl, double xr) {
	double mid = - B / (2 * A), xstar;
	if (mid < xl)
		xstar = xl;
	else if (mid > xr)
		xstar = xr;
	else
		xstar = mid;
	return A * sqr(xstar) + B * xstar + C;
}
