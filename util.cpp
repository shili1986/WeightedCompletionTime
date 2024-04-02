double inner_product(int len, double *x, double *y) {
	double res = 0;
	for (int i = 0; i < len; i++) res += x[i] * y[i];
	return res;
}

double sqr(double t){
	return t * t;
}

//Return the maximum value of Ax^2 + Bx + C over x in [xl, xr].
double max_quadratic(double A, double B, double C, double xl, double xr) {
	double mid = - B / (2 * A), value, ret;
	ret = A * sqr(xl) + B * xl + C;
	value = A * sqr(xr) + B * xr + C;
	if (value > ret) ret = value;
	if (xl < mid && mid < xr){
		value = A * sqr(mid) + B * mid + C;
		if (value > ret) ret = value;
	}
	return ret;
}
