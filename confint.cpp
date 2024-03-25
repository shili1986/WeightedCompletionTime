#include "confint.h"
#include "util.h"

#define ERR 1e-9

const double rho = 2;
const double t[5] = {1, rho, rho * rho, rho * rho * rho, 1e10};
double Ll, Lr;

void confint_t::init(int nst0, int nst1, int n0t1, int n1t1, int flexw, int flexl, int flexr){
	this->nst0 = nst0;
	this->nst1 = nst1;
	this->n0t1 = n0t1;
	this->n1t1 = n1t1;
	this->flexw = flexw;
	this->flexl = flexl;
	this->flexr = flexr;

	if (flexl == -1) flexlvalue = 0; else flexlvalue = t[flexl];
	flexrvalue = t[flexr];

	sump = nst0 * t[0] + (nst1 + n0t1 + n1t1) * t[1];
	sumps = nst0 * sqr(t[0]) + (nst1 + n0t1 + n1t1) * sqr(t[1]);

	v[0] = n0t1 * t[1]; v[1] = n1t1 * t[1]; v[2] = 0;
}

bool confint_t::check_alpha(double alpha, double mu, double *lambda, double save) {
	double A = 1-alpha;
	double B = -alpha * sump + mu;
	if (flexw >= 0) B -= lambda[flexw];
	double C, L;
	for (int i = 0; i < 2; i++){
		if (i == 0) L = Ll; else L = Lr;
		C = (1 - alpha/2) * sumps + 0.5 * sqr(L) - alpha/2 * sqr(sump) + mu * (sump - L) + 0.5 * inner_product(3, lambda, lambda) - inner_product(3, lambda, v) - save;

		if (max_quadratic(A, B, C, flexlvalue, flexrvalue) > ERR) return false;
	}
	return true;
}

void confint_t::print_to_file(FILE *fp){
	int i;
	for (i = 0; i < nst0; i++) fprintf(fp, "t0, none; ");
	for (i = 0; i < nst1; i++) fprintf(fp, "t1, none; ");
	for (i = 0; i < n0t1; i++) fprintf(fp, "t1, v0; ");
	for (i = 0; i < n1t1; i++) fprintf(fp, "t1, v1; ");
	if (flexl == -1) fprintf(fp, "[0, "); else fprintf(fp, "[t%d, ", flexl);
	if (flexr == 4) fprintf(fp, "infty], "); else fprintf(fp, "t%d], ", flexr);
	if (flexw == -1) fprintf(fp, "none\n"); else fprintf(fp, "v%d\n", flexw);
}

