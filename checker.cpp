#include "cstdlib"
#include "cstdio"
#include "cmath"
#include "memory"

#include "util.h"

#define ERR 1e-9

using namespace std;

struct confint_t {
	double sump, sumps,	v0, v1, l, r;
	int w;

	confint_t(){}
	confint_t(double sump, double sumps, double v0, double v1, int w, double l, double r) :sump(sump), sumps(sumps), v0(v0), v1(v1), w(w), l(l), r(r) {}

};

const double t0 = 1, t1 = 2, t2 = 4, t3 = 8;
const double infty = 1e10;
const int nci0 = 6;

confint_t ci0[nci0] = { 
	confint_t(0, 0, 0, 0, -1, 0, t1),		//[0, t1], none
	confint_t(0, 0, 0, 0, 1, t1, t2), 		//[t1, t2], v1
	confint_t(t1, t1*t1, 0, 0, -1, 0, t1), 	//t1, none; [0, t1], none
	confint_t(t1, t1*t1, 0, 0, 1, t1, t2),	//t1, none; [t1, t2], v1
	confint_t(0, 0, 0, 0, 2, t2, t3), 		//[t2, t3], v2
	confint_t(0, 0, 0, 0, -1, t3, infty)	//[t3, infty], none
};

const int nci1 = 23;
confint_t ci1[nci1] = {
	confint_t(0, 0, 0, 0, -1, 0, t0), 			//[0, t0], none
	confint_t(0, 0, 0, 0, 0, t0, t1), 			//[t0, t1], v0
	confint_t(0, 0, 0, 0, 1, t1, t2),			//[t1, t2], v1
	confint_t(t0, t0*t0, 0, 0, -1, 0, t0),		//t0, none; [0, t0], none
	confint_t(t0, t0*t0, 0, 0, 0, t0, t1),		//t0, none; [t0, t1], v0
	confint_t(t0, t0*t0, 0, 0, 1, t1, t2),		//t0, none; [t1, t2], v1
	confint_t(t0*2, t0*t0*2, 0, 0, -1, 0, t0),	//t0, none; t0, none; [0, t0], none
	confint_t(t0*2, t0*t0*2, 0, 0, 0, t0, t1),	//t0, none; t0, none; [t0, t1], v0
	confint_t(t0*2, t0*t0*2, 0, 0, 1, t1, t2),	//t0, none; t0, none; [t1, t2], v1
	confint_t(t0*3, t0*t0*3, 0, 0, -1, 0, t0),	//t0, none; t0, none; t0, none; [0, t0], none
	confint_t(t0*3, t0*t0*3, 0, 0, 0, t0, t1),	//t0, none; t0, none; t0, none; [t0, t1], v0
	confint_t(t1, t1*t1, t1, 0, -1, 0, t0),		//t1, v0; [0, t0], none
	confint_t(t1, t1*t1, 0, t1, -1, 0, t0), 	//t1, v1; [0, t0], none
	confint_t(t1, t1*t1, t1, 0, 0, t0, t1),		//t1, v0; [t0, t1], v0
	confint_t(t1, t1*t1, 0, t1, 0, t0, t1),		//t1, v1; [t0, t1], v0
	confint_t(t1, t1*t1, t1, 0, 1, t1, t2),		//t1, v0; [t1, t2], v1
	confint_t(t1, t1*t1, 0, t1, 1, t1, t2),		//t1, v1; [t1, t2], v1
	confint_t(t0+t1, t0*t0+t1*t1, t1, 0, -1, 0, t0),	//t0, none; t1, v0; [0, t0], none
	confint_t(t0+t1, t0*t0+t1*t1, 0, t1, -1, 0, t0),	//t0, none; t1, v1; [0, t0], none
	confint_t(t0+t1, t0*t0+t1*t1, t1, 0, 0, t0, t1),	//t0, none; t1, v0; [t0, t1], v0
	confint_t(t0+t1, t0*t0+t1*t1, 0, t1, 0, t0, t1),	//t0, none; t1, v1; [t0, t1], v0
	confint_t(0, 0, 0, 0, 2, t2, t3),					//[t2, t3], v2
	confint_t(0, 0, 0, 0, -1, t3, infty)				//[t3, infty], none
};

double parameters[10][8] = {
	{2.060000, 0.080000, 0.000000,
	2.072000, 0.314500, 0.336400, 0.082800,
	1.376228},
	{2.293595, 0.160000, 0.040000,
	2.220500, 0.315400, 0.364200, 0.160400,
	1.370445},
	{2.458214, 0.220000, 0.160000,
	2.508757, 0.317800, 0.444200, 0.320000,
	1.364426},
	{2.634649, 0.320000, 0.280000,
	2.720583, 0.322000, 0.528800, 0.453200,
	1.356049},
	{2.823747, 0.420000, 0.400000,
	2.874152, 0.325500, 0.594200, 0.547200,
	1.349022},
	{2.998133, 0.480000, 0.480000,
	2.961363, 0.327900, 0.631000, 0.598400,
	1.344238},
	{3.213319, 0.560000, 0.560000,
	3.150265, 0.329200, 0.688600, 0.675200,
	1.341530},
	{3.346480, 0.600000, 0.600000,
	3.343231, 0.329600, 0.736200, 0.733600,
	1.340912},
	{3.412558, 0.600000, 0.600000,
	3.404549, 0.327300, 0.739800, 0.738000,
	1.356413},
	{3.470883, 0.600000, 0.600000,
	3.507084, 0.300900, 0.732800, 0.732400,
	1.375000}
};

int main(){
	int i, j, k; 
	double Ll, Lr, L, A, B, C;
	double alpha, mu, lambda0, lambda1, lambda2;
	for (int i = 0; i < 10; i++){
		Ll = t1 * exp(log(2)*i/10);
		Lr = t1 * exp(log(2)*(i+1)/10);
		printf("************* L interval is [%lf, %lf] ***************\n", Ll, Lr);

		alpha = parameters[i][7];

		printf("checking the first subprogram\n");
		if (alpha + ERR <= 1.5 - 0.5 * sqr(t1/Lr)) {printf("wrong\n"); return 0;}

		printf("checking the second subprogram\n");
		mu = parameters[i][0]; lambda1 = parameters[i][1]; lambda2 = parameters[i][2]; 

		for (j = 0; j < nci0; j++) { 
			for (k = 0; k < 2; k++) {
				if (k == 0) L = Ll; else L = Lr;
				A = 1 - alpha;
				B = - alpha * ci0[j].sump + mu;
				if (ci0[j].w == 1) B -= lambda1;
				if (ci0[j].w == 2) B -= lambda2;
				C = (1 - alpha/2) * ci0[j].sumps - alpha/2 * sqr(ci0[j].sump) + mu * ci0[j].sump 
					- lambda1 * ci0[j].v1 + 0.5*sqr(L) 
					- mu * L - 0.5 * sqr(t0) + 0.5 * (sqr(lambda1) + sqr(lambda2));
				if (max_quadratic(A, B, C, ci0[j].l, ci0[j].r) > ERR) {printf("wrong\n"); return 0;}
			}
		}

		printf("checking the third subprogram\n");
		mu = parameters[i][3]; lambda0 = parameters[i][4];
		lambda1 = parameters[i][5]; lambda2 = parameters[i][6]; 

		for (j = 0; j < nci1; j++) {
			for (k = 0; k < 2; k++) {
				if (k == 0) L = Ll; else L = Lr;
				A = 1 - alpha;
				B = - alpha * ci1[j].sump + mu;
				if (ci1[j].w == 0) B -= lambda0;
				if (ci1[j].w == 1) B -= lambda1;
				if (ci1[j].w == 2) B -= lambda2;
				C = (1 - alpha/2) * ci1[j].sumps - alpha/2 * sqr(ci1[j].sump) + mu * ci1[j].sump 
					- lambda0 * ci1[j].v0 - lambda1 * ci1[j].v1 
					+ 0.5*sqr(L) - mu * L + 0.5 * (sqr(lambda0) + sqr(lambda1) + sqr(lambda2));
				if (max_quadratic(A, B, C, ci1[j].l, ci1[j].r) > ERR) {printf("wrong\n"); return 0;}
			}
		}
	}

	printf("correct\n");
	return 0;
}
