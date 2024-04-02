#include "cstdlib"
#include "cstdio"
#include "cmath"
#include "memory"
#include "util.h"

#define ERR 1e-9

using namespace std;

const double infty = 1e10;

/*
	config_t is the class for a size-configuration with one flexible element.
*/
struct config_t {
	double sum_a, sum_a_sqr, v0, v1, flex_l, flex_r;
	//sum_a: the sum of elements in g, excluding the flexible element.
	//sum_a_sqr: the sum of squares of elements in g, excluding the flexible element.
	//v0: v0(g), excluding the flexible element.
	//v1: v1(g), excluding the flexible element. 
	//flex_l and flex_r are respectively the lower and upper bounds of the flexible element.

	int flex_type;
	/* flex_type: the type of the flexible element. 
		0 means type 0, 
		1 means type 1, 
		2 means type 2, 
		-1 means type-s or type-b.
	 */ 

	config_t(){}
	config_t(double sum_a, double sum_a_sqr, double v0, double v1, int flex_type, double flex_l, double flex_r) :sum_a(sum_a), sum_a_sqr(sum_a_sqr), v0(v0), v1(v1), flex_type(flex_type), flex_l(flex_l), flex_r(flex_r) {}

};


const int nconfigs_sp2 = 6; //number of categories of configurations for the second sub-program

config_t configs_sp2[nconfigs_sp2] = { 
	config_t(0, 0, 0, 0, -1, 0, 2),		//[0, 2], type-s
	config_t(0, 0, 0, 0, 1, 2, 4), 		//[2, 4], type-1
	config_t(0, 0, 0, 0, 2, 4, 8), 		//[4, 8], type-2
	config_t(0, 0, 0, 0, -1, 8, infty),	//[8, infty], type-b
	config_t(2, 4, 0, 0, -1, 0, 2), 	//2, type-s; [0, 2], type-s
	config_t(2, 4, 0, 0, 1, 2, 4)		//2, type-s; [2, 4], type-1
}; 
//the 6 categories of configurations for the second sub-program. In the comments, elements are separated using semi-colons. Type-o means other types.

const int nconfigs_sp3 = 23; //number of categories of configurations for the third sub-program

config_t configs_sp3[nconfigs_sp3] = {
	config_t(0, 0, 0, 0, -1, 0, 1), 		//[0, 1], type-s
	config_t(0, 0, 0, 0, 0, 1, 2), 			//[1, 2], type-0
	config_t(0, 0, 0, 0, 1, 2, 4),			//[2, 4], type-1
	config_t(0, 0, 0, 0, 2, 4, 8),			//[4, 8], type-2
	config_t(0, 0, 0, 0, -1, 8, infty),		//[8, infty], type-b
	config_t(1, 1, 0, 0, -1, 0, 1),			//1, type-s; [0, 1], type-s
	config_t(1, 1, 0, 0, 0, 1, 2),			//1, type-s; [1, 2], type-0
	config_t(1, 1, 0, 0, 1, 2, 4),			//1, type-s; [2, 4], type-1
	config_t(2, 4, 2, 0, -1, 0, 1),			//2, type-0; [0, 1], type-s
	config_t(2, 4, 0, 2, -1, 0, 1), 		//2, type-1; [0, 1], type-s
	config_t(2, 4, 2, 0, 0, 1, 2),			//2, type-0; [1, 2], type-0
	config_t(2, 4, 0, 2, 0, 1, 2),			//2, type-1; [1, 2], type-0
	config_t(2, 4, 2, 0, 1, 2, 4),			//2, type-0; [2, 4], type-1
	config_t(2, 4, 0, 2, 1, 2, 4),			//2, type-1; [2, 4], type-1
	config_t(2, 2, 0, 0, -1, 0, 1),			//1, type-s; 1, type-s; [0, 1], type-s
	config_t(2, 2, 0, 0, 0, 1, 2),			//1, type-s; 1, type-s; [1, 2], type-0
	config_t(2, 2, 0, 0, 1, 2, 4),			//1, type-s; 1, type-s; [2, 4], type-1
	config_t(3, 5, 2, 0, -1, 0, 1),			//1, type-s; 2, type-0; [0, 1], type-s
	config_t(3, 5, 0, 2, -1, 0, 1),			//1, type-s; 2, type-1; [0, 1], type-s
	config_t(3, 5, 2, 0, 0, 1, 2),			//1, type-s; 2, type-0; [1, 2], type-0
	config_t(3, 5, 0, 2, 0, 1, 2),			//1, type-s; 2, type-1; [1, 2], type-0
	config_t(3, 3, 0, 0, -1, 0, 1),			//1, type-s; 1, type-s; 1, type-s; [0, 1], type-s
	config_t(3, 3, 0, 0, 0, 1, 2)			//1, type-s; 1, type-s; 1, type-s; [1, 2], type-0
}; 
// the 23 categories of configurations for the third sub-program.

struct parameters_t {
	double sp2_mu, sp2_lambda1, sp2_lambda2, sp3_mu, sp3_lambda0, sp3_lambda1, sp3_lambda2, alpha;
	//the 8 parameters, sp2 and sp3 standard for the 2nd and 3rd sub-programs respectively. 

	parameters_t(){}

	parameters_t(
		double sp2_mu, double sp2_lambda1, double sp2_lambda2, 
		double sp3_mu, double sp3_lambda0, double sp3_lambda1, double sp3_lambda2, 
		double alpha
	):
		sp2_mu(sp2_mu), sp2_lambda1(sp2_lambda1), sp2_lambda2(sp2_lambda2),
		sp3_mu(sp3_mu), sp3_lambda0(sp3_lambda0), sp3_lambda1(sp3_lambda1), sp3_lambda2(sp3_lambda2), 
		alpha(alpha) 
	{}
};

parameters_t paras[10] = {
	parameters_t(
		2.060000, 0.080000, 0.000000,
		2.072000, 0.314500, 0.336400, 0.082800,
		1.376228),
	parameters_t(
		2.293595, 0.160000, 0.040000,
		2.220500, 0.315400, 0.364200, 0.160400,
		1.370445),
	parameters_t(
		2.458214, 0.220000, 0.160000,
		2.508757, 0.317800, 0.444200, 0.320000,
		1.364426),
	parameters_t(
		2.634649, 0.320000, 0.280000,
		2.720583, 0.322000, 0.528800, 0.453200,
		1.356049),
	parameters_t(
		2.823747, 0.420000, 0.400000,
		2.874152, 0.325500, 0.594200, 0.547200,
		1.349022),
	parameters_t(
		2.998133, 0.480000, 0.480000,
		2.961363, 0.327900, 0.631000, 0.598400,
		1.344238),
	parameters_t(
		3.213319, 0.560000, 0.560000,
		3.150265, 0.329200, 0.688600, 0.675200,
		1.341530),
	parameters_t(
		3.346480, 0.600000, 0.600000,
		3.343231, 0.329600, 0.736200, 0.733600,
		1.340912),
	parameters_t(
		3.412558, 0.600000, 0.600000,
		3.404549, 0.327300, 0.739800, 0.738000,
		1.356413),
	parameters_t(
		3.470883, 0.600000, 0.600000,
		3.507084, 0.300900, 0.732800, 0.732400,
		1.375000)
};//the 8 parameters for each of the 10 cases.

int main(){
	int i, j, k; 
	double Ll, Lr, L, A, B, C;
	for (int i = 0; i < 10; i++){
		Ll = 2.0 * exp(log(2.0)*i/10);
		Lr = 2.0 * exp(log(2.0)*(i+1)/10);

		printf("************* L interval is [%lf, %lf] ***************\n", Ll, Lr);

		printf("checking the first sub-program\n");
		if (paras[i].alpha + ERR <= 1.5 - 0.5 * sqr(2.0/Lr)) {printf("wrong\n"); return 0;}

		printf("checking the second sub-program\n");

		for (j = 0; j < nconfigs_sp2; j++) { 
			for (k = 0; k < 2; k++) {
				if (k == 0) L = Ll; else L = Lr;
				A = 1 - paras[i].alpha;
				B = - paras[i].alpha * configs_sp2[j].sum_a + paras[i].sp2_mu;
				if (configs_sp2[j].flex_type == 1) B -= paras[i].sp2_lambda1;
				if (configs_sp2[j].flex_type == 2) B -= paras[i].sp2_lambda2;
				C = (1 - paras[i].alpha/2) * configs_sp2[j].sum_a_sqr - paras[i].alpha/2 * sqr(configs_sp2[j].sum_a) 
					+ paras[i].sp2_mu * configs_sp2[j].sum_a 
					- paras[i].sp2_lambda1 * configs_sp2[j].v1 + 0.5*sqr(L) 
					- paras[i].sp2_mu * L - 0.5 * sqr(1) + 0.5 * (sqr(paras[i].sp2_lambda1) + sqr(paras[i].sp2_lambda2));
				if (max_quadratic(A, B, C, configs_sp2[j].flex_l, configs_sp2[j].flex_r) > ERR) {
					printf("wrong\n"); return 0;
				}
			}
		}

		printf("checking the third sub-program\n");

		for (j = 0; j < nconfigs_sp3; j++) {
			for (k = 0; k < 2; k++) {
				if (k == 0) L = Ll; else L = Lr;
				A = 1 - paras[i].alpha;
				B = - paras[i].alpha * configs_sp3[j].sum_a + paras[i].sp3_mu;
				if (configs_sp3[j].flex_type == 0) B -= paras[i].sp3_lambda0;
				if (configs_sp3[j].flex_type == 1) B -= paras[i].sp3_lambda1;
				if (configs_sp3[j].flex_type == 2) B -= paras[i].sp3_lambda2;
				C = (1 - paras[i].alpha/2) * configs_sp3[j].sum_a_sqr - paras[i].alpha/2 * sqr(configs_sp3[j].sum_a) 
					+ paras[i].sp3_mu * configs_sp3[j].sum_a 
					- paras[i].sp3_lambda0 * configs_sp3[j].v0 - paras[i].sp3_lambda1 * configs_sp3[j].v1 
					+ 0.5*sqr(L) - paras[i].sp3_mu * L 
					+ 0.5 * (sqr(paras[i].sp3_lambda0) + sqr(paras[i].sp3_lambda1) + sqr(paras[i].sp3_lambda2));
				if (max_quadratic(A, B, C, configs_sp3[j].flex_l, configs_sp3[j].flex_r) > ERR) {
					printf("wrong\n"); return 0;
				}
			}
		}
	}

	printf("correct\n");
	return 0;
}
