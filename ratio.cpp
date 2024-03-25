#include "cstdlib"
#include "cstdio"
#include "cmath"
#include "memory"

#include "util.h"
#include "confint.h"

#define ERR 1e-9
#define C_FILE "configurations.txt"
#define P_FILE "parameters.txt"

using namespace std;

const double nintervals = 10;

double min_alpha(int nconfints, confint_t *confints, double mu, double *lambda, double save) {

	double alphal = 1, alphar = 1.5;

	while (alphar - alphal > 1e-5) {
		double alpha = (alphal + alphar) / 2;
		bool OK = true;
		for (int i = 0; i < nconfints; i++)
			if (!confints[i].check_alpha(alpha, mu, lambda, save)) {
				OK = false;
				break;
			}
		if (OK) alphar = alpha; else alphal = alpha;
	}
	return alphar;
}

double search_parameters(int nconfints, confint_t *confints, double mul, double mur, double mustep, double *lambdal, double *lambdar, double *lambdastep, double save, double &bestmu, double *bestlambda) {
	double alpha = 1.51, mu, lambda[3];
	for (mu = mul; mu < mur + ERR; mu += mustep)
		for (lambda[0] = lambdal[0]; lambda[0] < lambdar[0] + ERR; lambda[0] += lambdastep[0])
			for (lambda[1] = lambdal[1]; lambda[1] < lambdar[1] + ERR; lambda[1] += lambdastep[1])
				for (lambda[2] = lambdal[2]; lambda[2] < lambdar[2] + ERR; lambda[2] += lambdastep[2]){
					double a = min_alpha(nconfints, confints, mu, lambda, save);
					if (a < alpha) {
						alpha = a; 
						bestmu = mu;
						memcpy(bestlambda, lambda, sizeof(lambda[0]) * 3);
					}
					
				}
	return alpha;
}

/*void test(){
	L = 1;
	int nconfints = 0;
	confint_t confints[20];
	confints[nconfints++].init(0, 0, 0, 1, t[1], t[2]);

	double mu = 0.5, lambda[3] = {0, 0, 0};
	printf("min_alpha = %lf\n", min_alpha(nconfints, confints, mu, lambda, sqr(t[0])/2));
}*/

int main(){

//	test(); exit(0);


	int nconfints0, nconfints1;
	confint_t confints0[10], confints1[1000];

	nconfints0 = 0;
	confints0[nconfints0++].init(0, 0, 0, 0, -1, -1, 1);
	confints0[nconfints0++].init(0, 0, 0, 0, 1, 1, 2);
	confints0[nconfints0++].init(0, 1, 0, 0, -1, -1, 1);
	confints0[nconfints0++].init(0, 1, 0, 0, 1, 1, 2);
	confints0[nconfints0++].init(0, 0, 0, 0, 2, 2, 3);
	confints0[nconfints0++].init(0, 0, 0, 0, -1, 3, 4);

	nconfints1 = 0;
	for (int nt1 = 0; t[1] * nt1 < t[2]; nt1++)
		for (int nst0 = 0; t[0] * nst0 + t[1] * nt1 < t[2]; nst0++) {
			confints1[nconfints1++].init(nst0, 0, nt1, 0, -1, -1, 0);
			if (nt1 > 0) confints1[nconfints1++].init(nst0, 0, 0, nt1, -1, -1, 0);
			confints1[nconfints1++].init(nst0, 0, nt1, 0, 0, 0, 1);
			if (nt1 > 0) confints1[nconfints1++].init(nst0, 0, 0, nt1, 0, 0, 1);
			if (nst0 == 0 || t[0] * (nst0 - 1) + t[1] * (nt1 + 1) < t[2]) {
				confints1[nconfints1++].init(nst0, 0, nt1, 0, 1, 1, 2);
				if (nt1 > 0) confints1[nconfints1++].init(nst0, 0, 0, nt1, 1, 1, 2);
			}
		}
	confints1[nconfints1++].init(0, 0, 0, 0, 2, 2, 3);
	confints1[nconfints1++].init(0, 0, 0, 0, -1, 3, 4);
	
	FILE *fp = fopen(C_FILE, "w");
	fprintf(fp, "%d\n", nconfints0);
	for (int i = 0; i < nconfints0; i++) 
		confints0[i].print_to_file(fp);
	fprintf(fp, "%d\n", nconfints1);
	for (int i = 0; i < nconfints1; i++) 
		confints1[i].print_to_file(fp);
	fclose(fp);

	fp = fopen(P_FILE, "w");
	double avg = 0;

	for (int iL = 0; iL < nintervals; iL++){
		Ll = t[1] * exp(iL * log (rho)/nintervals);
		Lr = t[1] * exp((iL + 1) * log (rho)/nintervals);
		printf("**************    L interval = [%lf, %lf]    ***************\n", Ll, Lr);

		double bestmu, bestlambda[3], finalalpha, alpha;

		// interval 1 is full
		finalalpha = 1.5 - sqr(t[1]/Lr)/2;
		printf("interval 1 full: alpha = %lf\n", finalalpha);

		// interval 0 is full, interval 1 is not full

		double lambdal[3] = {0, 0, 0}, lambdar[3] = {0, 0.5 * t[1], 0.5 * t[2]};
		double lambdastep[3] = {0.1, 0.01 * t[1], 0.01 * t[2]};
		alpha = search_parameters(nconfints0, confints0, Ll/2, Lr, 0.01 * Ll, lambdal, lambdar, lambdastep, sqr(t[0])/2, bestmu, bestlambda);
		if (alpha > finalalpha) finalalpha = alpha;

		printf("interval 0 full interval 1 not full: alpha = %lf\n", alpha);
		fprintf(fp, "%lf %lf %lf\n", bestmu, bestlambda[1], bestlambda[2]);

		// intervals 0 and 1 are not full
		
		// coarse search
		lambdal[0] = 0; lambdal[1] = 0; lambdal[2] = 0; 
		lambdar[0] = 0.5 * t[0]; lambdar[1] = 0.5 * t[1]; lambdar[2] = 0.5 * t[2];
		lambdastep[0] = 0.1 * t[0]; lambdastep[1] = 0.1 * t[1]; lambdastep[2] = 0.1 * t[2];
		alpha = search_parameters(nconfints1, confints1, Ll/2, Lr, 0.1 * Ll, lambdal, lambdar, lambdastep, 0, bestmu, bestlambda);

		// fine-grained search
		for (int i = 0; i < 3; i++){
			lambdal[i] = bestlambda[i] - lambdastep[i];
			lambdar[i] = bestlambda[i] + lambdastep[i];
		}
		lambdastep[0] = 0.01 * t[0]; lambdastep[1] = 0.01 * t[1]; lambdastep[2] = 0.01 * t[2];
		alpha = search_parameters(nconfints1, confints1, bestmu - 0.1 * Ll, bestmu + 0.1 * Ll, 0.01 * Ll, lambdal, lambdar, lambdastep, 0, bestmu, bestlambda);

		for (int i = 0; i < 3; i++){
			lambdal[i] = bestlambda[i] - lambdastep[i];
			lambdar[i] = bestlambda[i] + lambdastep[i];
		}
		lambdastep[0] = 0.001 * t[0]; lambdastep[1] = 0.001 * t[1]; lambdastep[2] = 0.001 * t[2];
		alpha = search_parameters(nconfints1, confints1, bestmu - 0.01 * Ll, bestmu + 0.01 * Ll, 0.001 * Ll, lambdal, lambdar, lambdastep, 0, bestmu, bestlambda);

		for (int i = 0; i < 3; i++){
			lambdal[i] = bestlambda[i] - lambdastep[i];
			lambdar[i] = bestlambda[i] + lambdastep[i];
		}
		lambdastep[0] = 0.0001 * t[0]; lambdastep[1] = 0.0001 * t[1]; lambdastep[2] = 0.0001 * t[2];
		alpha = search_parameters(nconfints1, confints1, bestmu - 0.001 * Ll, bestmu + 0.001 * Ll, 0.0001 * Ll, lambdal, lambdar, lambdastep, 0, bestmu, bestlambda);
		if (alpha > finalalpha) finalalpha = alpha;

		printf("intervals 0 and 1 are not full: alpha = %lf\n", alpha);
		fprintf(fp, "%lf %lf %lf %lf\n", bestmu, bestlambda[0], bestlambda[1], bestlambda[2]);

		printf("final_alpha = %lf\n\n", finalalpha);
		fprintf(fp, "%lf\n", finalalpha); 

		avg += finalalpha;
	}

	avg = avg/nintervals;
	printf("average ratio = %lf\n", avg); 

	fclose(fp);
	return 0;
}

