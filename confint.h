#include "cstdio"
using namespace std;

extern const double rho, t[5];
extern double Ll, Lr;

class confint_t {
public:
	int nst0, nst1, n0t1, n1t1;
	int flexw, flexl, flexr;
	double flexlvalue, flexrvalue;
	
	double sump, sumps, v[3];
	
	void init(int nst0, int nst1, int n0t1, int n1t1, int flexw, int flexl, int flexr);

	bool check_alpha(double alpha, double mu, double *lambda, double save);

	void print_to_file(FILE *fp);
};


