#include <math.h>
#include "parallel.h"
#include "linear-program.h"

//#define run_linear

double MyFunction(double x) {
	return 4 / (1 + (x * x));
}



int main(int argc, char** argv) {
	integral_params integ_params;
	integ_params.a = 0;
	integ_params.b = 1;
	integ_params.steps_num = 10;

#ifdef run_linear
	LinearProgram(MyFunction, integ_params);
#else // run_linear
	ParallelProgram(MyFunction, integ_params);
#endif

}