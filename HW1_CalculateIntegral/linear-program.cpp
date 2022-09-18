#include "linear-program.h"

void LinearProgram(double(*function)(double x), integral_params integ_params)
{
	std::cout << "----------- Linear program start -----------\n";

	auto start = high_resolution_clock::now();

	double integral_result = CalcIntegral(function, integ_params.a, integ_params.b, integ_params.steps_num);

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	std::cout << "Integral answer = " << integral_result << std::endl;
	std::cout << "Total elapsed time = " << duration.count()/1000000. << std::endl;
	std::cout << "----------- Linear program end -----------\n\n\n\n";

}
