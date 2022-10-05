#ifndef INTEGRAL_H 
#define INTEGRAL_H

struct integral_params {
	double upp_lim;
	double low_lim;
	double iterations;
};

double CalcIntegral(double(*function)(double x), double a, double b, double steps_num);


#endif /* INTEGRAL_H */