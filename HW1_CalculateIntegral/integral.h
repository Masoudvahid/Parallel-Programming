#ifndef INTEGRAL_H 
#define INTEGRAL_H

struct integral_params {
	double a;
	double b;
	double steps_num;
};

double CalcIntegral(double(*function)(double x), double a, double b, double steps_num);


#endif /* INTEGRAL_H */