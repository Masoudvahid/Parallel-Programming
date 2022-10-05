#ifndef LINEAR_PROGRAM_H 
#define LINEAR_PROGRAM_H

#include <iostream>
#include <chrono>
#include "integral.h"
using namespace std::chrono;

void LinearProgram(double(*function)(double x), integral_params integ_params);


#endif /* LINEAR_PROGRAM_H */