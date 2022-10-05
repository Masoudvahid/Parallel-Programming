#ifndef PARALLEL_H 
#define PARALLEL_H

#include "mpi.h"
#include <iostream>
#include "integral.h"

void ParallelProgram(double(*function)(double x), integral_params integ_params);


#endif /* PARALLEL_H */