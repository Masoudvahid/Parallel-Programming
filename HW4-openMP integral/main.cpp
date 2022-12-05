#include <iostream>
#include <stdlib.h>
#include <omp.h> // OpenMP library
#include <bits/stdc++.h>

using std::cin;
using std::cout;
using std::endl;

void SaveResults(double speedup);
double Function(double x);

const int THREADS_NUM = 16;
double intgral_result = 0; // Result of integral as a shared memory variable
double N = 1'000'000;      // Number of steps

int upper_limit = 1; // Upper limit of integral
int lower_limit = 0; // Lower limit of integral

double dx = (upper_limit - lower_limit) / N; // Step size

int main(int argc, char **argv)
{
    /* -------- Sequential code -------- */
    clock_t sqe_start, sqe_end;
    sqe_start = clock();

    intgral_result += Function(upper_limit) + Function(lower_limit);

    for (int i = 0; i < N; i++)
    {
        intgral_result += 2 * Function(i * dx);
    }

    sqe_end = clock();
    double seq_elapsed = double(sqe_end - sqe_start) / double(CLOCKS_PER_SEC);

/* -------- Parallel code -------- */
#ifdef _OPENMP
    cout << "OPENMP is supported! -ver: " << _OPENMP << endl;
#endif

    intgral_result = 0;
    double integral_result_reduction = 0;

    double parallel_start = omp_get_wtime();

    // intgral_result += Function(upper_limit) + Function(lower_limit); // Adding the first and last parts of the

    int myid, num_procs, num_threads;
    num_procs = omp_get_max_procs(); //   getting the number of available computing cores
    printf("Num of processors = %d \n", num_procs);

    omp_set_num_threads(THREADS_NUM); // explicit setting amount of threads

    num_threads = omp_get_num_threads(); //   getting the number of working threads
    printf("Number of threads  = %d \n", num_threads);

    myid = omp_get_thread_num(); //   getting thread’s number

#pragma omp parallel shared(intgral_result, integral_result_reduction) private(myid) //   beginning of the parallel part
    {
        myid = omp_get_thread_num();
        double local_result = 0;
#pragma omp for
        for (int i = 0; i < (int)N; i++)
        {
            local_result += 2 * Function(i * dx);

        } //   the loop counter i must be private according to the logic, but if it is shared, then it is implicitly converted to private
#pragma omp for reduction(+:integral_result_reduction)
        for (int i = 0; i < (int)N; i++)
        {
            integral_result_reduction += 2 * Function(i * dx);
        }
        
#pragma omp critical
        cout << "integral calculated by thread <" << myid << "> is = " << local_result << endl;
        intgral_result += local_result;
    }

    cout << "\n\nIntegral result = " << (dx / 2) * intgral_result << endl;
    cout << "Integral result using reductions = " << (dx / 2) * integral_result_reduction << endl;

    double parallel_elapsed = omp_get_wtime() - parallel_start;
    double speedup = seq_elapsed / parallel_elapsed;
    cout << "Speedup = " << speedup << endl;

    SaveResults(speedup);
}

double Function(double x)
{
    return 4 / (1 + x * x);
}

void SaveResults(double speedup)
{
    std::ofstream output;
    std::string file_name = "P=" + std::to_string(THREADS_NUM) + " - " + std::to_string(speedup) + ".txt";
    output.open(file_name);
    output << THREADS_NUM << "\n"
           << speedup;
}
