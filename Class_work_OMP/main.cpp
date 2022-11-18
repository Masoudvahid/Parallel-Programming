#include<stdio.h>
#include<stdlib.h>
#include<omp.h>		// the OpenMP header file

int main(int argc, char* argv[]) {
#ifdef _OPENMP
	printf("OpenMP is supported! %d\n", _OPENMP);
#endif
	int a[10];
	int i = 0;
	int myid, num_procs, num_threads;

	num_procs = omp_get_num_procs(); // getting the number of available computing cores
	printf("Num of processors = %d\n", num_procs);

	omp_set_num_threads(2); // explicit setting amount of threads

	num_threads = omp_get_num_threads(); // getting the number of working threads

	printf("Number of threads = %d\n", num_threads);

	for (i = 0; i < 10; i++) {
		a[i] = 0;
	}

	myid = omp_get_thread_num(); // getting thread's number 
	printf("Consecutive part 1, myid = %d\n", myid);

#pragma omp parallel shared(a) private (myid, i) // beginning of the parallel part
	{
		myid = omp_get_thread_num();
		printf("Parallel part, myid = %d\n", myid);

		// !!! this is the place for "#pragma omp for"
	}// end of the parallel part

	myid = omp_get_thread_num();
	printf("Consecutive part 2, myid = %d\n", myid);

} // end of main function
