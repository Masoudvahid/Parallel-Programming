#include "parallel.h"

void ParallelProgram(double(*function)(double x), integral_params integ_params) {
	MPI_Init(NULL, NULL);

	double integral_result = 0;
	double final_answer = 0;

	MPI_Barrier(MPI_COMM_WORLD);
	double begin_time = MPI_Wtime();
	double total_elapsed = 0;

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	integral_params splitted_params;
	if (my_rank == 0) {
		//std::cout << "----------- Parallel program start -----------\n";
		for (double i = 1; i < world_size; i++)
		{
			double jump_coeff = (integ_params.low_lim - integ_params.upp_lim) / (world_size - 1.);
			splitted_params.upp_lim = jump_coeff * (i - 1);
			splitted_params.low_lim = jump_coeff * i;
			splitted_params.iterations = integ_params.iterations / (world_size - 1.);

			MPI_Send(&splitted_params.upp_lim, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			MPI_Send(&splitted_params.low_lim, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
			MPI_Send(&splitted_params.iterations, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
		}
	}

	if (my_rank != 0) {
		std::cout << "Rank " << my_rank << " -> ";
		MPI_Recv(&splitted_params.upp_lim, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		std::cout << "a = " << splitted_params.upp_lim << " --- ";
		MPI_Recv(&splitted_params.low_lim, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		std::cout << "b = " << splitted_params.low_lim << " --- ";
		MPI_Recv(&splitted_params.iterations, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		std::cout << "steps_num = " << splitted_params.iterations << std::endl;

		integral_result = CalcIntegral(function, splitted_params.upp_lim, splitted_params.low_lim, splitted_params.iterations);
		MPI_Send(&integral_result, 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);

		double end_time = MPI_Wtime() - begin_time;
		MPI_Send(&end_time, 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
	}

	if (my_rank == 0) {
		for (size_t i = 1; i < world_size; i++)
		{
			MPI_Recv(&integral_result, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			final_answer += integral_result;

			double elapsed = 0;
			MPI_Recv(&elapsed, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			total_elapsed += elapsed;
		}
		
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (my_rank == 0) {
		std::cout << "Integral answer = " << final_answer << std::endl;
		total_elapsed += MPI_Wtime() - begin_time; /////////// inserted without debug !!!!
		std::cout << "Total elapsed time = " << total_elapsed << std::endl;
	}

	MPI_Finalize();
	// QUESTIONS: why MPI_Barrier doesn  not work -> i want to print start and end of the program?
	// QUESTIONS: why mpi does not run between init and finalize?
	// QUESTIONS: why wrong formulation works correctly
	// QUESTIONS: why mpi runs the whole project multiple times
	// QUESTIONS:
	// QUESTIONS:
}