#include "parallel.h"

void ParallelProgram(double(*function)(double x), integral_params integ_params) {
	MPI_Init(NULL, NULL);

	double integral_result = 0;
	double final_answer = 0;

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
			double jump_coeff = (integ_params.b - integ_params.a) / (world_size - 1.);
			splitted_params.a = jump_coeff * (i - 1);
			splitted_params.b = jump_coeff * i;
			splitted_params.steps_num = integ_params.steps_num / (world_size - 1.);

			MPI_Send(&splitted_params.a, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			MPI_Send(&splitted_params.b, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
			MPI_Send(&splitted_params.steps_num, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD); // TODO: Why this line is not imortant?

	if (my_rank != 0) {
		std::cout << "Rank " << my_rank << " -> ";
		MPI_Recv(&splitted_params.a, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		std::cout << "a = " << splitted_params.a << " --- ";
		MPI_Recv(&splitted_params.b, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		std::cout << "b = " << splitted_params.b << " --- ";
		MPI_Recv(&splitted_params.steps_num, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		std::cout << "steps_num = " << splitted_params.steps_num << std::endl;

		integral_result = CalcIntegral(function, splitted_params.a, splitted_params.b, splitted_params.steps_num);
		MPI_Send(&integral_result, 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);

		double end_time = MPI_Wtime() - begin_time;
		MPI_Send(&end_time, 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD); // TODO: Why this line is not imortant?

	if (my_rank == 0) {
		for (size_t i = 1; i < world_size; i++)
		{
			MPI_Recv(&integral_result, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			final_answer += integral_result;

			double elapsed = 0;
			MPI_Recv(&elapsed, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			total_elapsed += elapsed;
		}
		total_elapsed += MPI_Wtime() - begin_time; /////////// inserted without debug !!!!
		std::cout << "Integral answer = " << final_answer << std::endl;
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