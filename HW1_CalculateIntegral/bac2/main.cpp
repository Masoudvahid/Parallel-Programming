#include <mpi.h>
#include <math.h>
#include <iostream>

double MyFunction(double x) {
	return 4 / (1 + (x * x));
}

struct integral_params {
	double upp_lim;
	double low_lim;
	double iterations;
};

double CalcIntegral(double a, double step_size, double steps_num) {
	double step_sizee = (step_size - a) / steps_num;
	double area = 0.0;
	for (int step_i = 0; step_i < steps_num; step_i++) {
		area += 2 * MyFunction(a + step_i * step_sizee);
	}
	return (area * step_sizee) / 2;
}



int main() {

	integral_params integ_params;
	integ_params.upp_lim = 0;
	integ_params.low_lim = 1;
	integ_params.iterations = 10;

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

	if (my_rank == 0) {

		double step_size = (integ_params.low_lim - integ_params.upp_lim) / (world_size - 1.);
		double step_size_remainder = fmod((integ_params.low_lim - integ_params.upp_lim), (world_size - 1.));
		step_size += step_size_remainder;
		std::cout << "<<<< " << step_size << " >>>>> ";


		double steps_num = integ_params.iterations / (world_size - 1.);
		double steps_num_remainder = fmod(integ_params.iterations, (world_size - 1.));
		steps_num += steps_num_remainder;

		for (double i = 1; i < world_size; i++)
		{
			double lower_limit = integ_params.upp_lim * (i - 1);

			MPI_Send(&lower_limit, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			MPI_Send(&step_size, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
			MPI_Send(&steps_num, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
		}
	}
	
	integral_params splitted_params;
	if (my_rank != 0) {
		std::cout << "Rank " << my_rank << " -> ";
		MPI_Recv(&splitted_params.upp_lim, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		std::cout << "<a = " << splitted_params.upp_lim << "> --- ";
		MPI_Recv(&splitted_params.low_lim, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		std::cout << "<b = " << splitted_params.low_lim << "> --- ";
		MPI_Recv(&splitted_params.iterations, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		std::cout << "<steps_num = " << splitted_params.iterations << ">" << std::endl;

		integral_result = CalcIntegral(splitted_params.upp_lim, splitted_params.low_lim, splitted_params.iterations);
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
		total_elapsed += MPI_Wtime() - begin_time;
		std::cout << "Total elapsed time = " << total_elapsed << std::endl;
	}

	MPI_Finalize();
}