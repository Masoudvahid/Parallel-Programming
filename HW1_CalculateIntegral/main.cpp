#include <mpi.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>


double MyFunction(double x) {
	return 4 / (1 + (x * x));
}

struct integral_params {
	double low_lim;
	double step_size;
	double steps_num;
};

struct integral_info {
	double low_lim;
	double up_lim;
	double sub_intervals;
};

double CalcIntegral(integral_params integ_params) {
	double area = 0.;
	for (int step_i = 0; step_i < integ_params.steps_num; step_i++) {
		area += 2 * MyFunction(integ_params.low_lim + step_i * integ_params.step_size);
	}
	return (area * integ_params.step_size) / 2;
}

void PrintOutput(int my_rank, integral_params integ_params) {
	std::cout << "Rank " << my_rank << " -> ";
	std::cout << "<low_lim = " << integ_params.low_lim << "> --- ";
	std::cout << "<step_size = " << integ_params.step_size << "> --- ";
	std::cout << "<steps_num = " << integ_params.steps_num << ">" << std::endl;
}

void SaveResults(int world_size, double sub_intervals, double total_elapsed) {
	std::ofstream output;
	std::string file_name = "P=" + std::to_string(world_size) + " - N=" + std::to_string((int)sub_intervals) + ".txt";
	output.open(file_name);
	output << world_size << "\n" << sub_intervals << "\n" << total_elapsed;
}

int main() {
	integral_info integ_info;
	integ_info.low_lim = 0;
	integ_info.up_lim = 1;
	integ_info.sub_intervals = 100000000;

	integral_params integ_params{};

	MPI_Init(NULL, NULL);

	double integral_result = 0;
	double final_answer = 0;

	MPI_Barrier(MPI_COMM_WORLD);
	double begin_time = 0;
	double total_elapsed = 0;

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	double jump = (int)integ_info.sub_intervals / world_size;
	int jump_residue = (int)integ_info.sub_intervals % world_size;

	double new_step_size = (integ_info.up_lim - integ_info.low_lim) / integ_info.sub_intervals;
	if (my_rank == 0) {
		begin_time = MPI_Wtime();
		integ_params.steps_num = jump;
		integ_params.step_size = new_step_size;
		PrintOutput(my_rank, integ_params);

		final_answer += CalcIntegral(integ_params);
		integ_params.low_lim += jump * new_step_size;

	}

	if (my_rank != 0) {
		integral_params splitted_params;
		double new_lower_limit = 0;
		if (my_rank < jump_residue) {
			if (my_rank == 1) {
				double new_lower_limit = integ_params.low_lim + (jump * new_step_size);
			}
			else
			{
				double new_lower_limit = integ_params.low_lim + (my_rank * new_step_size);
			}
			double new_step_size = (integ_info.up_lim - integ_info.low_lim) / integ_info.sub_intervals;

			integ_params.steps_num = jump + 1;
			integ_params.step_size = (jump + 1) * new_step_size;
		}
		else
		{

		}

		
		splitted_params.low_lim = new_lower_limit;



		integral_result = CalcIntegral(splitted_params);
		MPI_Send(&integral_result, 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);

		//double end_time = MPI_Wtime() - begin_time;
		//MPI_Send(&end_time, 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
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
		total_elapsed += (MPI_Wtime() - begin_time) / world_size;
		std::cout << "Time elapsed for < " << world_size << " > nodes = " << total_elapsed << std::endl;

		SaveResults(world_size, integ_info.sub_intervals, total_elapsed);
	}

	MPI_Finalize();



}