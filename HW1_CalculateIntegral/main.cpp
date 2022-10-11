#include <mpi.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>


double MyFunction(double x) {
	return 4 / (1 + (x * x));
}

double CalcIntegral(double low_lim, double step_size, int steps_num) {
	double area = 0.;
	for (int step_i = 0; step_i < steps_num; step_i++) {
		area += 2 * MyFunction(low_lim + step_i * step_size);
	}
	return (area * step_size) / 2;
}

void PrintOutput(int my_rank, double low_lim, double step_size, int steps_num) {
	std::cout << "Rank " << my_rank << " -> ";
	std::cout << "<low_lim = " << low_lim << "> --- ";
	std::cout << "<step_size = " << step_size << "> --- ";
	std::cout << "<steps_num = " << steps_num << ">" << std::endl;
}

void SaveResults(int world_size, double sub_intervals, double total_elapsed) {
	std::ofstream output;
	std::string file_name = "P=" + std::to_string(world_size) + " - N=" + std::to_string((int)sub_intervals) + ".txt";
	output.open(file_name);
	output << world_size << "\n" << sub_intervals << "\n" << total_elapsed;
}

int main() {
	double low_lim = 0;
	double up_lim = 1;
	double sub_intervals = 10;


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

	double step_size = (up_lim - low_lim) / sub_intervals;
	int steps_residue = (int)sub_intervals % world_size;

	if (my_rank == 0) {
		begin_time = MPI_Wtime();
		int steps_num = (int)sub_intervals / world_size;
		PrintOutput(my_rank, low_lim, step_size, steps_num);

		final_answer += CalcIntegral(low_lim, step_size, steps_num);
	}

	if (my_rank != 0) {
		if (my_rank <= steps_residue) {
			int steps_num = (int)sub_intervals / world_size + 1;
			low_lim += ((up_lim - low_lim) / sub_intervals) * ((int)sub_intervals / world_size);
			low_lim += (my_rank - 1) * step_size * steps_num;
			PrintOutput(my_rank, low_lim, step_size, steps_num);

			integral_result = CalcIntegral(low_lim, step_size, steps_num);
		}
		else {
			int steps_num = (int)sub_intervals / world_size;
			low_lim += ((up_lim - low_lim) / sub_intervals) * ((int)sub_intervals / world_size);

			low_lim += steps_residue != 0 ? (steps_residue)*step_size * (steps_num + 1) : (my_rank - 1) * step_size * (steps_num);
			PrintOutput(my_rank, low_lim, step_size, steps_num);

			integral_result = CalcIntegral(low_lim, step_size, steps_num);
		}

		MPI_Send(&integral_result, 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
	}

	if (my_rank == 0) {
		for (int i = 1; i < world_size; i++)
		{
			MPI_Recv(&integral_result, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			final_answer += integral_result;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (my_rank == 0) {
		std::cout << "Integral answer = " << final_answer << std::endl;
		fflush(stdout);
		total_elapsed += (MPI_Wtime() - begin_time);
		std::cout << "Time elapsed for < " << world_size << " > nodes = " << total_elapsed << std::endl;
		fflush(stdout);

		SaveResults(world_size, sub_intervals, total_elapsed);
	}

	MPI_Finalize();
}