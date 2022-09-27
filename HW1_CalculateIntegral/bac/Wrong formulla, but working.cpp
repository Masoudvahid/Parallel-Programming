#include <mpi.h>
#include <stdio.h>
#include <cmath>

double CalcIntegral(double(*function)(double x), double a, double b, int steps_num) {
	double step = (b - a) / steps_num;
	double area = 0.0;
	for (int step_i = 0; step_i < steps_num; step_i++) {
		area += function(a + (step_i + 0.5) * step) * step;
	}
	return area;
}

double MyFunction(double x) {
	return 4 / (1 + (x * x));
}

struct integral_params {
	double a;
	double b;
	double steps_num;
};

int main(int argc, char** argv) {
	double begin_time;
	double end_time;


	integral_params integ_params;
	integ_params.a = 0;
	integ_params.b = 1;
	integ_params.steps_num = 1000000;
	double integral_result = 0;

	int my_rank, world_size;

	MPI_Init(&argc, &argv);
	begin_time = MPI_Wtime();

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	integral_params splitted_params;
	if (my_rank == 0) {
		for (double i = 1; i < world_size; i++)
		{
			double jump_coeff = (i - 1.) / (world_size - 1.);
			splitted_params.a         =         integ_params.a + jump_coeff;
			splitted_params.b         =         integ_params.b - jump_coeff;
			splitted_params.steps_num = integ_params.steps_num / (world_size - 1.);

			MPI_Send(&splitted_params.a,		  1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			MPI_Send(&splitted_params.b,		  1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
			MPI_Send(&splitted_params.steps_num, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
		}
	}

	if (my_rank != 0) {
		printf("Rank %d  ->  ", my_rank);
		MPI_Recv(&splitted_params.a,         1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("a = %lf  ---  ", splitted_params.a);
		MPI_Recv(&splitted_params.b,		  1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("b = %lf  ---  ", splitted_params.b);
		MPI_Recv(&splitted_params.steps_num, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("steps_num = %lf\n\n", splitted_params.steps_num);
		integral_result += CalcIntegral(MyFunction, splitted_params.a, splitted_params.b, splitted_params.steps_num);
		//printf("%lf\n", integral_result);
	}

	 end_time = MPI_Wtime();



	MPI_Finalize();

	if (my_rank == 1) {
		printf("Integral answer = %lf\n", integral_result);
		printf("Total elapsed time = %lf\n", end_time - begin_time);

	}

}