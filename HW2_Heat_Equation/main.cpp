#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <mpi.h>

#define leftSideTemp 0.0
#define rightSideTemp 0.0

double time_tot = 0.1;          // total time 
double dt = 0.0002;             // time step
// double h = 0.02;                // step on the x-coordinate  /* double h = l / nPoints; */
int nPoints = 50000;
double l = 1.0;                 // length of the rod
double h = l / nPoints;                // step on the x-coordinate
double u_0 = 1.0;               // initial value
double pi = 3.14159265358;

void SaveResults(int world_size, double speedup, double nPoints) {
	std::ofstream output;
	std::string file_name = "P=" + std::to_string(world_size) + " - np=" + std::to_string(nPoints) + ".txt";
	output.open(file_name);
	output << world_size << "\n" << nPoints << "\n" << speedup;
}

int main(int argc, char** argv)
{
	MPI_Init(NULL, NULL);

	int beginIterator, endIterator, world_size, my_rank;

	int countOfTimePieces = time_tot / dt;

	double* u_prev = (double*)malloc(sizeof(double) * (nPoints + 2));  // first array for the numerical solution
	double* u_next = (double*)malloc(sizeof(double) * (nPoints + 2));  // first array for the numerical solution
	double* u_exact = (double*)malloc(sizeof(double) * (nPoints + 2)); // array for the exact solution

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	MPI_Barrier(MPI_COMM_WORLD);
	double begin_time = 0;
	double total_elapsed = 0;

	/* Set initial conditions. */
	for (int i = 1; i < nPoints + 1; i++)
	{
		u_prev[i] = u_next[i] = u_0;
	}

	/* Set the boundary conditions. */
	u_prev[0] = u_next[0] = leftSideTemp;
	u_prev[nPoints + 1] = u_next[nPoints + 1] = rightSideTemp;

	/***** Sequential program start *****/
	MPI_Barrier(MPI_COMM_WORLD);
	if (my_rank == 0) {
		begin_time = MPI_Wtime(); // Start the timer for sequential program

		double time = 0;
		while (time < time_tot) {
			for (int i = 1; i <= nPoints; i++) {
				// finite difference scheme
				u_next[i] = u_prev[i] + dt / (h * h) * (u_prev[i + 1] - 2 * u_prev[i] + u_prev[i - 1]);
			}
			for (int i = 1; i <= nPoints; ++i) {
				u_prev[i] = u_next[i];
			}
			time = time + dt;
		}
		printf("Numerical solution with sequential programming: \n");
		for (int i = 1; i <= nPoints; i++) {
			printf("%f ", u_next[i]);
		}
		printf("\n");
		printf("\n");

		//Reset all conditions
		for (int i = 1; i <= nPoints; i++) {
			u_prev[i] = u_0;
			u_next[i] = u_0;
		}
		u_prev[0] = 0;
		u_prev[nPoints + 1] = 0;
		u_next[0] = 0;
		u_next[nPoints + 1] = 0;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	double sequential_time = MPI_Wtime() - begin_time;
	/***** Sequential program end *****/

	/***** Parallel program start *****/
	int k = nPoints / world_size;
	int remainder = nPoints % world_size;
	int left, right, amount;

	// Set the left and right index for every processes
	MPI_Barrier(MPI_COMM_WORLD);
	begin_time = MPI_Wtime();
	if (my_rank < remainder) {
		left = my_rank * k + my_rank + 1;
		right = left + k;
		amount = k + 1;
	}
	else {
		left = my_rank * k + remainder + 1;
		right = left + k - 1;
		amount = k;
	}

	double time = 0;
	while (time < time_tot) {
		if (my_rank % 2 != 0) {
			if (my_rank > 0)
				MPI_Send(&u_prev[left], 1, MPI_DOUBLE, my_rank - 1, 1, MPI_COMM_WORLD);
			if (my_rank > 0)
				MPI_Recv(&u_prev[left - 1], 1, MPI_DOUBLE, my_rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (my_rank < world_size - 1)
				MPI_Recv(&u_prev[right + 1], 1, MPI_DOUBLE, my_rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (my_rank < world_size - 1)
				MPI_Send(&u_prev[right], 1, MPI_DOUBLE, my_rank + 1, 1, MPI_COMM_WORLD);
		}
		else {
			if (my_rank < world_size - 1)
				MPI_Recv(&u_prev[right + 1], 1, MPI_DOUBLE, my_rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (my_rank < world_size - 1)
				MPI_Send(&u_prev[right], 1, MPI_DOUBLE, my_rank + 1, 1, MPI_COMM_WORLD);
			if (my_rank > 0)
				MPI_Send(&u_prev[left], 1, MPI_DOUBLE, my_rank - 1, 1, MPI_COMM_WORLD);
			if (my_rank > 0)
				MPI_Recv(&u_prev[left - 1], 1, MPI_DOUBLE, my_rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		for (int i = left; i <= right; i++) {
			u_next[i] = u_prev[i] + dt / (h * h) * (u_prev[i + 1] - 2 * u_prev[i] + u_prev[i - 1]);
		}
		for (int i = left; i <= right; i++) {
			u_prev[i] = u_next[i];
		}
		time = time + dt;
	}

	// Gather the results using process 0
	if (my_rank == 0) {
		for (int i = 1; i < world_size; ++i) {
			if (i < remainder) {
				left = i * k + i + 1;
				amount = k + 1;
			}
			else {
				left = i * k + remainder + 1;
				amount = k;
			}
			MPI_Recv(&u_next[left], amount, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		// printing the numerical solution on the screen
		printf("Numerical solution with parallel programming: \n");
		for (int i = 1; i <= nPoints; i++) {
			printf("%f ", u_next[i]);
		}
		printf("\n");
		printf("\n");

	}
	else {
		MPI_Send(&u_next[left], amount, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	double parallel_time = MPI_Wtime() - begin_time;

	if (my_rank == 0) {
		// Print the time calculating
		printf("sequetial programming elapsed time: %f\n", sequential_time);
		printf("parallel programming elapsed time: %f\n", parallel_time);
		double speed_up = sequential_time / parallel_time;
		printf("Speed up: %f\n", speed_up);
		printf("Efficiency: %f\n", speed_up / world_size);
		
		// exact solution
		printf("Exact solution: \n");
		for (int i = 1; i <= nPoints; i++) {
			double x = i * h;
			double sum = 0;
			for (int j = 0; j < 5; j++) {
				double a = exp(-pi * pi * (2 * j + 1) * (2 * j + 1) * time_tot) * sin(pi * (2 * j + 1) * x / l) / (2 * j + 1);
				sum = sum + 4 * u_0 * a / pi;
			}
			u_exact[i] = sum;
			// printing the exact solution on the screen 
			printf("%f  ", u_exact[i]);
		}
		printf("\n");

		/***** Save output *****/
		SaveResults(world_size, speed_up, nPoints);
	}
	/***** Parallel program end *****/



	/*
	if (world_size > 1)
	{
		int mainLenghForOneProcces = (nPoints - 1) / (world_size - 1);

		int residualLenghForOneProcces = (nPoints - 1) % (world_size - 1);

		if (my_rank != 0)
		{
			// Procces rank isn't 0

			if ((my_rank - 1) < residualLenghForOneProcces)
			{
				beginIterator = (my_rank - 1) * (mainLenghForOneProcces + 1) + 1;
				endIterator = beginIterator + mainLenghForOneProcces + 1;
			}
			else
			{
				beginIterator = (my_rank - 1) * mainLenghForOneProcces + 1 + residualLenghForOneProcces;
				endIterator = beginIterator + mainLenghForOneProcces;
			}

			if (my_rank == world_size - 1)
			{
				endIterator = nPoints;
			}

			// Integrate over time.
			for (int i = 0; i < countOfTimePieces; i++)
			{

				if ((my_rank % 2 == 1) &&
					(my_rank != world_size - 1) &&
					(my_rank != 1))
				{
					MPI_Send(&u_prev[endIterator - 1], 1, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD);
					MPI_Recv(&u_prev[endIterator], 1, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

					MPI_Recv(&u_prev[beginIterator - 1], 1, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Send(&u_prev[beginIterator], 1, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD);
				}
				else
				{
					if ((my_rank % 2 == 0) &&
						(my_rank != world_size - 1) &&
						(my_rank != 1))
					{
						MPI_Recv(&u_prev[beginIterator - 1], 1, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						MPI_Send(&u_prev[beginIterator], 1, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD);

						MPI_Send(&u_prev[endIterator - 1], 1, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD);
						MPI_Recv(&u_prev[endIterator], 1, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					}
					else
					{
						if (my_rank == 1)
						{
							MPI_Send(&u_prev[endIterator - 1], 1, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD);
							MPI_Recv(&u_prev[endIterator], 1, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						}
						else
						{
							if (my_rank == world_size - 1)
							{
								MPI_Recv(&u_prev[beginIterator - 1], 1, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Send(&u_prev[beginIterator], 1, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD);
							}
						}
					}
				}
				for (int j = beginIterator; j < endIterator; j++)
				{
					u_next[j] = u_prev[j] + dt / h / h * (u_prev[j - 1] - 2.0 * u_prev[j] + u_prev[j + 1]);
				}

				u_exact = u_prev;
				u_prev = u_next;
				u_next = u_exact;
			}

			MPI_Send(u_next + beginIterator, endIterator - beginIterator, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
		else
		{
			//Data collection.
			// Procces rank is 0
			int i;

			for (i = 1; i < world_size; i++)
			{
				if ((i - 1) < residualLenghForOneProcces)
				{
					beginIterator = (i - 1) * mainLenghForOneProcces + i;
					endIterator = beginIterator + mainLenghForOneProcces + 1;
				}
				else
				{
					beginIterator = (i - 1) * mainLenghForOneProcces + 1 + residualLenghForOneProcces;
					endIterator = beginIterator + mainLenghForOneProcces + 1;
				}

				if (i == world_size - 1)
				{
					endIterator = nPoints;
				}

				MPI_Recv(u_next + beginIterator, endIterator - beginIterator, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}

			// Output on display.
			for (i = 0; i < nPoints + 1; i++)
			{
				printf("%lf\n", u_next[i]);
			}

			total_elapsed += (MPI_Wtime() - begin_time);
			printf("Time elapsed for < %ld > = %lf \n", world_size, total_elapsed);
			fflush(stdout);

		}
	}
	else
	{
		// Integrate over time.
		for (int i = 0; i < countOfTimePieces; i++)
		{
			for (int i = 1; i < nPoints; i++)
			{
				u_next[i] = u_prev[i] + dt / h / h * (u_prev[i - 1] - 2.0 * u_prev[i] + u_prev[i + 1]);
			}

			u_exact = u_prev;
			u_prev = u_next;
			u_next = u_exact;
		}

		Output on display.
		for (int i = 0; i < nPoints + 1; i++) {
			printf("%lf \n", u_next[i]);
		}

		total_elapsed += (MPI_Wtime() - begin_time);
		printf("Time elapsed for < %ld > = %lf \n", world_size, total_elapsed);
		fflush(stdout);

	}

	*/

	free(u_prev);
	free(u_next);

	MPI_Finalize();

	return 0;
}