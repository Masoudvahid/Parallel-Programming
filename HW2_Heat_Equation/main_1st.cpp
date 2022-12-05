#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

using namespace std;

double time_tot = 0.1; // total time
double dt = 0.00002;   // time step
double h = 1;          // step on the x-coordinate
double l = 1;          // length of the rod
double u_0 = 1;        // initial value
double pi = 3.14159265358;

int main(int argc, char *argv[])
{
    int i, m;
    double time, x, sum, a;
    int nPoints = 11;
    h = l / (nPoints + 1);
    //    dt = 0.25*h*h;
    double u_prev[nPoints + 2];  // first array for the numerical solution
    double u_next[nPoints + 2];  // second array for the numerical solution
    double u_exact[nPoints + 2]; // array for the exact solution

    int myrank;
    int size;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Initial conditions
    for (i = 1; i <= nPoints; i++)
    {
        u_prev[i] = u_0;
        u_next[i] = u_0;
    }

    // Boundary conditions
    u_prev[0] = 0;
    u_prev[nPoints + 1] = 0;
    u_next[0] = 0;
    u_next[nPoints + 1] = 0;

    // main loop on time
    time = 0;
    while (time < time_tot)
    {
        for (i = 1; i <= nPoints; i++)
        {
            // finite difference scheme
            u_next[i] = u_prev[i] + dt / (h * h) * (u_prev[i + 1] - 2 * u_prev[i] + u_prev[i - 1]);
        }
        for (i = 1; i <= nPoints; i++)
        {
            u_prev[i] = u_next[i];
        }
        time = time + dt;
    }

    double begin, end, total;

    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0)
    {
        begin = MPI_Wtime();
    }

    double *my_local = (double *)malloc(sizeof(double) * (nPoints + 2));
    double *my_local_temp = (double *)malloc(sizeof(double) * (nPoints + 2));

    int k, n;
    for (k = 1; k <= nPoints; k++)
    {
        my_local[k] = my_local_temp[k] = 1;
    }

    my_local[0] = my_local_temp[0] = 0;
    my_local[nPoints + 1] = my_local_temp[nPoints + 1] = 0;

    int *cutter = (int *)malloc(sizeof(int) * (size + 1));
    cutter[0] = 1;
    cutter[size] = nPoints + 1;
    for (int i = 1; i < size; i++)
    {
        cutter[i] = cutter[i - 1] + (nPoints / size) + ((i - 1) < ((nPoints % size) - 2));
    }

    // main loop on time
    int N = time_tot / dt;
    for (n = 0; n < N; n++)
    {
        for (k = cutter[myrank]; k < cutter[myrank + 1]; k++)
        {
            // finite difference scheme
            my_local_temp[k] = my_local[k] + dt / (h * h) * (my_local[k + 1] - 2 * my_local[k] + my_local[k - 1]);
        }

        if (myrank < size - 1)
        {
            MPI_Send(my_local + cutter[myrank + 1] - 1, 1, MPI_DOUBLE, myrank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(my_local + cutter[myrank + 1], 1, MPI_DOUBLE, myrank + 1, 0, MPI_COMM_WORLD, &status);
        }

        if (myrank > 0)
        {
            MPI_Recv(my_local + cutter[myrank] - 1, 1, MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD, &status);
            MPI_Send(my_local + cutter[myrank], 1, MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD);
        }

        double *t = my_local;
        my_local = my_local_temp;
        my_local_temp = t;
    }

    if (size > 1)
    {
        if (myrank == 0)
        {
            for (int i = 1; i < size; i++)
            {
                MPI_Recv(my_local_temp + cutter[i], cutter[i + 1] - cutter[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
            }
        }
        else
        {
            MPI_Send(my_local_temp + cutter[myrank], cutter[myrank + 1] - cutter[myrank], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0)
    {
        end = MPI_Wtime();
        total = end - begin;
        printf("\nTotal time is %.5f\n", total);
    }

    if (myrank == 0)
    {
        for (int i = 1; i <= nPoints; i++)
        {
            printf("%f  ", my_local_temp[i]);
        }
        printf("\n");
        printf("\n");
    }

    if (myrank == 0)
    {
        // exact solution
        printf("Exact solution: \n");
        for (int i = 1; i <= nPoints; i++)
        {
            x = i * h;
            sum = 0;
            for (int m = 0; m < 5; m++)
            {
                a = exp(-pi * pi * (2 * m + 1) * (2 * m + 1) * time_tot) * sin(pi * (2 * m + 1) * x / l) / (2 * m + 1);
                sum = sum + 4 * u_0 * a / pi;
            }
            u_exact[i] = sum;
            // printing the exact solution on the screen
            printf("%f  ", u_exact[i]);
        }
    }

    free(my_local);
    free(my_local_temp);

    MPI_Finalize();
    return 0;
}