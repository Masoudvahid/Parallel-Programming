#include <iostream>
#include <cmath>
#include <time.h>
#include <chrono>
#include <pthread.h>
#include <bits/stdc++.h>
#include <semaphore.h>


#define PI 3.14159265359
#define NUM_THREADS 1
double POINTS_NUM = 100000000;

void *MonteCarloEstimate(void *param);
double MyFunction(double x, double y);
void SaveResults(int world_size, double total_elapsed);

struct timespec begin, end; // requires key –lrt  !
double elapsed;

double result = 0;

sem_t sem;

int main()
{
    /* --- initialize semaphore --- */
    sem_init(&sem, 0, 1);

    /* sequential approach */
    clock_t sqe_start, sqe_end;
    sqe_start = clock();

    double total_sum = 0;
    double estimate = 0;
    struct timespec time2;
    // clock_gettime(CLOCK_REALTIME, &time2);
    unsigned int seed_x = 122;
    unsigned int seed_y = 123;
    unsigned int seed_z = 124;
    double rand_x = 0;
    double rand_y = 0;
    double rand_z = 0;
    for (size_t i = 0; i < POINTS_NUM; i++)
    {
        rand_x = PI * (float)rand_r(&seed_x) / RAND_MAX;
        rand_y = (float)rand_r(&seed_y) / RAND_MAX;
        rand_z = 2 * (float)rand_r(&seed_z) / RAND_MAX;
        if (rand_y <= std::sin(rand_x) && rand_z <= rand_x * std::sin(rand_x))
        {
            total_sum += 1;
        }
        estimate = 2 * PI * total_sum / POINTS_NUM;
    }
    sqe_end = clock();
    double time_taken = double(sqe_end - sqe_start) / double(CLOCKS_PER_SEC);

    /* ------- parallel approach  ------- */
    clock_gettime(CLOCK_REALTIME, &begin);

    pthread_t pthr[NUM_THREADS];
    void *integral_answer;

    for (int i = 0; i < NUM_THREADS; i++)
    {
        int rc = pthread_create(&pthr[i], NULL, MonteCarloEstimate, nullptr);
        if (rc)
            printf("ERROR; return code from pthread_create() is %d \n", rc);
    }

    for (int i = 0; i < NUM_THREADS; i++)
    {
        int rc = pthread_join(pthr[i], &integral_answer);
        result += *(double *)integral_answer;
        free(integral_answer);
        if (rc)
            printf("ERROR; return code from pthread_join() is %d \n", rc);
    }
    result /= 2;
    std::cout << "integ res = " << result * 2 * PI / POINTS_NUM << std::endl;

    clock_gettime(CLOCK_REALTIME, &end);

    elapsed = end.tv_sec - begin.tv_sec; // time in seconds
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
    std::cout << "Elapsed time: " << elapsed << std::endl;
    double speedup = elapsed / time_taken;
    std::cout << "Speed up =  " << speedup << std::endl;

    /* --- kill semaphore --- */
    sem_destroy(&sem); // releasing of the semaphore
    SaveResults(NUM_THREADS, speedup);
}

void SaveResults(int num_threads, double speedup)
{
    std::ofstream output;
    std::string file_name = "P=" + std::to_string(NUM_THREADS) + ".txt";
    output.open(file_name);
    output << NUM_THREADS << "\n"
           << speedup;
}

void *MonteCarloEstimate(void *param)
{
    double V_0 = 0;
    double *rand_x = (double*)malloc(sizeof(double));
    double *rand_y = (double*)malloc(sizeof(double));
    double *rand_z = (double*)malloc(sizeof(double));

    double *estimate = (double *)malloc(sizeof(double));

    struct timespec time;
    unsigned int *seed_x = (unsigned int *)malloc(sizeof(unsigned int));
    unsigned int *seed_y = (unsigned int *)malloc(sizeof(unsigned int));
    unsigned int *seed_z = (unsigned int *)malloc(sizeof(unsigned int));

    clock_gettime(CLOCK_REALTIME, &time);
    *seed_x = time.tv_nsec;
    clock_gettime(CLOCK_REALTIME, &time);
    *seed_y = time.tv_nsec;
    clock_gettime(CLOCK_REALTIME, &time);
    *seed_z = time.tv_nsec;

    for (size_t i = 0; i < POINTS_NUM / NUM_THREADS; i++)
    {
        *rand_x = ((double)rand_r(seed_x) / RAND_MAX) * PI;
        *rand_y = (double)rand_r(seed_y) / RAND_MAX;
        *rand_z = ((double)rand_r(seed_z) / RAND_MAX) * 2;
        if (*rand_y <= std::sin(*rand_x) && *rand_z <= (*rand_x) * (*rand_y))
        {
            V_0++;
        }
    }
    *estimate = V_0;
    /* --- critical section --- */
    sem_wait(&sem);
    result += V_0;
    sem_post(&sem);

    /* ---- deaclocating memory ---- */
    free(rand_x);
    free(rand_y);
    free(rand_z);
    free(seed_x);
    free(seed_y);
    free(seed_z);

    return (void *)estimate;
}
