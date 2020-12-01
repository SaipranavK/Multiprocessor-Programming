 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <unistd.h>


#define MAX_SIZE 4096
#define BILLION  1000000000.0

typedef double matrix[MAX_SIZE][MAX_SIZE];

int	N;		            // matrix size	
matrix	A;		        // matrix A
double	b[MAX_SIZE];	// vector b 
double	y[MAX_SIZE];	// vector y 
int chunksize;          // chunk size 
int number_of_threads;

void work(void);
void Init_Matrix(void);
void Print_Matrix(void);



int main()
{
    int i;
    struct timespec start, end;
    printf("ENRER THE SIZE:\n");
    scanf("%d",&N);
    printf("Enter number of threads:\n");
    scanf("%d",&number_of_threads);
	
    Init_Matrix();	
    clock_gettime(CLOCK_REALTIME, &start);
	work();
	clock_gettime(CLOCK_REALTIME, &end);
	Print_Matrix();

	double time_spent = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / BILLION;
	printf("SIZE:%d\n",N);
	printf("THREAS:%d\n",number_of_threads);
	printf("Time elpased is %f seconds\n", time_spent);

}



void work(void)   /* Gaussian Elimination algorithm */
{
    int i, j, k;
	
    for (k = 0; k < N; k++) {   
		/* Creates a new group of threads */
		#pragma omp parallel for num_threads(number_of_threads) schedule(dynamic, chunksize)            
		for (j = k+1; j < N; j++)
			{                                 
			A[k][j] = A[k][j] / A[k][k];
			}                          

		y[k] = b[k] / A[k][k];	
		A[k][k] = 1.0;
		
		
		#pragma omp parallel for num_threads(number_of_threads) schedule(dynamic, chunksize) collapse(2)              
		for (i = k+1; i < N; i++) {		                           								       
			for (j = k+1; j < N; j++)		                  //colapse for 2 for loop
				A[i][j] = A[i][j] - A[i][k]*A[k][j];              // Elimination step 
		}

		#pragma omp parallel for num_threads(number_of_threads) schedule(dynamic, chunksize)
	-	for (i = k+1; i < N; i++) {		
			b[i] = b[i] - A[i][k]*y[k];			
			A[i][k] = 0.0;
		}
    }
}


void
Init_Matrix()
{
    int i, j;
 
    printf("\nsize      = %dx%d ", N, N);
    printf("Initializing matrix...");
 
	chunksize = N/number_of_threads;   // chunk size of each chunk that is allocated to a thread by omp
	
   	
    
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++) {
			if (i == j)                                      // diagonal dominance 
				A[i][j] = 5.0;
			else
				A[i][j] = 2.0;
		}
	}
    

    for (i = 0; i < N; i++) {                                    // Initialize vectors b and y
		b[i] = 2.0;
		y[i] = 1.0;
    }

    printf("init done \n\n");
	Print_Matrix();
}

void Print_Matrix()
{
    int i, j;
 
    printf("Matrix A:\n");
    for (i = 0; i < N; i++) {
    	printf("[");
		for (j = 0; j < N; j++)
	    	printf(" %5.3f,", A[i][j]);
		printf("]\n");
    }
    printf("Vector b:\n[");
    for (j = 0; j < N; j++)
		printf(" %5.3f,", b[j]);
    printf("]\n");
    printf("Vector y:\n[");
    for (j = 0; j < N; j++)
		printf(" %5.3f,", y[j]);
    printf("]\n");
    printf("\n\n");
}



