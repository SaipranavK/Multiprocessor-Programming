/*
*****************************
*Denim Deshmukh      *
*Koyyada Sai Pranav  *
*****************************
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <pthread.h>

#define _REENTRANT
#define BILLION  1000000000.0
#define MAXN 2050

int N;

volatile float A[MAXN][MAXN] ,B[MAXN],X[MAXN];

int NumThreads;

pthread_t Threads[_POSIX_THREAD_THREADS_MAX];
pthread_mutex_t Mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t Global_CountLock = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t NextIter = PTHREAD_COND_INITIALIZER;

int Norm,G_CurrentRow,Global_Count;


int minimum(a,b)
{
  if(a<b)
  {
    return a;
  }
  else
    return b;
}


void init_input()
{
	int i,j;
	printf("\nIntitializing matrix.......\n");
	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
		if (i == j) /* diagonal dominance */
		    A[i][j] = 5.0;
		else
		    A[i][j] = 2.0;
	    }
	}
	for (i = 0; i < N; i++) {
	B[i] = 2.0;
	X[i] = 1.0;
    }

}

void print_matrix_A_B()
{
  int i,j;
  printf("\n Matrix A:\n");
  printf("[");
     for (i = 0; i < N; i++)
     {
      /* code */
      for (j = 0; j < N; j++)
      {
        /* code */
        printf("%6.3f ", A[i][j]);
        if(j<N-1)
        {
          printf(", ");
        }
        else
        {
          printf("\n\t");
        }
                        // if j < N-1 print , else matrix over end next line 
      }
      }
     printf("]\n");
     printf("\nMatrix B: \n [");
     for (j = 0; j < N; j++)
     {
      /* code */
      printf("%6.3f ", B[j]);
      if(j<N-1)
        {
          printf(", ");
        }
        else
        {
          printf("]\n");
        }
                    // if j < N-1 print , else matrix over end next line 
     }
     
}



 void print_X()
{
  int i;
  printf("\nX = [");
     for (i = 0; i < N; i++)
     {
      /* code */
      printf("%6.3f  ", X[i]);
      if (i<N-1)
        {
          /* code */
          printf(", ");
        }
      else
      {
        printf("]\n");
      } 
     }
     
}

int get_block(int *myRow, int myNorm)
{
 int chunkSize;
 pthread_mutex_lock(&Mutex);
 *myRow = G_CurrentRow;
 chunkSize = (*myRow < N) ? (N-myNorm-1)/(2*NumThreads)+1 : 0;
 G_CurrentRow += chunkSize; 
 pthread_mutex_unlock(&Mutex);
 return chunkSize;
 }


void sync(int *myNorm)
{
 
 pthread_mutex_lock(&Global_CountLock);  
 if (Global_Count == 0)
    {
     Norm++;
     Global_Count = NumThreads-1;
     G_CurrentRow = Norm+1;
     pthread_cond_broadcast(&NextIter);
     }
 else
    {
     Global_Count--;
     pthread_cond_wait(&NextIter, &Global_CountLock);
     }

 *myNorm = Norm;

 pthread_mutex_unlock(&Global_CountLock);
 }

 void *gauss(void *required_unused)
{
 int myRow = 0, row, col;

 int myNorm = 0;

 float ratio;
 int chunkSize;

 while (myNorm < N-1)
 {
  while (chunkSize = get_block(&myRow, myNorm))
  {  

   for (row = myRow; row < (minimum(N, myRow+chunkSize)); row++)
   {
    ratio = A[row][myNorm]/A[myNorm][myNorm];
    for (col = myNorm; col < N; col++)
    A[row][col] -= A[myNorm][col]*ratio;
    B[row] -= B[myNorm]*ratio;
    }
   }

  sync(&myNorm);
  }
 }


void create_threads()
{
  int i;
  for (i = 0; i < NumThreads; i++)
  {
    pthread_create(&Threads[i], NULL, gauss, NULL);
  }
}


void wait_for_threads()
{
  int i;
  for (i = 0; i < NumThreads; i++)
  {
    pthread_join(Threads[i], NULL);
  }
}

int main()
{
  printf("\nEnter (Size of matrix) Limit 2050 N:\n");
  scanf("%d",&N);
  printf("\nEnter (Number of Threads) Limit 60 ThreadNum:\n");
  scanf("%d",&NumThreads);
	int row, col;
  struct timespec start, end;
	init_input();
  G_CurrentRow = Norm+1;
 	Global_Count = NumThreads-1;
 	clock_gettime(CLOCK_REALTIME, &start);
  create_threads();

 	wait_for_threads();
  clock_gettime(CLOCK_REALTIME, &end);
  print_matrix_A_B();

//Back substitution extra
 for (row = N-1; row >= 0; row--)
 {
  X[row] = B[row];
  for (col = N-1; col > row; col--)
  X[row] -= A[row][col]*X[col];
  X[row] /= A[row][row];
  }

 print_X();
 double t_t = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / BILLION;
 printf("SIZE:%d THREADS:%d\n",N,NumThreads );
  printf("Time to compute  in %f seconds:\n", t_t);
 
 return 0;

}

