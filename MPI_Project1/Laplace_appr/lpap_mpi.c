/***
Project 1 - Task 2
Parallel implementation of Laplace approximation using 
Successive-Over-Relaxation(SOR) numerical method
author : Koyyada Sai Pranav 
***/


/* Imports */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<mpi.h>


/* Default Declarations */
#define MAX_SIZE 4096   /* Maximum allowed dimensionality for matrix init */	
#define MSGINIT 0		
#define MASTER 1		/* Messaage Tags */
#define SLAVE 2
#define ODD 1			/* Red-black */
#define EVEN 0


/* CLI args */
int mtype;     		/* Message Type */
int N;				/* Matrix dimensionality */
int maxnum;			/* Maximum allowed dimensionality for matrix init */
char *Init;			/* Matrix init type */
double difflimit;	/* Stopping Condition */	
double w;			/* Relaxation Factor (Omega) */
int PRINT;			/* Debug or output matrix */

/* Matrix initialisation with +2 for boundary tuples and attrs */
static double A[MAX_SIZE+2][MAX_SIZE+2];

/* to capture status of message passing */
MPI_Status status;

/* forward declarations */
int seqwork();			/* Seqential implementation if only one processor is available */					 
int work(int,int);		/* Parallel implementation with multiple processors */			 		
void Init_Matrix();		/* Initialise matrix */
void Print_Matrix();	/* Print matrix */
void Init_Default();	/* initialise defaults */
int Read_Options(int, char **); /* Take CLI params from user */


/* Program Start */
int main(int argc, char **argv)
{
	int iter, processor_rank, np;  /* iter - iterations, np - number of processors */

	/* MPI */
	MPI_Init(&argc, &argv);	
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &processor_rank);

	printf("\nProcessor = %d/%d\n", processor_rank,np-1);

	double t_start = MPI_Wtime();  /* Capture start time */
	
	Init_Default();					/* Defaults initialised */
	Read_Options(argc,argv);		/* Overwrite Defaults with CLI params from user if any */

	MPI_Barrier(MPI_COMM_WORLD);    //Barrier for Synchronization
	if(np == 1)				/* if only one processor is available trigger seqential work */
	{
		if(processor_rank == 0)		/* Master Node */
		{
			Init_Matrix(processor_rank); 	/* Matrix initialised */
			iter = seqwork(); 				/* Trigger Sequential work */

			if(PRINT == 1)
			{
				Print_Matrix();
			}

			printf("Number of iterations = %d\n", iter );
		}
	}

	else if(np>1)		/* if more than one processor is available trigger parallel work*/
	{
		if(processor_rank == 0)   /* Master Node */
		{
			Init_Matrix(processor_rank);
			iter = work(processor_rank, np);

			if(PRINT == 1)
			{
				Print_Matrix();
			}

			printf("Number of iterations = %d\n", iter );
		}
		else
		{
			work(processor_rank, np);
		}
	}

	double t_end = MPI_Wtime(); /* Capture completion time */
	MPI_Finalize();

	printf("\nTime Elapsed in Processor %d/%d = %f\n", processor_rank, np-1, t_end - t_start);

	return 1;

}

/* Sequential Implementation */
int seqwork()
{
	int m,n;		/* Looping vars */
	int iteration = 0; /* Count iterations for covergence */
	int finished = 0;	/* Convergence check var */
	int turn = EVEN;	/* red-black approach */

	double maxi, sum, prevmax_even, prevmax_odd; /* Convergence check vars */


	prevmax_odd = 0.0;
	prevmax_even = 0.0;

	while(!finished)
	{
		if(turn == EVEN)			/* Check iteration */
		{
			turn = ODD;			/* Set next iteration */

			for(m = 1; m < N+1; m++) 
			{
				for(n = 1; n < N+1; n++)
				{
					if(((m+n)%2) == 0)	/* Check cell for red - even */	
					{
						A[m][n] = (1 - w) * A[m][n] + w *(A[m-1][n] + A[m+1][n] + A[m][n - 1] + A[m][n + 1])/4;  /* SOR */
					}
				}
			}

			maxi = -999999.0;
			for(m = 1; m < N + 1; m++)
			{
				sum = 0.0;
				for(n = 1; n < N+1; n++)
				{
					sum += A[m][n];
				}
				if(sum > maxi)
				{
					maxi = sum;
				}
			}	
			if(fabs(maxi - prevmax_even) <= difflimit) /* Stop if converged */
				finished = 1;
			
			prevmax_even = maxi;	/* Initialise as current maxi */	
		}

		else if(turn == ODD)			/* Check iteration */
		{
			turn = EVEN;				/* Set next iteration */

			for(m = 1; m < N+1; m++)
			{
				for(n = 1; n < N+1; n++)
				{
					if(((m+n)%2) == 1)		/* Check cell for black - odd */
					{
						A[m][n] = (1 - w) * A[m][n] + w *(A[m-1][n] + A[m+1][n] + A[m][n - 1] + A[m][n + 1])/4;	/* SOR */
					}
				}

			}

			maxi = -999999.0;
			for(m = 1; m < N + 1; m++)
			{
				sum = 0.0;
				for(n = 1; n < N+1; n++)
				{
					sum += A[m][n];
				}
				if(sum > maxi)
				{
					maxi = sum;
				}
			}

			if(fabs(maxi - prevmax_odd) <= difflimit)	/* Stop if converged */
				finished = 1;

			prevmax_odd = maxi;			/* Initialise as current maxi */
		}

		iteration++;		/* Increment iteration */

		if (iteration > 100000) 	/* If iteration count is beyond 1000000 */
		{   
			/* exit if we don't converge fast enough */    
			printf("Max number of iterations reached! Exit!\n");    
			finished = 1;
		} 
	}

	return iteration;
}

/* Parallel implementation */
int work(int rank, int p)
{
		MPI_Barrier(MPI_COMM_WORLD); //Barrier for Synchronization
		
		int cols = N + 2;       			/* Total number of columns to be processed by individual processor */
		int parts = N/p;					/* Possible partition of rows */
		int nparts = parts;					/* Var to partition matrix for slave nodes */ 
		int rows = parts + 2;				/* Total number of rows to be processed by individual processor */
		int iteration = 0;					/* Count iterations for covergence   */
		int finished = 0;					/* Convergence check vars */
		int turn = EVEN;					/* red-black approach */

		/* Looping variables */				
		int m, n, i, j, x, y;				
		int dest, src;
		
		/* Convergence check vars */						
		double maxi;
		double sum;
		double prevmax_even, prevmax_odd;
		prevmax_even = 0.0;
		prevmax_odd = 0.0;
		
		if(rank == 0)		/* Master Node */
		{
			for(dest = 1; dest < p; dest++)
			{
				for(i = 0; i < rows; i++)
				{					
					MPI_Send(&A[nparts + i][0], cols, MPI_DOUBLE, dest, MSGINIT, MPI_COMM_WORLD); /* Send parts to slave nodes */
				}
				nparts += parts;
			}				
		}
		
		if(rank != 0)		/* Slave Node */
		{
			for(j = 0; j<rows; j++)
			{
				MPI_Recv(&A[j][0], cols, MPI_DOUBLE, 0, MSGINIT, MPI_COMM_WORLD, &status); /* Receive respective part from Master */
			}
		}
		
		
		while (!finished) 
		{
			if (rank != 0)  /* Each node sends the top row shared with previous nodes except master node */
			{	
				mtype = SLAVE;
				MPI_Send(&A[1][0], cols, MPI_DOUBLE, rank-1 ,mtype, MPI_COMM_WORLD);
				mtype = MASTER;
				MPI_Recv(&A[0][0], cols, MPI_DOUBLE, rank-1, mtype, MPI_COMM_WORLD, &status);
			}
	
			if (rank != p-1) 
			{				/* Each node sends the bottom row shared with next nodes */
				mtype = SLAVE;
				MPI_Recv(&A[parts+1][0], cols , MPI_DOUBLE, rank+1, mtype, MPI_COMM_WORLD, &status);
				mtype = MASTER;
				MPI_Send(&A[parts][0], cols, MPI_DOUBLE, rank+1, mtype, MPI_COMM_WORLD);
				
			}
	
			if(turn == EVEN)		/* Check Iteration */
			{	
				turn = ODD;			/* Set next iteration */

				for(m = 1; m < parts+1; m++)
				{
					for(n = 1; n < N+1; n++)
					{
						if(((m+n) % 2) == 0)		/* Check cell for red - even */
						{
							A[m][n] = (1 - w) * A[m][n] + w *(A[m-1][n] + A[m+1][n] + A[m][n - 1] + A[m][n + 1])/4; 
						}
					}
				}
				
				/* Sum in each individual node */
				maxi = -999999.0;
				for (m = 1; m < parts + 1; m++) {
					sum = 0.0;
					for (n = 1; n < N+1; n++) {
						sum += A[m][n];
					}
					if(sum > maxi)
						maxi = sum;
				}
			
				double maxsum1; /* Store max across all nodes */

				MPI_Allreduce(&maxi, &maxsum1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); /* Maximum sum across all the nodes */
				
				if(fabs(maxsum1 - prevmax_even) <= difflimit)	/* Stop if converged */
					finished = 1;
				
				prevmax_even = maxsum1;	/* Initialise to current max across all nodes */
				
			}

			else if (turn == ODD)		/* Check iteration */
			{
				turn = EVEN;			/* Set next iteration */ 			

				for(m = 1; m < parts+1; m++)
				{
					for(n = 1; n < N+1; n++)
					{
						if(((m+n) % 2) == 1)		/* Check cell for black - odd */
						{
							A[m][n] = (1 - w) * A[m][n] + w *(A[m-1][n] + A[m+1][n] + A[m][n - 1] + A[m][n + 1])/4; /* SOR */ 
						}
					}
				}
				
				/* Sum in each individual node */
				maxi = -999999.0;
				for (m = 1; m < parts + 1; m++) {
					sum = 0.0;
					for (n = 1; n < N+1; n++) {
						sum += A[m][n];
					}
					if(sum > maxi)
						maxi = sum;
				}
			
				double maxsum2;	/* Store max across all nodes */

				MPI_Allreduce(&maxi, &maxsum2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); /* Maximum sum across all the nodes */
				
				if(fabs(maxsum2 - prevmax_odd) <= difflimit) /* Stop if converged */
					finished = 1;
				
				prevmax_odd = maxsum2;	/* Initialise to current max across all nodes */
				
			}
			
			iteration++;	/* Increment iterations  */
					
			if (iteration > 100000) /* If iteration count is beyond 1000000 */
			{
				/* exit if we don't converge fast enough */
				printf("Max number of iterations reached! Exit!\n");
				finished = 1;
			}

		}
		
		/* Sending results from Slave to the Master node */
		if (rank != 0)
		{
			mtype = SLAVE;
			for(x = 1; x < parts+1; x++)
			{
				MPI_Send(&A[x][0], N+2, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD);	
			}
		}
		
		
		int recvparts = parts+1;	/* parts to map the position of returning result rows in matrix A in Master Node */
		
		/* Receiving at Master node */
		if(rank == 0)
		{
			mtype = SLAVE;
			for(src = 1; src < p; src++)
			{
				for(y = 0; y<parts; y++)
				{
					MPI_Recv(&A[recvparts][0], N+2, MPI_DOUBLE, src, mtype, MPI_COMM_WORLD, &status);
				    recvparts++;
				}
			}	
		}
		
	return iteration;
}
 
void
Init_Matrix(int rank)
{
    int i, j, dmmy;
	
	printf("\nsize      = %dx%d ",N,N);
	printf("\nmaxnum    = %d \n",maxnum);
	printf("difflimit = %.7lf \n",difflimit);
	printf("Init	  = %s \n",Init);
	printf("w	  = %f \n\n",w);
	printf("Initializing matrix...");
	
    /* Initialize all grid elements, including the boundary */
    for (i = 0; i < N+2; i++) {
	for (j = 0; j < N+2; j++) {
	    A[i][j] = 0.0;
	}
    }
    if (strcmp(Init,"count") == 0) {
	for (i = 1; i < N+1; i++){
	    for (j = 1; j < N+1; j++) {
			A[i][j] = (double)i/2;
	    }
	}
    }
    if (strcmp(Init,"rand") == 0) {
	for (i = 1; i < N+1; i++){
	    for (j = 1; j < N+1; j++) {
			A[i][j] = (rand() % maxnum) + 1.0;
	    }
	}
    }
    if (strcmp(Init,"fast") == 0) {
	for (i = 1; i < N+1; i++){
	    dmmy++;
	    for (j = 1; j < N+1; j++) {
		dmmy++;
		if ((dmmy%2) == 0)
		    A[i][j] = 1.0;
		else
		    A[i][j] = 5.0;
	    }
	}
    }

    /* Set the border to the same values as the outermost rows/columns */
    /* fix the corners */
    A[0][0] = A[1][1];
    A[0][N+1] = A[1][N];
    A[N+1][0] = A[N][1];
    A[N+1][N+1] = A[N][N];
    /* fix the top and bottom rows */
    for (i = 1; i < N+1; i++) {
	A[0][i] = A[1][i];
	A[N+1][i] = A[N][i];
    }
    /* fix the left and right columns */
    for (i = 1; i < N+1; i++) {
	A[i][0] = A[i][1];
	A[i][N+1] = A[i][N];
    }

    printf("done in node : %d", rank);
    if (PRINT == 1)
		Print_Matrix();
}

void
Print_Matrix()
{
    int i, j;
 
    for (i=0; i<N+2 ;i++){
	for (j=0; j<N+2 ;j++){
	    printf(" %f",A[i][j]);
	}
	printf("\n");
    }
    printf("\n");
}


void 
Init_Default()
{
    N = 2048;
    difflimit = 0.00001*N;
    Init = "rand";
    maxnum = 15.0;
    w = 0.5;
    PRINT = 0;
}

int
Read_Options(int argc, char **argv)
{
    char    *prog;
 
    prog = *argv;
    while (++argv, --argc > 0)
	if (**argv == '-')
	    switch ( *++*argv ) {
	    case 'n':
		--argc;
		N = atoi(*++argv);
		difflimit = 0.00001*N;
		break;
	    case 'h':
		printf("\nHELP: try sor -u \n\n");
		exit(0);
		break;
	    case 'u':
		printf("\nUsage: sor [-n problemsize]\n");
		printf("           [-d difflimit] 0.1-0.000001 \n");
		printf("           [-D] show default values \n");
		printf("           [-h] help \n");
		printf("           [-I init_type] fast/rand/count \n");
		printf("           [-m maxnum] max random no \n");
		printf("           [-P print_switch] 0/1 \n");
		printf("           [-w relaxation_factor] 1.0-0.1 \n\n");
		exit(0);
		break;
	    case 'D':
		printf("\nDefault:  n         = %d ", N);
		printf("\n          difflimit = %f ", difflimit);
		printf("\n          Init      = rand" );
		printf("\n          maxnum    = 5 ");
		printf("\n          w         = 0.5 \n");
		printf("\n          P         = 0 \n\n");
		exit(0);
		break;
	    case 'I':
		--argc;
		Init = *++argv;
		break;
	    case 'm':
		--argc;
		maxnum = atoi(*++argv);
		break;
	    case 'd':
		--argc;
		difflimit = atof(*++argv);
		break;
	    case 'w':
		--argc;
		w = atof(*++argv);
		break;
	    case 'P':
		--argc;
		PRINT = atoi(*++argv);
		break;
	    default:
		printf("%s: ignored option: -%s\n", prog, *argv);
		printf("HELP: try %s -u \n\n", prog);
		break;
	    } 
}
