/*   Multiplication of square matrices by the use of  ALgo_Fox's algorithm 

 * Arguments:
 * n FILE1_name FILE1_name_input FILE2_name_input FILE3_name_output
 * n = oDR of the matrices [n*n matrix]
 FLOAT number input/output

*/

 /*-ln for sqrt link
 //mpicc -o <out_put_file_name> FOX_modif.c -lm       
 // mpirun -n 1,4,9,16 ./temp21 8 mdata1.txt mdata2.txt m_out.txt
*/




//TO compile and run
//mpicc -o <output_file_name> file_name.c -lm
//mpirun -n 1,4,9,16 ./temp <size_of_matrix> <input_matrix_file_name1> <input_matrix_file_name2> <output_matrix_file_name>
//mpd& is useful run when file cant be executed for kraken only




#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>


/*STRUCT GRID*/typedef struct {
/*STRUCT GRID*/    int       p;         /* Number of processor in GRID   */
/*STRUCT GRID*/    MPI_Comm  comm;      /* Grid communicator*/
/*STRUCT GRID*/    MPI_Comm  row_comm;  /* ROw communicator      */
/*STRUCT GRID*/    MPI_Comm  col_comm;  /* coloum communicator      */
/*STRUCT GRID*/    int       q;         /* q=sqrt(p) number of processor in a row of coum of grid */
/*STRUCT GRID*/    int       my_row_number;    /* My row number                */
/*STRUCT GRID*/    int       my_col_number;    /* My column number             */
/*STRUCT GRID*/    int       my_rank;   /* My rank in the grid comm     */
/*STRUCT GRID0*/} GRD_IFO_T;

//#define MAX 65000
/*STRUCT BLOCK*/#define MAX 131072 
/*STRUCT BLOCK*/typedef struct {
/*STRUCT BLOCK*/    int n_block;
/*STRUCT BLOCK*/    #define oDR(A) ((A)->n_block) /* MACRO FUNCTION DEFINED */
/*STRUCT BLOCK*/    float  entries[MAX];
/*STRUCT BLOCK*/    #define P_Entery(A,i,j) (*(((A)->entries) + ((A)->n_block)*(i) + (j)))  /* MACRO FUNCTION DEFINED */
/*STRUCT BLOCK*/} LOCAL_MTX_T;

/* Function Declarations */
LOCAL_MTX_T*  local_mtx_Allocation(int n_block);
void             Free_local_matrix(LOCAL_MTX_T** local_mtx_A);
void             reading_matrix_from_file(char* prompt, LOCAL_MTX_T* local_mtx_A,GRD_IFO_T* grid, int n);
void             Write_MTX_to_FILE(char* title, LOCAL_MTX_T* local_mtx_A, GRD_IFO_T* grid, int n);
void             Set_to_zero(LOCAL_MTX_T* local_mtx_A);
void             Local_matrix_multiplcation(LOCAL_MTX_T* local_mtx_A, LOCAL_MTX_T* local_mtx_B, LOCAL_MTX_T* local_mtx_C);
void             Blt_mtx_type(LOCAL_MTX_T* local_mtx_A);
MPI_Datatype     local_matrix_mpi_t;
LOCAL_MTX_T*  tmp_mtx;

/*********************************************************/

main(int argc, char* argv[]) {
/**/    int              p; /*total number of processor */
/**/    int              my_rank; 
/**/    GRD_IFO_T      grid; 
/**/    LOCAL_MTX_T*  local_mtx_A; 
/**/    LOCAL_MTX_T*  local_mtx_B;
/**/    LOCAL_MTX_T*  local_mtx_C;
/**/    int              n; /* matrix SIZE */
/**/    int              n_block; /*LATER INITALIZE BLOCK SIZE */ /*SIZE OF MATRIX / SQRT(NUMBER OF PROCESSOR) */
/**/    char input_file_name_1[255];
/**/    char input_file_name_2[255];
/**/    char output_file_name[255];
/**/    //double      wtime_overhead;
/**/    double      start, finish;
/**/    double      t_time;
/**/    int i;
/**/   // clock_t v_begin = clock();
/**/    
/**/    void Setup_grid(GRD_IFO_T*  grid);
/**/    void ALgo_Fox(int n, GRD_IFO_T* grid, LOCAL_MTX_T* local_mtx_A,LOCAL_MTX_T* local_mtx_B, LOCAL_MTX_T* local_mtx_C);
/**/
/**///start MPI
/**/    MPI_Init(&argc, &argv);
/**/    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
/**/
/**/    Setup_grid(&grid);
/**/
/**/    if (my_rank == 0) {
/**/      if (argc >4)
/**/      {
/**/   	    n=atoi(argv[1]);    /* atoi string to int*/
/**/        strcpy(input_file_name_1,argv[2]);
/**/        strcpy(input_file_name_2,argv[3]);
/**/        strcpy(output_file_name,argv[4]);
/**/      }
/**/      else
/**/      {
/**/        printf("Usage: %s n file1 file2 file3 ....\n",argv[0]);
/**/        printf("WRITE COMMAND AS mpirun -n (1,4,9,16) ./out n input_file_name_1 input_file_name_2 output_file_name \n");
/**/        exit(1);
/**/      }
/**/    }
/**/
/**/    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
/**/    MPI_Bcast(input_file_name_1,sizeof(input_file_name_1)+1,MPI_CHAR,0,MPI_COMM_WORLD);
/**/    MPI_Bcast(input_file_name_2,sizeof(input_file_name_2)+1,MPI_CHAR,0,MPI_COMM_WORLD);
/**/    MPI_Bcast(output_file_name,sizeof(output_file_name)+1,MPI_CHAR,0,MPI_COMM_WORLD);
/**/    
/**/    n_block = n/grid.q;   /*SIZE OF MATRIX / SQRT(NUMBER OF PROCESSOR) */
/**/
/**/    local_mtx_A = local_mtx_Allocation(n_block);
/**/    oDR(local_mtx_A) = n_block;
/**/    reading_matrix_from_file(input_file_name_1, local_mtx_A, &grid, n);
/**/    //Write_MTX_to_FILE("We read A =", local_mtx_A, &grid, n);
/**/    local_mtx_B = local_mtx_Allocation(n_block);
/**/    oDR(local_mtx_B) = n_block;
/**/    reading_matrix_from_file(input_file_name_2, local_mtx_B, &grid, n);
/**/    //Write_MTX_to_FILE("We read B =", local_mtx_B, &grid, n);
/**/    Blt_mtx_type(local_mtx_A);
/**/    tmp_mtx = local_mtx_Allocation(n_block);
/**/    local_mtx_C = local_mtx_Allocation(n_block);
/**/    oDR(local_mtx_C) = n_block;
/**/    start = MPI_Wtime();
/**/ 
/**/    /* MAIN CALCULATION ALgo_Fox calling */
/**/    ALgo_Fox(n, &grid, local_mtx_A, local_mtx_B, local_mtx_C);
/**/    /* MAIN CALCULATION ALgo_Fox calling over */
/**/
/**/    finish = MPI_Wtime();
/**/    t_time = finish - start;
/**/    printf("RANK:%d TIME:%f\n",my_rank ,t_time);
/**/    Write_MTX_to_FILE(output_file_name, local_mtx_C, &grid, n);
/**/    //Free memory
/**/    Free_local_matrix(&local_mtx_A);//Free memory
/**/    Free_local_matrix(&local_mtx_B);//Free memory
/**/    Free_local_matrix(&local_mtx_C);//Free memory
/**/
/**/    MPI_Finalize();
/**/    //END MPI
/**/    //clock_t v_end = clock();
/**/    //double time_spent = (double)(v_end - v_begin) / CLOCKS_PER_SEC;
/**/    //printf("TOTAL TIME TO EXECUTE:%f SECONDS\n",time_spent);

}  /* main */


/*****************FUUNCTION DEFINITIONS***************************/
void Setup_grid(GRD_IFO_T*  grid) {
    int old_rank;
    int dimensions[2]; /* 2d GRID*/ /* dimensions of Cartesian grid . */
    int wrap_around[2]; /* for rotation*/
    int coordinates[2];
    int free_coords[2];

    /* Set up Global Grid Information */
    MPI_Comm_size(MPI_COMM_WORLD, &(grid->p));
    MPI_Comm_rank(MPI_COMM_WORLD, &old_rank);

    /* We assume p is a perfect square 1 4 9 16 etc  */
    grid->q = (int) sqrt((double) grid->p);
    dimensions[0] = dimensions[1] = grid->q; /*specifying the number of processes in each dimension*/

    /*  second dimension circular shift */
    
    wrap_around[0] = wrap_around[1] = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, wrap_around, 1, &(grid->comm));
    MPI_Comm_rank(grid->comm, &(grid->my_rank));
    MPI_Cart_coords(grid->comm, grid->my_rank, 2, coordinates);
    grid->my_row_number = coordinates[0];
    grid->my_col_number = coordinates[1];

    /*  MPI_Cart_sub can be used to partition the communicator group into subgroups */
    /*  communication row */
    free_coords[0] = 0; 
    //(0,1)
    free_coords[1] = 1;
/* free cord for ith dimension of (0,1) is kept in the subgrid is kept or not*/

    MPI_Cart_sub(grid->comm, free_coords, &(grid->row_comm));
    /*  col communication */
    free_coords[0] = 1; 
    //(1,0)
    free_coords[1] = 0;
    /* free cord for ith dimension (1,0) is kept in the subgrid or not */
    MPI_Cart_sub(grid->comm, free_coords, &(grid->col_comm)); 
} /* Setup_grid */


/*********************************************************/
void ALgo_Fox( int n,GRD_IFO_T* grid ,LOCAL_MTX_T*  local_mtx_A,LOCAL_MTX_T*  local_mtx_B/**/,LOCAL_MTX_T*  local_mtx_C) {

    LOCAL_MTX_T*  temp_A; /*  sub-    *//* matrix store sopace of A used during */ 
                             /* the current stage       */
    int              stage;
    int              bcast_root;
    int              n_block;  /* n/sqrt(p)               */
    int              source; /*from*/
    /*to*/
    int              dest;
    MPI_Status       status;
/*set block size*/
    n_block = n/grid->q;
    Set_to_zero(local_mtx_C);

    /* circular shift of B addresses Calculate*/  
    source = (grid->my_row_number + 1) % grid->q;
    dest = (grid->my_row_number + grid->q - 1) % grid->q;

    /* tenp Set for stor for the b_cast block of A */
    temp_A = local_mtx_Allocation(n_block);

/*fox for*/    for (stage = 0; stage < grid->q; stage++) {
/*fox for*/        bcast_root = (grid->my_row_number + stage) % grid->q;
/*fox for*/        if (bcast_root == grid->my_col_number) {
/*fox for*/            MPI_Bcast(local_mtx_A, 1, local_matrix_mpi_t,bcast_root, grid->row_comm);
/*fox for*/            Local_matrix_multiplcation(local_mtx_A, local_mtx_B,local_mtx_C);
/*fox for*/        } else {
/*fox for*/            MPI_Bcast(temp_A, 1, local_matrix_mpi_t,bcast_root, grid->row_comm);
/*fox for*/            Local_matrix_multiplcation(temp_A, local_mtx_B,local_mtx_C);
/*fox for*/        }
/*fox for*/        MPI_Sendrecv_replace(local_mtx_B, 1, local_matrix_mpi_t,dest, 0, source, 0, grid->col_comm, &status);
/*fox for*/    } 
    
} 

/*********************************************************/

/*********************************************************/


/*********************************************************/
LOCAL_MTX_T* local_mtx_Allocation(int local_oDR) {

    LOCAL_MTX_T* temp;
    temp = (LOCAL_MTX_T*) malloc(sizeof(LOCAL_MTX_T));
    if (temp == NULL) printf("ERROR MATRIX ALLOCATION iiii0!!!!");
    return temp;

}  /* local_mtx_Allocation */


/*********************************************************/

/*********************************************************/


/*********************************************************/
/* Reading and distribute matrix:  
 *     for  row in global_matrix,
 *         for every grid column 
 *             read a block of n_block floats on process 0
 *             and send them to the appropriate process. and recv on processor
        if processor_0:
            send
        else:
            recv
 */
void reading_matrix_from_file(
         char            file_name[255]   /* in  */, 
         LOCAL_MTX_T*  local_mtx_A  /* out */,
         GRD_IFO_T*     grid     /* in  */,
         int              n        /* in  */) {

    int        mat_row, mat_col;
    int        grid_row, grid_col;
    int        dest;
    int        coords[2];
    float*     temp;
    MPI_Status status;
    FILE *fp; 
    fp=fopen(file_name,"r");
    if (grid->my_rank == 0) {
        temp = (float*) malloc(oDR(local_mtx_A)*sizeof(float));
        if (temp == NULL) printf("Error->allocaiton error in reading matrix function\n");
        printf("%s\n", file_name);
        fflush(stdout);
        for (mat_row = 0;  mat_row < n; mat_row++) {  /* for every matrix row*/
            grid_row = mat_row/oDR(local_mtx_A);            
            coords[0] = grid_row;
            for (grid_col = 0; grid_col < grid->q; grid_col++) { /* for every grid col */
                coords[1] = grid_col;
                MPI_Cart_rank(grid->comm, coords, &dest); 
                if (dest == 0) {
                    for (mat_col = 0; mat_col < oDR(local_mtx_A); mat_col++)
                        fscanf(fp,"%f", (local_mtx_A->entries)+mat_row*oDR(local_mtx_A)+mat_col);/*from fileinpt */
                } else {
                    for(mat_col = 0; mat_col < oDR(local_mtx_A); mat_col++)
                        fscanf(fp,"%f", temp + mat_col); /*from fileinpt */
                    MPI_Send(temp, oDR(local_mtx_A), MPI_FLOAT, dest, 0,grid->comm); /* distribute to other processors*/
                }
            }
        }
        fclose(fp);
        free(temp);
    } else {
        for (mat_row = 0; mat_row < oDR(local_mtx_A); mat_row++) 
            MPI_Recv(&P_Entery(local_mtx_A, mat_row, 0), oDR(local_mtx_A), MPI_FLOAT, 0, 0, grid->comm, &status); /*RECV*/
    }
                     
}  /* reading_matrix_from_file */


/*********************************************************/

/*********************************************************/
void Write_MTX_to_FILE(
         char            file_name[255]    /* in  */,  
         LOCAL_MTX_T*  local_mtx_A  /* out */,
         GRD_IFO_T*     grid     /* in  */,
         int              n        /* in  */) {
    int        mat_row, mat_col;
    int        grid_row, grid_col;
    int        source;
    int        coords[2];
    float*     temp;
    FILE *fp;
    MPI_Status status;
    fp=fopen(file_name,"w");//OPEN to write
    if (grid->my_rank == 0) {//by root
        temp = (float*) malloc(oDR(local_mtx_A)*sizeof(float));
        if(temp == NULL) printf("Error ->Error in Write_MTX_to_FILE FUNCTION iii0V!!! \n");
//        printf("%s\n", title);
        for (mat_row = 0;  mat_row < n; mat_row++) {
            grid_row = mat_row/oDR(local_mtx_A);
            coords[0] = grid_row;
            for (grid_col = 0; grid_col < grid->q; grid_col++) {
                coords[1] = grid_col;
                MPI_Cart_rank(grid->comm, coords, &source);
                if (source == 0) {
                    for(mat_col = 0; mat_col < oDR(local_mtx_A); mat_col++)
                        fprintf(fp,"%f ", P_Entery(local_mtx_A, mat_row, mat_col)); /*to file out*/
                } else {
                    MPI_Recv(temp, oDR(local_mtx_A), MPI_FLOAT, source, 0,grid->comm, &status);
                    for(mat_col = 0; mat_col < oDR(local_mtx_A); mat_col++)
                        fprintf(fp,"%f ", temp[mat_col]); /*to file out*/
                }
            }
            fprintf(fp,"\n");
        }
	fclose(fp);
        free(temp);
    } else {
        for (mat_row = 0; mat_row < oDR(local_mtx_A); mat_row++) 
            MPI_Send(&P_Entery(local_mtx_A, mat_row, 0), oDR(local_mtx_A), MPI_FLOAT, 0, 0, grid->comm);
    }
                     
}  /* Write_MTX_to_FILE */


/*********************************************************/

/*********************************************************/
void Blt_mtx_type(LOCAL_MTX_T*  local_mtx_A  /* in */) {
    MPI_Datatype  temp_mpi_t;
    int           block_lengths[2];
    MPI_Aint      displacements[2];
    MPI_Datatype  typelist[2];
    MPI_Aint      start_address;
    MPI_Aint      address;

    /*The simplest datatype constructor is MPI_Type_contiguous,
    which allows replication of a datatype
    into contiguous locations.[https://www.open-mpi.org/doc/v4.0/man3/MPI_Type_contiguous.3.php]*/

    MPI_Type_contiguous(oDR(local_mtx_A)*oDR(local_mtx_A), MPI_FLOAT, &temp_mpi_t);

    block_lengths[0] = block_lengths[1] = 1;
   
    typelist[0] = MPI_INT;
    typelist[1] = temp_mpi_t;

    MPI_Address(local_mtx_A, &start_address);
    MPI_Address(&(local_mtx_A->n_block), &address);
    displacements[0] = address - start_address;
    
    MPI_Address(local_mtx_A->entries, &address);
    displacements[1] = address - start_address;

    MPI_Type_struct(2, block_lengths, displacements,typelist, &local_matrix_mpi_t);
    MPI_Type_commit(&local_matrix_mpi_t); 
}  /* Blt_mtx_type */


/*********************************************************/

/*********************************************************/
void Local_matrix_multiplcation(
         LOCAL_MTX_T*  local_mtx_A  /* in  */,
         LOCAL_MTX_T*  local_mtx_B  /* in  */, 
         LOCAL_MTX_T*  local_mtx_C  /* out */) {
    int i, j, k;
/*
for i < SIZE:
    fo j < SIZE:
        for k SIE:
            c[i][j]=c[i][j]+a[i][k]*b[k][j]
*/
    for (i = 0; i < oDR(local_mtx_A); i++)
        for (j = 0; j < oDR(local_mtx_A); j++)
            for (k = 0; k < oDR(local_mtx_B); k++)
                P_Entery(local_mtx_C,i,j) = P_Entery(local_mtx_C,i,j) 
                    + P_Entery(local_mtx_A,i,k)*P_Entery(local_mtx_B,k,j);

}  /* Local_matrix_multiplcation */




/*********************************************************/

/*********************************************************/
void Set_to_zero(LOCAL_MTX_T*  local_mtx_A) {

    int p, q; /*row,col*/

    for (p = 0; p < oDR(local_mtx_A); p++) /*r,c*/
        for (q = 0; q < oDR(local_mtx_A); q++)  /*r,c*/
            P_Entery(local_mtx_A,p,q) = 0.0;  

}  /* Set_to_zero */


/*********************************************************/

/*********************************************************************/
void Free_local_matrix(LOCAL_MTX_T** local_mtx_A_ptr  /* in/out */) {
    free(*local_mtx_A_ptr);
}  /* Free_local_matrix */

