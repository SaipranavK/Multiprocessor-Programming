parts = n/p;
npart = parts;
if processor = master
{
	for slave processors {MPI_Send(part_matrix);} /* Master send part of matrix */ 	
}
else if processors = slave {MPI_Recv(part_matrix);} /* Recieve part of matrix from master */
while(!converged)
{
	if processors = slave{
		MPI_Send(matrix[1][0]); /* Share top with previous node*/ 
		MPI_Recv(matrix[0][0]);
	}
	if not last processor{
		MPI_Recv(matrix[parts+1][0]); /* Share top with previous node*/
		MPI_Send([parts][0]);
	}
	if iteration is odd{
		if matrix[cell][cell] == odd{
			A[cell][cell] = (1 - w) * A[cell][cell] + w *(A[cell-1][cell] + A[cell+1][cell] + A[cell][cell - 1] + A[cell][cell + 1])/4; /* SOR */
		}
		MPI_Reduce(MaxAcrossAll); /*find max across all nodes*/

		if((MaxAcrossAll - prevoddMax) <= difflimit){converged;}
		prevoddMax = MaxAcrossAll;
		iterations++;
		if iterations > 1000000 stop;
	}
	else if iterations is even{
		//compute same as odd but for even cells
		.
		if((MaxAcrossAll - prevevenMax) <= difflimit){converged;}
		prevevenMax = MaxAcrossAll;
		.
	}
}
if processors = slave{MPI_Send(matrix);} /* Send back matrix on convergence */
else if processors = master{MPI_Recv(matrix);} /* Recieve from slave */
