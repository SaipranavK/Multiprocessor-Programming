QuickSort(A,low,high)
{
	if(low<high)
	{
		pivot_index = partition(A,high,end);
		QuickSort(A,low,pivot_index-1); //parallelize
		QuickSort(A,pivot_index+1,end); //parallelize
	}
}