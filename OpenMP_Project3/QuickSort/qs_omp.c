/***************************************************************************
 *
 * Sequential version of Quick sort
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include<omp.h>


#define KILO (1024)
#define MEGA (1024*1024)
#define MAX_ITEMS (5)

#define swap(v, a, b) {unsigned tmp; tmp=v[a]; v[a]=v[b]; v[b]=tmp;}

static int *v;

static void
print_array(void)
{
    int i;

    for (i = 0; i < MAX_ITEMS; i++)
        printf("%d ", v[i]);
    printf("\n");
}

static void
init_array(void)
{
    int i;

    v = (int *) malloc(MAX_ITEMS*sizeof(int));
    for (i = 0; i < MAX_ITEMS; i++)
        v[i] = rand();
}

static unsigned
partition(int *v, unsigned low, unsigned high, unsigned pivot_index)
{
    /* move pivot to the bottom of the vector */
    if (pivot_index != low)
        swap(v, low, pivot_index);

    pivot_index = low;
    low++;

    /* invariant:
     * v[i] for i less than low are less than or equal to pivot
     * v[i] for i greater than high are greater than pivot
     */

    /* move elements into place */
    while (low <= high) {
        if (v[low] <= v[pivot_index])
            low++;
        else if (v[high] > v[pivot_index])
            high--;
        else
            swap(v, low, high);
    }

    /* put pivot back between two groups */
    if (high != pivot_index)
        swap(v, pivot_index, high);
    return high;
}

static void
quick_sort(int *v, unsigned low, unsigned high)
{
    unsigned pivot_index;
    
    /* no need to sort a vector of zero or one element */
    if (low >= high)
        return;

    /* select the pivot value */
    pivot_index = (low+high)/2;

    /* partition the vector */
    pivot_index = partition(v, low, high, pivot_index);

    /* sort the two sub arrays */
    if (low < pivot_index)
    {
        #pragma omp task firstprivate(v,low,pivot_index) //Quicksort subarray parallely
        {
            quick_sort(v, low, pivot_index-1);
        }
    }

    if (pivot_index < high) 
    {
        #pragma omp task firstprivate(v,low,pivot_index) //Quicksort subarray parallely
        {
            quick_sort(v, pivot_index+1, high);
        }
        
    }
}

int
main(int argc, char **argv)
{
    omp_set_num_threads(8);
    double start_time, run_time;
    
    printf("\nArray Initialised:");
    init_array();
    
    print_array();

    start_time = omp_get_wtime();
    #pragma omp parallel
    {

        int id = omp_get_thread_num();
        int nthrds = omp_get_num_threads();
        #pragma omp single nowait

        quick_sort(v, 0, MAX_ITEMS-1);

    }
    run_time = omp_get_wtime() - start_time;

    printf("\nArray sorted:");
    print_array();
  
    printf("\nTime Elapsed: %lf seconds\n ",run_time);
    return 1;
}
