#include<stdio.h>
#include <stdlib.h>
#include <time.h>

#define SIZE 40000

void mergeSort(int arr[],int,int);     // Prototype functions for mergeSort and merge functions.
void merge(int arr[],int,int,int);



int main(){
    
    clock_t start_time;  // clock_t is a data type for measuring processor time in clock ticks.
    clock_t end_time;
    double total_time;

    
    
    // Initialize the array to store random numbers
    int arr[SIZE];

    // Seed the random number generator
    srand(time(0));

    // Generate random numbers and store them in the array
    for (int i = 0; i < SIZE; i++) {
        arr[i] = rand();
    }
    printf("\n");

    start_time=clock();

    mergeSort(arr,0,SIZE-1); // Implementing Merge Sort.

    end_time=clock();   // Ending the clock here as merge sort function is implemented.

    total_time=((double)(end_time - start_time)) / CLOCKS_PER_SEC;  //Measuring the time taken.
    printf("Total elapsed time is %f seconds", total_time);
    
    return 0;
}


void merge(int arr[], int beg, int mid, int end) // Function using the conquer approach to merge and sort the divided sub-arrays.
{ 
    int i=beg, j=mid+1, index=beg, k; 
    int temp[SIZE];
    while ((i<=mid) &&(j<=end))   // Copying the elements into a temporary array after sorting them.
    {
          if (arr[i]<arr[j])
          {
            temp[index]=arr[i];
            i++;
          }
          else 
          {
            temp[index]=arr[j];
            j++;
          }
          index++;
    }
    if (i>mid)               // Copying the remaining elements of the right subarray if there exists any.
    {
        while (j<=end)
        {
            temp[index]=arr[j];
            j++;
            index++;
        }
    }
    else
    {
        while (i<=mid)      // Copying the remaining elements of the left subarray if there exists any.
        {
            temp[index]=arr[i];
            i++;
            index++;
        }
    }
   
    for (k=beg;k<index;k++)        // Finally copying the sorted elements of the temporary back into the original array.
    {
        arr[k]=temp[k];
    }
}


void mergeSort(int arr[], int beg, int end) // Function using the divide approach to divide the array into subarrays recursively until all subarrays have 1 element.
{ 
    if (beg < end) { 
       
        int mid = beg + (end - beg) / 2; 
  
        // Sort first and second halves.
        mergeSort(arr, beg, mid); 
        mergeSort(arr, mid + 1, end); 

        // Merge the sorted halves.
        merge(arr, beg, mid, end); 
    } 
} 