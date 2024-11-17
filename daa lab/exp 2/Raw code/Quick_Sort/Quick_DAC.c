#include<stdio.h>
#include <stdlib.h>
#include <time.h>

#define SIZE 60000


void swap(int* , int* );             // Prototype functions for swap, parition and quick sort functions.
int partition(int arr[], int , int );
void quickSort(int arr[], int , int );


int main(){


    clock_t start_time;  // clock_t is a data type for measuring processor time in clock ticks.
    clock_t end_time;
    double total_time;

    // Starting the clock to measure time.
    
    // Initialize the array to store random numbers
    int arr[SIZE];

    // Seed the random number generator
    srand(time(0));

    // Generate random numbers and store them in the array
    for (int i = 0; i < SIZE; i++) {
        arr[i] = rand();
    }
    printf("\n");

    start_time=clock();  // Starting the clock to measure time.

    quickSort(arr,0,SIZE-1); // Implementing Quick Sort.

    end_time=clock();   // Ending the clock here as Quick sort function is implemented.

    total_time=((double)(end_time - start_time)) / CLOCKS_PER_SEC;  //Measuring the time taken.
    printf("Total elapsed time is %f seconds", total_time);
    
    return 0;
}


// Function to swap two elements
void swap(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Quick Sort partition function
int partition(int arr[], int low, int high) {
    int pivot = arr[high]; // Taking the last element as pivot
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
        if (arr[j] < pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

// Quick Sort function
void quickSort(int arr[], int low, int high) {
    if (low < high) {
        int pi = partition(arr, low, high);
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}