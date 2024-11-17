#include<stdio.h>
#include <stdlib.h>
#include <time.h>

#define SIZE 60000


struct Node* insertIterative(struct Node* , int );  // Prototype functions for creating and and inserting iteratively in a binary tree.
struct Node* createNode(int );

// Structure for a node in the binary search tree
struct Node {
    int data;
    struct Node* left;
    struct Node* right;
};

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

    // Create an empty binary search tree (root node is NULL)
    struct Node* root = NULL;

    start_time=clock();     // Starting the clock to measure time.

    // Insert numbers from the array into the binary search tree iteratively.
    for (int i = 0; i < SIZE; i++) {
        root = insertIterative(root, arr[i]);
    }

    end_time=clock();
    total_time=((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    
    printf("Total elapsed time is %f seconds", total_time);

    
    return 0;
}




// Function to create a new node
struct Node* createNode(int data) {
    struct Node* newNode = (struct Node*)malloc(sizeof(struct Node));
    newNode->data = data;
    newNode->left = NULL;
    newNode->right = NULL;
    return newNode;
}



// Function to insert a node iteratively in a binary search tree
struct Node* insertIterative(struct Node* root, int data) {
    struct Node* newNode = createNode(data);

    if (root == NULL) {
        return newNode; // If the tree is empty, return the new node as root
    }

    struct Node* current = root;
    struct Node* parent = NULL;

    while (current != NULL) {
        parent = current;
        if (data < current->data) {
            current = current->left;
        } else if (data > current->data) {
            current = current->right;
        } else {
            // Duplicates not allowed, return the root unchanged
            return root;
        }
    }

    // Insert the new node at the correct position
    if (data < parent->data) {
        parent->left = newNode;
    } else {
        parent->right = newNode;
    }

    return root; // Return the unchanged root pointer
}


