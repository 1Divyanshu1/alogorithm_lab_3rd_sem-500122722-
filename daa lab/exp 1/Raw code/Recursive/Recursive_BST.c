#include<stdio.h>
#include <stdlib.h>
#include <time.h>

#define SIZE 60000

struct Node* insertRecursive(struct Node* , int);  // Prototype functions for creating and and inserting iteratively in a binary tree.
struct Node* createNode(int);

// Structure for a node in the binary search tree.
struct Node {
    int data;
    struct Node* left;
    struct Node* right;
};

int main(){

    clock_t start_time;  // clock_t is a data type for measuring processor time in clock ticks.
    clock_t end_time;
    double total_time;

    start_time=clock(); 
    
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

    // Insert numbers from the array into the binary search tree iteratively
    for (int i = 0; i < SIZE; i++) {
        root = insertRecursive(root, arr[i]);
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

// Function to insert a node recursively in a binary search tree
struct Node* insertRecursive(struct Node* root, int data) {
    // Base case: If the tree is empty, return the new node as the root
    if (root == NULL) {
        return createNode(data);
    }

    // Recursive case: Traverse the tree to find the correct position
    if (data < root->data) {
        root->left = insertRecursive(root->left, data);  // Insert in the left subtree
    } else if (data > root->data) {
        root->right = insertRecursive(root->right, data); // Insert in the right subtree
    }

    // Return the unchanged root pointer
    return root;
}