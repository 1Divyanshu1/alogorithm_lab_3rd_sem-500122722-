
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void printMatrix(int** matrix, int n);
void fillMatrix(int** matrix, int n);
void traditionalMultiply(int** A, int** B, int** C, int n);
void strassenMultiply(int** A, int** B, int** C, int n);
int** allocateMatrix(int n);
void freeMatrix(int** matrix, int n);
void addMatrix(int** A, int** B, int** result, int size);
void subtractMatrix(int** A, int** B, int** result, int size);

int main() {
    int n;
    clock_t end_time, start_time;
    double Traditional_time, Strassen_time;

    printf("Enter the size of the matrix (must be a power of 2): ");
    scanf("%d", &n);

    // Ensure that n is a power of 2
    if ((n & (n - 1)) != 0) {
        printf("Matrix size must be a power of 2!\n");
        return -1;
    }

    srand(time(NULL));  // Seed the random number generator

    int** A = allocateMatrix(n);
    int** B = allocateMatrix(n);
    int** C = allocateMatrix(n);

    fillMatrix(A, n);
    fillMatrix(B, n);

    // Traditional Matrix Multiplication
    start_time = clock();
    traditionalMultiply(A, B, C, n);
    end_time = clock();
    Traditional_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Time taken by Traditional Multiplication: %f seconds\n", Traditional_time);

    // Strassen's Matrix Multiplication
    start_time = clock();
    strassenMultiply(A, B, C, n);
    end_time = clock();
    Strassen_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Time taken by Strassen's Multiplication: %f seconds\n", Strassen_time);

    // Free allocated matrices
    freeMatrix(A, n);
    freeMatrix(B, n);
    freeMatrix(C, n);

    return 0;
}

// Function to allocate a matrix of size n x n
int** allocateMatrix(int n) {
    int** matrix = (int**)malloc(n * sizeof(int*));
    if (matrix == NULL) {
        printf("Memory allocation failed for matrix rows!\n");
        exit(1);
    }
    for (int i = 0; i < n; i++) {
        matrix[i] = (int*)malloc(n * sizeof(int));
        if (matrix[i] == NULL) {
            printf("Memory allocation failed for matrix columns!\n");
            exit(1);
        }
    }
    return matrix;
}

// Function to free the allocated matrix
void freeMatrix(int** matrix, int n) {
    for (int i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// Traditional matrix multiplication
void traditionalMultiply(int** A, int** B, int** C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = 0;
            for (int k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// Add two matrices
void addMatrix(int** A, int** B, int** result, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
}

// Subtract two matrices
void subtractMatrix(int** A, int** B, int** result, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
}

// Strassen's matrix multiplication
void strassenMultiply(int** A, int** B, int** C, int n) {
    // Base case: Use traditional multiplication for small matrices
    if (n <= 128) {  // Increased base case size to 128
        traditionalMultiply(A, B, C, n);
        return;
    }

    int newSize = n / 2;

    // Allocate sub-matrices
    int** A11 = allocateMatrix(newSize);
    int** A12 = allocateMatrix(newSize);
    int** A21 = allocateMatrix(newSize);
    int** A22 = allocateMatrix(newSize);
    int** B11 = allocateMatrix(newSize);
    int** B12 = allocateMatrix(newSize);
    int** B21 = allocateMatrix(newSize);
    int** B22 = allocateMatrix(newSize);

    int** P1 = allocateMatrix(newSize);
    int** P2 = allocateMatrix(newSize);
    int** P3 = allocateMatrix(newSize);
    int** P4 = allocateMatrix(newSize);
    int** P5 = allocateMatrix(newSize);
    int** P6 = allocateMatrix(newSize);
    int** P7 = allocateMatrix(newSize);

    int** C11 = allocateMatrix(newSize);
    int** C12 = allocateMatrix(newSize);
    int** C21 = allocateMatrix(newSize);
    int** C22 = allocateMatrix(newSize);

    int** tempA = allocateMatrix(newSize);
    int** tempB = allocateMatrix(newSize);

    // Dividing matrices into 4 sub-matrices
    for (int i = 0; i < newSize; i++) {
        for (int j = 0; j < newSize; j++) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + newSize];
            A21[i][j] = A[i + newSize][j];
            A22[i][j] = A[i + newSize][j + newSize];

            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + newSize];
            B21[i][j] = B[i + newSize][j];
            B22[i][j] = B[i + newSize][j + newSize];
        }
    }

    // Calculate P1 to P7
    addMatrix(A11, A22, tempA, newSize);
    addMatrix(B11, B22, tempB, newSize);
    strassenMultiply(tempA, tempB, P1, newSize);  // P1 = (A11 + A22) * (B11 + B22)

    addMatrix(A21, A22, tempA, newSize);
    strassenMultiply(tempA, B11, P2, newSize);  // P2 = (A21 + A22) * B11

    subtractMatrix(B12, B22, tempB, newSize);
    strassenMultiply(A11, tempB, P3, newSize);  // P3 = A11 * (B12 - B22)

    subtractMatrix(B21, B11, tempB, newSize);
    strassenMultiply(A22, tempB, P4, newSize);  // P4 = A22 * (B21 - B11)

    addMatrix(A11, A12, tempA, newSize);
    strassenMultiply(tempA, B22, P5, newSize);  // P5 = (A11 + A12) * B22

    subtractMatrix(A21, A11, tempA, newSize);
    addMatrix(B11, B12, tempB, newSize);
    strassenMultiply(tempA, tempB, P6, newSize);  // P6 = (A21 - A11) * (B11 + B12)

    subtractMatrix(A12, A22, tempA, newSize);
    addMatrix(B21, B22, tempB, newSize);
    strassenMultiply(tempA, tempB, P7, newSize);  // P7 = (A12 - A22) * (B21 + B22)

    // Calculate C11, C12, C21, C22
    addMatrix(P1, P4, tempA, newSize);
    subtractMatrix(tempA, P5, tempB, newSize);
    
    addMatrix(tempB, P7, C11, newSize);  // C11
