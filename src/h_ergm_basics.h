#include <R.h>
#include <stdio.h>

double epsilon, maximum;

double Epsilon();
/*
output: smallest double number which can be used in computations
*/

double ln(double x);

double e(double x);

double Stirling(int n);
/*
input: integer n
output: Stirling's approximation of n! on log scale:
- Stirling's approximation: proportional to sqrt(n) (n / e)^n 
- Stirling's approximation on log scale: 0.5 log(n) + n log(n / e) = (n + 0.5) log(n) - n
*/

int* I(int d);
/*
input: dimension of vector
output: vector
*/

int** II(int d1, int d2);
/*
input: order of matrix
output: matrix
*/

double* D(int d);
/*
input: dimension of vector
output: vector
*/

double** DD(int d1, int d2);
/*
input: order of matrix
output: matrix
*/

void Set_I_I(int d, int *vector1, int *vector2);
/*
input: two vectors
output: vector1 = vector2
*/

void Set_II_II(int d1, int d2, int **matrix1, int **matrix2);
/*
input: two matrices
output: matrix 1 = matrix2
*/

void Set_D_D(int d, double *vector1, double *vector2);
/*
input: two vectors
output: vector1 = vector2
*/

void Set_DD_DD(int d1, int d2, double **matrix1, double **matrix2);
/*
input: two matrices
output: matrix 1 = matrix2
*/

void Get_Column(int d, double *vector, double **matrix, int column);
/*
input: dimension of vector, matrix, column 
output: vector
*/ 

void Set_Column(int d, double **matrix, int column, double *vector);
/*
input: dimension of vector, vector, column 
output: vector
*/ 

void Set_Column_Column(int d, double **matrix1, int column, double **matrix2);
/*
input: dimension of columns, matrix, column 
output: matrix
*/ 

void Scale(int d1, int d2, double **matrix1, double s, double **matrix2);
/* 
input: order of matrix, matrix, scale
output: scaled matrix
*/ 

void Print_I(int d, int *vector);

void Print_II(int d1, int d2, int **matrix);

void Print_D(int d, double *vector);

void Print_DD(int d1, int d2, double **matrix);

