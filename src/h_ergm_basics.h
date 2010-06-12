#include <stdio.h>
#include <R.h>
#include <math.h>
#include <Rmath.h>

double** Lower_Half_Matrix(int n);
/* 
input: number of rows of matrix
output: pointer to the upper half of n x n matrix
note: 
- number of elements is n (n - 1) / 2
- elements are 0
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

double** Scale(int d1, int d2, double **matrix, double scale);
/* 
input: order of matrix, matrix, scale
output: scaled matrix
*/ 

void Print_I(int d, int *vector);

void Print_II(int d1, int d2, int **matrix);

void Print_D(int d, double *vector);

void Print_DD(int d1, int d2, double **matrix);

