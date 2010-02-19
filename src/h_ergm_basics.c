#include "h_ergm_basics.h"

double Epsilon()
/*
output: smallest double number which can be used in computations
*/
{
  double e, x;
  x = 1.0;
  do 
    {
    e = x;
    x = x / 2.0;
    }
  while ((1.0 + x) != 1.0);
  return e;
}

double ln(double x)
{
  double y;
  if (x < epsilon) y = log(epsilon);
  else if (x > maximum) y = log(maximum);
  else y = log(x);
  return y;
}

double e(double x)
{
  double y;
  if (x < log(epsilon)) y = epsilon;
  else if (x > log(maximum)) y = maximum;
  else y = exp(x);
  return y;
}

double Stirling(int n)
/*
input: integer n
output: Stirling's approximation of n! on log scale:
- Stirling's approximation: proportional to sqrt(n) (n / e)^n 
- Stirling's approximation on log scale: 0.5 log(n) + n log(n / e) = (n + 0.5) log(n) - n
*/
{
  double s;
  if (n == 0) s = 0; /* log(0!) = log(1) = 0 */
  else s = ((n + 0.5) * ln(n)) - n;
  return s;
}

int* I(int d)
/*
input: dimension of vector
output: vector
*/
{
  int *vector;
  vector = (int*) S_alloc(d,sizeof(int));
  return vector;
}

int** II(int d1, int d2)
/*
input: order of matrix
output: matrix
*/
{
  int i; 
  int **matrix;
  matrix = (int**) S_alloc(d1,sizeof(int*));
  for (i = 0; i < d1; i++)
    {
    matrix[i] = (int*) S_alloc(d2,sizeof(int)); 
    }
  return matrix;
}

double* D(int d)
/*
input: dimension of vector
output: vector
*/
{
  double *vector;
  vector = (double*) S_alloc(d,sizeof(double));
  return vector;
}

double** DD(int d1, int d2)
/*
input: order of matrix
output: matrix
*/
{
  int i; 
  double **matrix;
  matrix = (double**) S_alloc(d1,sizeof(double*));
  for (i = 0; i < d1; i++)
    {
    matrix[i] = (double*) S_alloc(d2,sizeof(double)); 
    }
  return matrix;
}

void Set_I_I(int d, int *vector1, int *vector2)
/*
input: two vectors
output: vector1 = vector2
*/
{
  int i;
  for (i = 0; i < d; i++)
    {
    vector1[i] = vector2[i];
    }
}

void Set_II_II(int d1, int d2, int **matrix1, int **matrix2)
/*
input: two matrices
output: matrix 1 = matrix2
*/
{
  int i, j;
  for (i = 0; i < d1; i++)
    {
    for (j = 0; j < d2; j++)
      {
      matrix1[i] = matrix2[i];
      }
    }
}

void Set_D_D(int d, double *vector1, double *vector2)
/*
input: two vectors
output: vector1 = vector2
*/
{
  int i;
  for (i = 0; i < d; i++)
    {
    vector1[i] = vector2[i];
    }
}

void Set_DD_DD(int d1, int d2, double **matrix1, double **matrix2)
/*
input: two matrices
output: matrix 1 = matrix2
*/
{
  int i, j;
  for (i = 0; i < d1; i++)
    {
    for (j = 0; j < d2; j++)
      {
      matrix1[i] = matrix2[i];
      }
    }
}

void Get_Column(int d, double *vector, double **matrix, int column)
/*
input: dimension of vector, matrix, column 
output: vector
*/ 
{
  int i;
  for (i = 0; i < d; i++)
    {
    vector[i] = matrix[i][column];
    }
}

void Set_Column(int d, double **matrix, int column, double *vector)
/*
input: dimension of vector, vector, column 
output: vector
*/ 
{
  int i;
  for (i = 0; i < d; i++)
    {
    matrix[i][column] = vector[i];
    }
}

void Set_Column_Column(int d, double **matrix1, int column, double **matrix2)
/*
input: dimension of columns, matrix, column 
output: matrix
*/ 
{
  int i;
  for (i = 0; i < d; i++)
    {
    matrix1[i][column] = matrix2[i][column];
    }
}

void Scale(int d1, int d2, double **matrix1, double s, double **matrix2)
/* 
input: order of matrix, matrix, scale
output: scaled matrix
*/ 
{
  int i, j;
  for (i = 0; i < d1; i++)
    {
    for (j = 0; j < d2; j++)
      {
      matrix2[i][j] = matrix1[i][j] * s;
      }
    }
}

void Print_I(int d, int *vector)
{
  int i;
  Rprintf("\n");
  for (i = 0; i < d; i++) Rprintf(" %i",vector[i]);
  Rprintf("\n");
}

void Print_II(int d1, int d2, int **matrix)
{
  int i, j;
  Rprintf("\n");
  for (i = 0; i < d1; i++)
    {
    for (j = 0; j < d2; j++) Rprintf(" %i",matrix[i][j]);
    Rprintf("\n");
    }
}

void Print_D(int d, double *vector)
{
  int i;
  Rprintf("\n");
  for (i = 0; i < d; i++) Rprintf(" %8.4f",vector[i]);
  Rprintf("\n");
}

void Print_DD(int d1, int d2, double **matrix)
{
  int i, j;
  Rprintf("\n");
  for (i = 0; i < d1; i++)
    {
    for (j = 0; j < d2; j++) Rprintf(" %8.4f",matrix[i][j]);
    Rprintf("\n");
    }
}

