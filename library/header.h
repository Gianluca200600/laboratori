/* 
** Main Header
*/

// Standard libraries
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>


// Definitions
#define PI 3.141592653589793
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define ABS(a) ((a) >= 0 ? (a) : -(a))
#define SQR(a) ((a) * (a))


// Structures
typedef double REAL;

typedef struct {
    REAL Re;
    REAL Im;
}COMPLEX;

typedef struct {
    REAL *v;
    int n;
}VECTOR;

typedef struct {
    REAL **A;
    int n_row;
    int n_col;
}MATRIX;

typedef void (*FUN) (
    VECTOR *input,
    VECTOR *output,
    void *param
);

typedef void (*DFUN) (
    VECTOR *input,
    MATRIX *output,
    void *param
);

typedef void (*FUNt) (
    REAL t,
    VECTOR *input,
    VECTOR *output,
    void *param
);

typedef void (*DFUNt) (
    REAL t,
    VECTOR *input,
    MATRIX *output,
    void *param
);

typedef struct {
    MATRIX *A;
    VECTOR *b;
    VECTOR *c;
}RKM;


// Libraries
#include "blas.h"
