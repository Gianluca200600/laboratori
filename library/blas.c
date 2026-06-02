/*
** Basic Linear Algebra Subprograms (BLAS).
*/

#include "header.h"


// FUNCTIONS

// Vector
VECTOR *get_vec(const int n) {

    VECTOR *v = malloc(sizeof(VECTOR));
    v->n = n;
    v->v = malloc(n*sizeof(REAL));
    return v;

}


void set_vec(REAL a, VECTOR *v){

    for(int i=0; i<v->n; i++)
        v->v[i] = a;

}


void print_vec(VECTOR *v) {

    for(int i=0; i<v->n; i++)
        printf("%f\n", v->v[i]);

}


void free_vec(VECTOR *v) {

    free(v->v);
    free(v);

}


void prod_scal_vec(REAL a, VECTOR *v, VECTOR *out) {

    for(int i=0; i<v->n; i++)
        out->v[i] = a*v->v[i];

}


void lin_comb_vec(REAL a, VECTOR *v, REAL b, VECTOR *w, VECTOR *out) {
        
    for(int i=0; i<v->n; i++)
        out->v[i] = a*v->v[i] + b*w->v[i];
        
}


REAL prod_vec_vec(VECTOR *v, VECTOR *w) {
    
    // dimension check
    if(v->n != w->n){
        fprintf(stderr, "Incompatible vector dimensions\n");
        exit(1);
    }

    REAL out = 0;

    for(int i=0; i<v->n; i++)
        out += v->v[i] * w->v[i];

    return out;

}

REAL nrm2_vec(VECTOR *v) {

    return sqrt(prod_vec_vec(v, v));

}


REAL nrminf_vec(VECTOR *v) {

    REAL max = 0;

    for(int i=0; i<v->n; i++)
        max = MAX(max, ABS(v->v[i]));

    return max;

}


// Matrix
MATRIX *get_mat(int n_row, int n_col) {

    MATRIX *A = malloc(sizeof(MATRIX));
    A->n_row = n_row;
    A->n_col = n_col;
    A->A = malloc(n_row*sizeof(REAL*));

    for(int i=0; i<n_row; i++)
        A->A[i] = malloc(n_col*sizeof(REAL));

    return A;

}


void set_mat(REAL a, MATRIX *A) {

    for(int i=0; i<A->n_row; i++)
        for(int j=0; j<A->n_col; j++)
            A->A[i][j] = a;

}


void print_mat(MATRIX *A) {
    
    for(int i=0; i<A->n_row; i++){
        for(int j=0; j<A->n_col; j++)
            printf("%f\t", A->A[i][j]);
        printf("\n");
    }

}


void free_mat(MATRIX *A) {

    for(int i=0; i<A->n_row; i++)
        free(A->A[i]);
    free(A->A);
    free(A);

}


void id_mat(MATRIX *A) {

    if(A->n_row != A->n_col){
        fprintf(stderr, "Input matrix is not square\n");
        exit(1);
    }

    set_mat(0.0, A);

    for(int i=0; i<A->n_row; i++)
        A->A[i][i] = 1.0;

}


void prod_scal_mat(REAL a, MATRIX *A, MATRIX *out) {

    for(int i=0; i<A->n_row; i++)
        for(int j=0; j<A->n_col; j++)
            out->A[i][j] = a*A->A[i][j];

}


void lin_comb_mat(REAL a, MATRIX *A, REAL b, MATRIX *B, MATRIX *out) { 
    
    if((A->n_row != B->n_row) || (A->n_col != B->n_col)){
        printf("\naxpby_vec: input matrices with different dimension\n");
        exit(1);
    }

    for(int i=0; i<A->n_row; i++)
        for(int j=0; j<A->n_col; j++)
            out->A[i][j] = a*A->A[i][j] + b*B->A[i][j];

}


// Operations
void prod_mat_vec(MATRIX *A, VECTOR *v, VECTOR *out) {

    if( (v->n != A->n_col) || (out->n != A->n_row)){
        printf("\nprod_mat_vec: input data with different dimension\n");
        exit(1);
    }

    for(int i=0; i<A->n_row; i++)
        for(int j=0; j<A->n_col; j++)
            out->v[i] += A->A[i][j]*v->v[j];

}


void prod_vec_mat(VECTOR *v, MATRIX *A, VECTOR *out) {

    if( (v->n != A->n_row) || (out->n != A->n_col)){
        printf("\nprod_vec_mat: input data with different dimension\n");
        exit(1);
    }

    for(int i=0; i<A->n_col; i++)
        for(int j=0; j<A->n_row; j++)
            out->v[i] += v->v[j]*A->A[j][i];

}