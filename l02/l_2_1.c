/*
** 2.1. Fattorizzazione LU con pivoting parziale
*/

#include "../library/header.h"

int main() {

    // Variables
    int dim = 4;
    MATRIX *A = get_mat(dim, dim);
    PERM *P = get_perm(dim);

    // A
    for(int i=0; i<dim; i++)
        A->A[i] = malloc(dim*sizeof(REAL));
    A->A[0][0] = 2.0; A->A[0][1] = 1.0; A->A[0][2] = 1.0; A->A[0][3] = 0.0;
    A->A[1][0] = 4.0; A->A[1][1] = 3.0; A->A[1][2] = 3.0; A->A[1][3] = 1.0;
    A->A[2][0] = 8.0; A->A[2][1] = 7.0; A->A[2][2] = 9.0; A->A[2][3] = 5.0;
    A->A[3][0] = 6.0; A->A[3][1] = 7.0; A->A[3][2] = 9.0; A->A[3][3] = 8.0;

    // LU factorization with pivoting
    lup(A, P);
    printf("\n\nLU factorization with pivoting:\n\n");
    for(int i=0; i<dim; i++)
        printf(" %d\n", P->perm[i]);
    print_mat(A);

    free_mat(A);
    free_perm(P);

    return 0;

}

