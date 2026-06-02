/*
** Functions for solving nonlinear systems of equations.
*/

#include "header.h"


// Permutations
PERM*  get_perm(int d) {

    PERM *P = malloc(sizeof(PERM));
    P->dim = d;
    P->perm = malloc(d*sizeof(int));

    return P;

}


void  set_perm(PERM *P) {

    for(int i=0; i<P->dim; i++)
        P->perm[i] = i;

}


void print_perm(PERM *P) {
    
    printf("\n");

    for(int i=0; i<P->dim; i++)
        printf("%d\t", P->perm[i]);

    printf("\n");

}


void  free_perm(PERM *P) {

    free(P->perm);
    free(P);
    
}


// LU factorization with partial pivoting
void lup(MATRIX *M, PERM *P) {

    REAL temp;
    int i_max, temp_max;
    int d = M->n_row;

    if(d != M->n_col || d != P->dim) {
        printf("\nlup: incompatible input data\n");
        exit(1);
    }

    set_perm(P);

    // columns loop
    for(int j=0; j<d-1; j++) {

        i_max = j;
        for(int i=j+1; i<d; i++)
            i_max = fabs(M->A[i][j]) > fabs(M->A[i_max][j]) ? i : i_max;

        temp_max = P->perm[j];
        P->perm[j] = P->perm[i_max];
        P->perm[i_max] = temp_max; // update permutation vector

        for(int k=0; k<M->n_col; k++) {

            temp = M->A[j][k];
            M->A[j][k] = M->A[i_max][k];
            M->A[i_max][k] = temp; // exchange rows j and i_max in A

        }

        for(int i=j+1; i<d; i++) {

            M->A[i][j] = M->A[i][j] / M->A[j][j]; // L

            for(int k = j+1; k<d; k++)
                M->A[i][k] = M->A[i][k] - M->A[i][j]*M->A[j][k]; // U

        }

    }

}



// Gauss method for linear systems
void lup_solve(MATRIX *LU, PERM *P, VECTOR *b, VECTOR *x) {

    VECTOR *temp = get_vec(b->n);
    int d = LU->n_row;

    if((d != LU->n_col) || (d != P->dim) || (d != b->n)) {
        printf("\nlup_solve: incompatible input data\n");
        exit(1);
    }

    for(int i=0; i<d; i++)
        temp->v[i] = b->v[P->perm[i]]; // permute entries in b and store in temp

    for(int i=0; i<d; i++) {

        b->v[i] = temp->v[i];
        
        for(int j=0; j<i; j++)
            b->v[i] -= LU->A[i][j]*b->v[j];
    
    } // Solve L*temp = b by forward substitution

    for(int i=d-1; i>=0; i--) { 
    
        b->v[i] /= LU->A[i][i];
        
        for(int j=i+1; j<d; j++)
            b->v[i] -= LU->A[i][j]*b->v[j];
    
    } // Solve U*x = temp by backward substitution

    free_vec(temp);

}


// Newton's method
REAL newton_solve(NONLIN_DATA *data, bool info) {

    int d = data->sol->n;

    MATRIX *DF = get_mat(d, d);
    VECTOR *res = get_vec(d); // residual
    VECTOR *upd = get_vec(d); // update
    PERM *P = get_perm(d);

    data->F(data->sol, res, data->param); // Initial residual
    REAL res_nrm = nrm2_vec(res);

    // Newton iteration
    int k = 0;

    if(info)
        printf("\nNEWTON'S METHOD\n");
    
    while(res_nrm > data->tol && k < data->kmax) {

        if(info)
            printf("iter %2d\t res = %.2e\n", k, res_nrm);

        data->DF(data->sol, DF, data->param); 
        lup(DF, P); 
        lup_solve(DF, P, res, upd); // Solve linear system A*upd = res

        lin_comb_vec(-1.0, upd, 1.0, data->sol, data->sol);
        data->F(data->sol, res, data->param);
        res_nrm = nrm2_vec(res); 

        k++;

    }

    if(info)
        printf("iter %2d\t res = %.2e\n", k, res_nrm);

    free_mat(DF);
    free_vec(res);
    free_vec(upd);
    free_perm(P);

    return res_nrm;

}


// Chord method
REAL chord_solve(NONLIN_DATA *data, bool info) {

    int d = data->sol->n;

    MATRIX *DF = get_mat(d, d);
    VECTOR *res = get_vec(d); // residual
    VECTOR *upd = get_vec(d); // update
    PERM *P = get_perm(d);

    data->F(data->sol, res, data->param); // Initial residual
    REAL res_nrm = nrm2_vec(res);

    if(res_nrm < data->tol)
        return res_nrm;

    data->DF(data->sol, DF, data->param); 
    lup(DF, P);

    // Chord iteration
    int k = 0;

    if(info)
        printf("\nCHORD METHOD\n");

    while(res_nrm > data->tol && k < data->kmax) {

        if(info)
            printf("iter %2d\t res = %.2e\n", k, res_nrm);

        lup_solve(DF, P, res, upd); // Solve linear system A*upd = res
        lin_comb_vec(-1.0, upd, 1.0, data->sol, data->sol); // New approximate solution

        data->F(data->sol, res, data->param); // New residual
        res_nrm = nrm2_vec(res); // New residual norm

        k++;

    }

    if(info)
        printf("iter %2d\t res = %.2e\n", k, res_nrm);

    free_mat(DF);
    free_vec(res);
    free_vec(upd);
    free_perm(P);

    return res_nrm;

}