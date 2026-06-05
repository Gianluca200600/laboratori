/*
** Solution of ordinary differential equations
*/

#include "header.h"


// get_uniform_grid -> Generate uniform grid
VECTOR *get_uniform_grid(REAL t0, REAL T, int N) {

    VECTOR *grid = get_vec(N+1);
    REAL tau = (T - t0) / N;

    for(int n=0; n<N+1; n++)
        grid->v[n] = t0 + tau*n;

    return grid;

}


// max_error_grid -> Max error at grid points
REAL max_error_grid(FUNt sol, MATRIX *approx, VECTOR *grid, int dim, VECTOR *v, void *param) {

    REAL err, max_err = 0.0;
    VECTOR *sol_n = get_vec(dim);
    VECTOR *app_n = get_vec(dim);

    // grid loop
    for(int n=0; n<grid->n; n++) {

        app_n->v = approx->A[n];
        sol(grid->v[n], v, sol_n, param);

        lin_comb_vec(-1.0, app_n, 1.0, sol_n, sol_n);
        err = nrm2_vec(sol_n);
        max_err = MAX(max_err, err);

    }

    return max_err;

}


void cp_from_std(VECTOR *K, VECTOR **k, int s){

    // Space dimension
    int d = k[0]->n;

    // Loop over the entries
    for(int i=0; i<d*s; i++)
        K->v[i] = k[i/d]->v[i%d];
 }


void cp_to_std(VECTOR *K, VECTOR **k, int s){

    // Space dimension
    int d = k[0]->n;

    // Loop over the entries
    for(int i=0; i<d*s; i++)
        k[i/d]->v[i%d] = K->v[i];
 }


// fe_solve -> Forward Euler method
void fe_solve(ODE_DATA *data, APPROX process) {

    VECTOR *approx = get_vec(data->v->n);
    VECTOR *approx_prev = get_vec(data->v->n);
    REAL tau; // step size
    REAL tn = data->t0; // time at the current step

    prod_scal_vec(1.0, data->v, approx);
    tau = (data->T - data->t0) / data->N;

    // time loop
    for(int n=0; n < data->N; n++) {

        prod_scal_vec(1.0, approx, approx_prev);
        data->f(tn, approx_prev, approx, data->param);

        prod_scal_vec(tau, approx_prev, approx_prev);
        lin_comb_vec(1.0, approx_prev, tau, approx, approx);

        tn += tau;

        if(process)
            process(tn, tau, approx);
    }

    free_vec(approx);
    free_vec(approx_prev);

}


// me_solve -> Modified Euler method
void me_solve(ODE_DATA *data, APPROX process) {

    VECTOR *approx = get_vec(data->v->n);
    VECTOR *approx_prev = get_vec(data->v->n);
    VECTOR *K = get_vec(data->v->n);
    REAL tau; // step size
    REAL tn = data->t0; // time at the current step

    prod_scal_vec(1.0, data->v, approx);
    tau = (data->T - data->t0) / data->N;

    // time loop
    for(int n=0; n < data->N; n++) {

        prod_scal_vec(1.0, approx, approx_prev);
        data->f(tn, approx_prev, K, data->param);

        lin_comb_vec(1.0, approx_prev, tau / 2.0, K, K);

        data->f(tn + tau / 2.0, K, approx, data->param);

        lin_comb_vec(1.0, approx_prev, tau/2.0, approx, approx);

        tn += tau;

        if(process)
            process(tn, tau, approx);

    }

    free_vec(approx);
    free_vec(approx_prev);
    free_vec(K);

}


// erk_solve -> Explicit Runge-Kutta method
void erk_solve(ODE_DATA_ERK *data, APPROX process) {

    int d = data->v->n; // dimension of the problem
    int s = data->mtd->b->n; // number of stages
    VECTOR *approx = get_vec(d);
    VECTOR *approx_prev = get_vec(d);
    VECTOR *Y = get_vec(d);

    VECTOR **K = malloc(s*sizeof(VECTOR *));
    for(int i=0; i<s; i++)
        K[i] = get_vec(d);

    REAL tau; // step size
    REAL tn = data->t0; // time at the current step

    prod_scal_vec(1.0, data->v, approx);
    tau = (data->T - data->t0) / data->N;

    // time loop
    for(int n=0; n < data->N; n++) {

        prod_scal_vec(1.0, approx, approx_prev);

        // stages loop
        for(int i=0; i<s; i++) {

            prod_scal_vec(1.0, approx_prev, Y);

            // Previous stages loop
            for(int j=0; j<i; j++)
                lin_comb_vec(tau*data->mtd->A->A[i][j], K[j], 1.0, K[i], Y);

            data->f(tn + tau*data->mtd->c->v[i], Y, K[i], data->param);

            lin_comb_vec(tau*data->mtd->b->v[i], K[i], 1.0, approx, approx);
        }

        tn += tau;

        if(process)
            process(tn, tau, approx);
    }

    free_vec(approx);
    free_vec(approx_prev);
    free_vec(Y);
    for(int i=0; i<s; i++)
        free_vec(K[i]);
    free(K);

}


// Backward Euler
typedef struct {
    FUNt f;
    DFUNt df;
    REAL *tn;
    REAL *tau;
    VECTOR *Un;
    void *param;
}BE_DATA;

void be_fun(VECTOR *input, VECTOR *output, void *param) {

    BE_DATA *data = (BE_DATA*) param;

    VECTOR *z = get_vec(data->Un->n);
    VECTOR *eval = get_vec(data->Un->n);

    lin_comb_vec(1.0, input, 0.0, output, output);
    lin_comb_vec(1.0, data->Un, *data->tau, input, z);

    data->f(*data->tn + *data->tau, z, eval, data->param);

    lin_comb_vec(-1.0, eval, 1.0, output, output);

    free_vec(z);
    free_vec(eval);

}


void be_der(VECTOR *input, MATRIX *output, void *param) {

    BE_DATA *data = (BE_DATA*) param;

    VECTOR *z =get_vec(data->Un->n);
    MATRIX *eval = get_mat(data->Un->n, data->Un->n);

    id_mat(output);

    lin_comb_vec(1.0, data->Un, *data->tau, input, z);

    data->df(*data->tn + *data->tau, z, eval, data->param);
    lin_comb_mat(-(*data->tau), eval, 1.0, output, output);

    free_vec(z);
    free_mat(eval);

}


void be_solve(ODE_DATA_BE *data, APPROX process) {

    VECTOR *approx = get_vec(data->v->n);
    VECTOR* K = get_vec(data->v->n);

    REAL tn = data->t0;
    REAL tau = (data->T - data->t0) / data->N;

    lin_comb_vec(1.0, data->v, 0.0, approx, approx);

    BE_DATA be_data = {
        data->f,
        data->df,
        &tn,
        &tau,
        approx,data->param
    };

    NONLIN_DATA nonlin_data = {
        &be_fun,
        be_der,
        K,
        data->kmax,
        data->tol,
        &be_data
    };

    // time loop
    for(int n=0; n<data->N; n++) {

        data->f(tn, approx, K, data->param);
        newton_solve(&nonlin_data, false);
        lin_comb_vec(tau, K, 1.0, approx, approx);

        tn += tau;

        if(process)
            process(tn, tau, approx);

    }

    free_vec(approx);
    free_vec(K);

}


// Implicit Runge-Kutta method
typedef struct {
    FUNt f; 
    DFUNt Df;
    REAL *tn; 
    REAL *tau; 
    VECTOR *Un;
    RKM *mtd;
    void *param;
}IRK_DATA;

void irk_fun(VECTOR *input, VECTOR *output, void *param) {

    IRK_DATA *data = (IRK_DATA*) param;

    int s = data->mtd->b->n;
    int d = data->Un->n;

    VECTOR *z = get_vec(d);
    VECTOR *eval = get_vec(d);

    lin_comb_vec(1.0, input, 0.0, output, output);

    // stadia loop
    for(int i=0; i<s; i++) {

        lin_comb_vec(1.0, data->Un, 0.0, z, z);

        for(int j=0; j<s*d; j++) 
            z->v[j%d] += *data->tau* data->mtd->A->A[i][j/d] * input->v[j];

        data->f(*data->tn + *data->tau*data->mtd->c->v[i], z, eval, data->param);

        // Add contribution to output
        for(int j=0; j<d; j++)
            output->v[i*d+j] -= eval->v[j];

    }

    free_vec(z);
    free_vec(eval);

}

void irk_der(VECTOR *input, MATRIX *output, void *param) {

    IRK_DATA *data = (IRK_DATA*) param;

    int s = data->mtd->b->n;
    int d =data->Un->n;

    VECTOR *z = get_vec(d);
    MATRIX *eval = get_mat(d, d);

    id_mat(output);

    // stadia loop
    for(int i=0; i<s; i++) {

        lin_comb_vec(1.0, data->Un, 0.0, z, z);

        for(int j=0; j<s*d; j++)
            z->v[j%d] += *data->tau * data->mtd->A->A[i][j/d] * input->v[j];
            
        data->Df(*data->tn + *data->tau * data->mtd->c->v[i], z, eval, data->param);

        for(int j=0; j<d; j++)
            for(int k=0; k<s*d; k++)
                output->A[i*d+j][k] -= *data->tau *data->mtd->A->A[i][k/d] *eval->A[j][k%d];

    }

    free_vec(z);
    free_mat(eval);

}

void irk_solve(ODE_DATA_IRK *data, APPROX process) {

    int s = data->mtd->b->n;
    int d = data->v->n;

    VECTOR *approx = get_vec(d);
    VECTOR *K = get_vec(s*d);
    VECTOR **k = malloc(s*sizeof(VECTOR *));

    for(int i=0; i<s; i++)
        k[i] = get_vec(d);

    REAL tn = data->t0;
    REAL tau = (data->T - data->t0) / data->N;

    lin_comb_vec(1.0, data->v, 0.0, approx, approx);

    IRK_DATA irk_data = {
        data->f,
        data->df,
        &tn,
        &tau,
        approx,
        data->mtd,
        data->param
    };

    NONLIN_DATA non_lin_data = {
        &irk_fun,
        &irk_der,
        K,
        data->kmax,
        data->tol,
        &irk_data
    };

    // time loop
    for(int n=0; n<data->N; n++) {

        // Newton method
        for(int i=0; i<s; i++)
            data->f(tn, approx, k[i], data->param);

        cp_from_std(K, k, s);

        newton_solve(&non_lin_data, false);

        cp_to_std(K, k, s);

        for(int i=0; i<s; i++)
            lin_comb_vec(tau*data->mtd->b->v[i], k[i], 1.0, approx, approx);

        tn += tau;

        if(process)
            process(tn, tau, approx);

    }

    free_vec(approx);
    free_vec(K);
    for(int i=0; i<s; i++)
        free_vec(k[i]);
    free(k);

}