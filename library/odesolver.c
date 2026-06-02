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

        lin_comb(-1.0, app_n, 1.0, sol_n, sol_n);
        err = vec_norm2(sol_n);
        max_err = MAX(max_err, err);

    }

    return max_err;

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
        lin_comb(1.0, approx_prev, tau, approx, approx);

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

        lin_comb(1.0, approx_prev, tau / 2.0, K, K);

        data->f(tn + tau / 2.0, K, approx, data->param);

        lin_comb(1.0, approx_prev, tau/2.0, approx, approx);

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
