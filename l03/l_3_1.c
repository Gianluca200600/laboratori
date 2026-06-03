/*
** 3.1. Crescita esponenziale.
**
** Note:
** - IVP u'=lambda*u u(0)=1
** - lambda in [-5, 5]
** - T = 5
** - Metodo di Eulero Esplicito -> FE (forward)
** - Metodo di Eulero Modificato -> ME
*/

#include "../library/header.h"

EXP_GROWTH_DATA exp_growth_data;
VECTOR *v;
REAL T;
REAL err;

// Approximation processing: error at final step
void process(REAL t, REAL tau, VECTOR *approx) {

    // final time check
    if (ABS(t - T) < 1e-13) {

        VECTOR *sol = get_vec(approx->n);

        exp_growth_solve(t, v, sol, &exp_growth_data);
        lin_comb_vec(-1.0, approx, 1.0, sol, sol);
        err = nrm2_vec(sol);

        free_vec(sol);
    }
}

int main() {

    FILE *f = fopen("l_3_1.txt", "w");

    v = get_vec(1);
    REAL t0 = 0.0;
    REAL T = 5.0;
    v->v[0] = 1.0;

    ODE_DATA ode_data = {

        &exp_growth,
        v,
        t0,
        T,
        10000,
        &exp_growth_data

    };

    REAL lambda = exp_growth_data.lambda;

    // growth rate loop
    for (lambda = -5.0; lambda <= 5.0; lambda += 0.1) {

        err = 0.0;
        fe_solve(&ode_data, &process);
        fprintf(f, "%lf %lf\n", lambda, err);

        err = 0.0;
        me_solve(&ode_data, &process);
        fprintf(f, "%lf %lf\n", lambda, err);
    }

    free_vec(v);
    fclose(f);

    return 0;

}