/*
** 4.1. Crescita logistica con metodi di Runge-Kutta espliciti
**
** Note:
** - IVP u′ = λ*(1− u)u  u(0) = 0.1
** - λ = 0.5
** - T = 10
** - passi: N = 10 * 2^k, k = 0, 1, ..., 15
** - FE
** - ME
** - IE
** - RL
** - RK3
** - RK4
*/

#include "../library/header.h"


LOG_GROWTH_DATA log_growth_data;
VECTOR *v;
REAL err;

// Maximum error
void process(REAL t, REAL tau, VECTOR *approx) {

    VECTOR *sol = get_vec(approx->n);

    log_growth_solve(t, v, sol, &log_growth_data);
    lin_comb_vec(-1.0, approx, 1.0, sol, sol);
    err = MAX(err, nrm2_vec(sol));

}


int main() {

    FILE *f = fopen("l_3_2.txt", "w");

    v = get_vec(1);  
    v->v[0] = 0.1;
    REAL t0 = 0.0;
    REAL T = 10.0;    

    log_growth_data = (LOG_GROWTH_DATA){0.5, 1.0};

    ODE_DATA data = {
        &log_growth,
        v,
        t0,
        T,
        10,
        &log_growth_data,
    };

    // steps loop
    for(int k = 0; k <= 15; k++) {

        // FE
        err = 0.0;
        fe_solve(&data, &process);
        fprintf(f, "%6d\t %.2e\t", data.N, err);

        // ME
        err = 0.0;
        me_solve(&data, &process);
        fprintf(f, "%.2e\n", err);

        // steps update
        data.N *= 2;

    }

    free_vec(v);
    fclose(f);

    return 0;

}