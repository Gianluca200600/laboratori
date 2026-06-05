/*
** 6.1. Metodo di Eulero implicito.
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

    FILE *f = fopen("l_6_1.txt", "w");

    v = get_vec(1);  
    v->v[0] = 0.1;
    REAL t0 = 0.0;
    REAL T = 10.0;    

    log_growth_data = (LOG_GROWTH_DATA){0.5, 1.0};

    ODE_DATA ode_data_fe = {
        &log_growth,
        v,
        t0,
        T,
        10,
        &log_growth_data
    };

    ODE_DATA_BE ode_data_be = {
        &log_growth,
        &log_growth_der,
        v,
        t0,
        T,
        10,
        1e-10,
        30,
        &log_growth_data        
    };

    // steps loop
    for(int k = 0; k <= 15; k++) {

        // FE
        err = 0.0;
        fe_solve(&ode_data_fe, &process);
        fprintf(f, "%6d\t %.2e\t", ode_data_fe.N, err);

        // BE
        err = 0.0;
        be_solve(&ode_data_be, &process);
        fprintf(f, "%.2e\n", err);

        // steps update
        ode_data_fe.N *= 2;
        ode_data_be.N *= 2;

    }

    free_vec(v);
    fclose(f);

    return 0;

}