/*
** 7.1. Metodi di Radau.
*/

#include "../library/header.h"


REAL err;
VECTOR *v;
LOG_GROWTH_DATA log_growth_data;

void process(REAL t, REAL tau, VECTOR *approx) {

    VECTOR *sol = get_vec(approx->n);
    log_growth_solve(t, v, sol, &log_growth_data);
    lin_comb_vec(-1.0, approx, 1.0, sol, sol);
    
    err = MAX(err, nrm2_vec(sol));

}

int main() {

    FILE *f = fopen("l_7_1.txt", "w");

    v = get_vec(1);
    v->v[0] = 0.1;

    REAL t0 = 0.0;
    REAL T = 10.0;

    log_growth_data = (LOG_GROWTH_DATA) {
        0.5,
        1.0
    };

    ODE_DATA_IRK data = {
        &log_growth,
        &log_growth_der,
        v,
        t0,
        T,
        10,
        NULL,
        1e-14,
        15,
        &log_growth_data
    };

    // steps loop
    for(int k=0; k<=15; k++) {

        fprintf(f, "%6d", data.N);

        // stadia loop
        for(int s=0; s<3; s++) {

            err = 0.0;
            data.mtd = get_rkm_ra(s+1);

        }

    }

    fclose(f);
    free_vec(v);

    return 0;

}