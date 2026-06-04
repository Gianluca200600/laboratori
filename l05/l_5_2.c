/*
** 5.2. Limitazione sul passo di metodi espliciti.
*/

#include "../library/header.h"


REAL T;
FILE *f;
REAL tol = 1e-13;

void process(REAL t, REAL tau, VECTOR *approx) {

    if(ABS(t-T) < 1e-13) {

        REAL norm = nrm2_vec(approx);
        fprintf(f, "%.2e\t", norm);

    }

}

int main() {

    f = fopen("l_5_2.txt", "w");

    REAL t0 = 0.0;
    T = 1.0;
    int m = 1;

    HEAT_DATA heat_data = {
        m
    };

    ODE_DATA_ERK data = {
        &heat,
        NULL,
        t0,
        T,
        1,
        &heat_data,
        get_rkm_rk4()
    };

    // dimension loop
    for(int k=0; k<10; k++) {

        heat_data.m ++;

        data.v = get_vec(heat_data.m-1);
        set_vec(0.0, data.v);

        fprintf(f, "%4d \t", heat_data.m-1);

        // option one
        data.N = (int) (4.0*(T-t0)*heat_data.m*heat_data.m);
        erk_solve(&data, &process);

        fprintf(f, "\t");

        // option two
        data.N = (int) ((T-t0)*heat_data.m*heat_data.m);
        erk_solve(&data, &process);

        fprintf(f, "\n");
        free_vec(data.v);
    }

    fclose(f);

    return 0;

}