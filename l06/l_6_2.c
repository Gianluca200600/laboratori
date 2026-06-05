/*
** 6.2. Brussellator-2.
*/

#include "../library/header.h"


FILE *f;
REAL T;
REAL tol = 1e-13;

void process(REAL t, REAL tau, VECTOR *approx) {

    if (ABS(t-T) < tol) 
        fprintf(f, "%11.8f\t %11.8f", approx->v[0], approx->v[1]);

}


int main() {

    f = fopen("l_6_2.txt", "w");

    VECTOR *v = get_vec(2);

    REAL t0 = 0.0;
    T = 100.0;

    v->v[0] = 2.0;
    v->v[1] = 2.0;

    BRUSS_DATA bruss_data = {
        30.0,
        3.0
    };

    ODE_DATA ode_data_fe = {
        &bruss, 
        v, 
        t0, 
        T, 
        100,
        &bruss_data
    };

    ODE_DATA_BE ode_data_be = {
        &bruss, 
        &bruss_der,
        v, 
        t0, 
        T, 
        100,
        1e-10,
        30,
        &bruss_data
    };

    // steps loop
    for(int k=0; k<=12; k++) {

        fprintf(f, "%6d\t", ode_data_fe.N);

        fe_solve(&ode_data_fe, &process);
        fprintf(f, "\t");

        be_solve(&ode_data_be, &process);
        fprintf(f, "\n");

        ode_data_fe.N *=2;
        ode_data_be.N *=2;

    }

    fclose(f);
    free_vec(v);

    return 0;

}