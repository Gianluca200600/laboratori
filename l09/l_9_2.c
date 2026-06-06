/*
** 9.2. Equazioni di Hodgkin-Huxley.
*/

#include "../library/header.h"


FILE *f;
REAL maxV; // max of the potential

void process(REAL t, REAL tau, VECTOR *approx){

    fprintf(f, "%11.6f\t %.4e\t %.4e\t %.4e\t %.4e\n", t, approx->v[0], approx->v[1], approx->v[2], approx->v[3]);
    maxV = MAX(maxV, approx->v[0]);

}

int main() {

    f = fopen("l_9_2.txt", "w");

    VECTOR *v = get_vec(4);
    v->v[0] = 5.0;
    v->v[1] = 5.2934e-2;
    v->v[2] = 5.9611e-1;
    v->v[3] = 3.1768e-1;

    REAL t0 = 0.0;
    REAL T = 10.0;

    HOD_HUX_DATA hod_hux_data = {
        120.0,
        36.0,
        0.3,
        115.0,
        -12.0,
        10.6
    };

    ODE_DATA_ERK data = {
        &hodhux, 
        v, 
        t0,
        T,
        1000,
        get_rkm_rk4(),
        &hod_hux_data
    };

    fprintf(f, "%11.6f\t %.4e\t %.4e\t %.4e\t %.4e\n", t0, v->v[0], v->v[1], v->v[2], v->v[3]);
    
    maxV = v->v[0];
    erk_solve(&data, &process);

    printf("\nMaximum of v: %.6f\n", maxV);

    fclose(f);
    free_vec(v);

    return 0;

}