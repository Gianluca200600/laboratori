/*
** 8.2. Un integrale primo quadratico.
*/

#include "../library/header.h"


FILE *f;
HARM_OSC_DATA harm_osc_data;

void process(REAL t, REAL tau, VECTOR *approx) {

    fprintf(f, "%11.6f\t %11.6f\n", t, SQR(harm_osc_data.omega*approx->v[0])/2.0 + SQR(approx->v[1])/2.0);

}

int main() {

    VECTOR *v = get_vec(2);
    v->v[0] = PI / 4.0;
    v->v[1] = 0.0;

    REAL t0 = 0.0;
    REAL T = 10.0 * PI;
    int N = 1000;

    harm_osc_data = (HARM_OSC_DATA) {
        PI
    };

    // FE
    f = fopen("l_8_2_me.txt", "w");

    ODE_DATA_ERK data_me = {
        &harm_osc,
        v,
        t0,
        T,
        N,
        get_rkm_me(),
        &harm_osc_data
    };

    fprintf(f, "%11.6f\t %11.6f\n", t0, SQR(harm_osc_data.omega*v->v[0])/2.0 + SQR(v->v[1])/2.0);
    erk_solve(&data_me, &process);
    fclose(f);

    // RK3
    f = fopen("l_8_2_rk3.txt", "w");

    ODE_DATA_ERK data_rk3 = {
        &harm_osc,
        v,
        t0,
        T,
        N,
        get_rkm_rk3(),
        &harm_osc_data
    };

    fprintf(f, "%11.6f\t %11.6f\n", t0, SQR(harm_osc_data.omega*v->v[0])/2.0 + SQR(v->v[1])/2.0);
    erk_solve(&data_rk3, &process);
    fclose(f);

    // RADAU 2
    f = fopen("l_8_2_ra2.txt", "w");

    ODE_DATA_IRK data_ra2 = {
        &harm_osc,
        &harm_osc_der,
        v,
        t0,
        T,
        N,
        1e-13,
        15,
        get_rkm_ra(2),
        &harm_osc_data
    };

    fprintf(f, "%11.6f\t %11.6f\n", t0, SQR(harm_osc_data.omega*v->v[0])/2.0 + SQR(v->v[1])/2.0);
    irk_solve(&data_ra2, &process);
    fclose(f);

    // Gauss 1
    f = fopen("l_8_2_ga1.txt", "w");

    ODE_DATA_IRK data_ga1 = {
        &harm_osc,
        &harm_osc_der,
        v,
        t0,
        T,
        N,
        1e-13,
        15,
        get_rkm_ga(1),
        &harm_osc_data
    };

    fprintf(f, "%11.6f\t %11.6f\n", t0, SQR(harm_osc_data.omega*v->v[0])/2.0 + SQR(v->v[1])/2.0);
    irk_solve(&data_ga1, &process);
    fclose(f);

    free_vec(v);

    return 0;

}