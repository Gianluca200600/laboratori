/*
** 8.1. Un integrale primo lineare.
**
** Note:
** - IVP 
** s' = -beta * s * i
** i' = (beta*s - gamma) * i
** r' = gamma * i
** - Dati iniziali: s(0) = 2000
**                  i(0) = 450
**                  r(0) = 50
** - T = 40 non compreso
** - param: beta = 2.4*e^-4
**          gamma = 1.0 / 6.0
** - Metodi: EE, R1, G1
*/

#include "../library/header.h"


FILE *f;
VECTOR *v;

void process(REAL t, REAL tau, VECTOR *approx) {

    REAL e0 = v->v[0] + v->v[1] + v->v[2];
    REAL en = approx->v[0] + approx->v[1] + approx->v[2];

    fprintf(f, "%11.6f\t %.2e\n", t, (e0 - en) / e0);

}

int main() {

    v = get_vec(3);
    v->v[0] = 2000.0;
    v->v[1] = 450.0;
    v->v[2] = 50.0;

    REAL t0 = 0.0;
    REAL T = 40.0;

    SIR_DATA sir_data = {2.4e-4, 1.0 / 6.0};

    // FE
    f = fopen("l_8_1_fe.txt", "w");

    ODE_DATA_ERK data_fe = {
        &sir,
        v,
        t0,
        T,
        100,
        get_rkm_fe(),
        &sir_data
    };

    fprintf(f, "%11.6f\t %.2e\n\n", t0, 0.0);
    erk_solve(&data_fe, &process);

    fclose(f);

    // RADAU 1
    f = fopen("l_8_1_ra1.txt", "w");

    ODE_DATA_IRK data_ra1 = {
        &sir,
        &sir_der,
        v,
        t0,
        T,
        100,
        1e-13,
        15,
        get_rkm_ra(1),
        &sir_data
    };

    fprintf(f, "%11.6f\t %.2e\n\n", t0, 0.0);
    irk_solve(&data_ra1, &process);

    fclose(f);

    // GAUSS 1
    f = fopen("l_8_1_ga1.txt", "w");

    ODE_DATA_IRK data_ga1 = {
        &sir,
        &sir_der,
        v,
        t0,
        T,
        100,
        1e-13,
        15,
        get_rkm_ga(1),
        &sir_data
    };

    fprintf(f, "%11.6f\t %.2e\n\n", t0, 0.0);
    irk_solve(&data_ga1, &process);
    
    fclose(f);

    free_vec(v);

    return 0;

}
