/*
** 4.3. Un modello epidemiologico
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
** - metodo a scelta
*/

#include "../library/header.h"


FILE *f;

void process(REAL t, REAL tau, VECTOR *approx) {

    fprintf(f, "%11.6f\t %11.6f\t %11.6f\t %11.6f\n", t, approx->v[0], approx->v[1], approx->v[2]);

}

int main() {

    f = fopen("l_4_3.txt", "w");

    VECTOR *v = get_vec(3);
    v->v[0] = 2000.0;
    v->v[1] = 450.0;
    v->v[2] = 50.0;

    REAL t0 = 0.0;
    REAL T = 40.0;

    SIR_DATA sir_data = {2.4e-4, 1.0 / 6.0};
    ODE_DATA_ERK data = {
        &sir,
        v,
        t0,
        T,
        1600,
        &sir_data,
        get_rkm_rk4()
    };

    fprintf(f, "%11.6f\t %11.6f\t %11.6f\t %11.6f\n\n", t0, v->v[0], v->v[1], v->v[2]);

    erk_solve(&data, &process);

    fclose(f);
    free_vec(v);

    return 0;

}