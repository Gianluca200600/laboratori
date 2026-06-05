/*
** 7.3. Equazioni di Lotka-Volterra.
** Notes:
** - IVP: (u1)' = α*u1 −β*u1*u2
**        (u2)' = β*u1*u2 −γ*u2
**   with iv u1(0) = v1 and u2(0) = v2
** - α =1.1
** - β =0.4
** - γ =0.1
** - (v1,v2) = (0.3,3.3)
*/

#include "../library/header.h"


FILE *f;

void process(REAL t, REAL tau, VECTOR *approx){

    fprintf(f, "%.5f\t %.5f\t %.5f\n", t, approx->v[0], approx->v[1]);

}

int main() {

    f = fopen("l_7_3.txt", "w");

    VECTOR *v = get_vec(2);
    v->v[0] = 80.2;
    v->v[1] = 2.0;

    REAL t0 = 0.0;
    REAL T = 1000.0;

    LOTVOL_DATA lotvol_data = {
        1.1,
        0.4,
        0.1
    };

    ODE_DATA_IRK data = {
        &lotvol,
        &lotvol_der,
        v,
        t0,
        T,
        20000,
        1e-14,
        15,
        get_rkm_ra(3),
        &lotvol_data
    };

    irk_solve(&data, &process);

    fclose(f);
    free_vec(v);

    return 0;

}