/*
** 9.1. Oregonator.
** 
** Note:
** - IVP: (u1)' = 77.27(u2 + u1(1 − 8.375*1e^-6u1 −u2))
**        (u2)' = (u3 −(1+u1)u2)/77.27
**        (u3)' = 0.161(u1 − u3)
** - T = 900
** - v = (v1,v2,v3) ∈ R3 | vi ≥ 0
*/

#include "../library/header.h"


FILE *f;

void process(REAL t, REAL tau, VECTOR *approx){

    // Print approximation on file
    fprintf(f, "%11.6f\t %.6e\t %.6e\t %.6e\n", t, approx->v[0], approx->v[1], approx->v[2]);
}

int main() {

    VECTOR *v = get_vec(3);
    v->v[0] = 1.0;
    v->v[1] = 1.0;
    v->v[2] = 1.0;

    REAL t0 = 0.0;
    REAL T = 900.0;

    ODE_DATA_IRK data = {
        &oregon, 
        &oregon_der,
        v, 
        t0,
        T,
        90000,
        1e-13,
        15,
        get_rkm_ra(3),
        NULL
    };

    f = fopen("l_9_1.txt", "w");

    fprintf(f, "%11.6f\t %.6e\t %.6e\t %.6e\n", t0, v->v[0], v->v[1], v->v[2]);
    irk_solve(&data, &process);

    fclose(f);
    free_vec(v);

    return 0;

}