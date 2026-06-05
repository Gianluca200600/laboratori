/*
** 6.3. Equazione di Allen-Cahn.
*/

#include "../library/header.h"


FILE *f;

void process(REAL t, REAL tau, VECTOR *approx) {

    for(int i=0; i<approx->n; i++)
        fprintf(f, "%.6f\t", approx->v[i]);

    fprintf(f, "%.6f\n", t+tau);

}

int main() {

    f = fopen("l_6_3.txt", "w");

    int d = 100;
    REAL t0 = 0.0;
    REAL T = 10.0;

    VECTOR *v = get_vec(d-1);

    for(int i=0; i<d-1; i++)
        v->v[i] = 2.0*(i+1) / ((REAL) d) -1.0;

    ALLENCAHN_DATA allencahn_data = {
        0.01
    };

    ODE_DATA_BE ode_data_be = {
        &allencahn, 
        &allencahn_der, 
        v,
        t0,
        T,
        500,
        1e-10,
        30,
        &allencahn_data
    };

    // initial value
    for(int i=0; i<d-1; i++)
        fprintf(f, "%.6f; ", v->v[i]);
    fprintf(f, "%.6f\n", t0);

    be_solve(&ode_data_be, &process);

    fclose(f);
    free_vec(v);

    return 0;

}
