/*
** 7.2. Equazione di Prothero-Robinson.
**
** Notes:
** - u′ = λ(u−g)+g′ con u(0)=v
** - λ = −200
** - g(t) = sin(2πt)
** - T = 1
** - v = 1
** - la soluzione esatta è u(t) = e^(λt)(v − g(0)) + g(t).
*/

#include "../library/header.h"


REAL err;
VECTOR *v;
PROROB_DATA prorob_data;
FILE *f;

void process(REAL t, REAL tau, VECTOR *approx) {

    VECTOR *sol = get_vec(approx->n);
    log_growth_solve(t, v, sol, &prorob_data);
    lin_comb_vec(-1.0, approx, 1.0, sol, sol);
    
    err = MAX(err, nrm2_vec(sol));
    fprintf(f, "%.5f\t %.5f", t, approx->v[0]);

}

int main() {

    f = fopen("l_7_2.txt", "w");

    v = get_vec(1);
    v->v[0] = 1.0;

    REAL t0 = 0.0;
    REAL T = 1.0;

    prorob_data = (PROROB_DATA) {
        -200.0
    };

    ODE_DATA_ERK data = {
        &prorob,
        v,
        t0,
        T,
        60,
        get_rkm_rk4(),
        &prorob_data
    };

    fprintf(f, "%.5f\t %.5f\n", t0, v->v[0]);

    err = 0.0;
    erk_solve(&data, &process);

    // error on terminal
    printf("\nN = %d, err = %.2e\n\n", data.N, err);

    fclose(f);
    free_vec(v);

    return 0;

}
