/*
** 4.2. Crescita logistica con metodi di Runge-Kutta espliciti
**
** Note:
** - IVP u′ = -(1-t)^p * u  u(0) = -exp{1/(1+p)}
** - p in {-0.5, 0.5, 1.5, 1500}
** - T = 1 non compreso
** - passi: N = 10 * 2^k, k = 0, 1, ..., 15
** - metodo a scelta
*/

#include "../library/header.h"


REAL err;
VECTOR *v;
POWERt_DATA powert_data;

// Max error
void process(REAL t, REAL tau, VECTOR *approx) {

    VECTOR *sol = get_vec(approx->n);
    powert_solve(t, v, sol, &powert_data);
    lin_comb_vec(-1.0, approx, 1.0, sol, sol);
    err = MAX(err, nrm2_vec(sol));

}

int main() {

    REAL err_prev;

    v = get_vec(1);
    REAL t0 = 0.0;
    REAL T = 0.9999999999;
    REAL p[] = {-0.5, 0.5, 1.5, 1500.0};

    ODE_DATA_ERK data = {
        &powert,
        v,
        t0,
        T,
        10,
        &powert_data,
        get_rkm_rk4()
    };

    // parameters loop
    for(int i=0; i<sizeof(p) / sizeof(p[0]); i++) {

        powert_data.p = p[i];
        v->v[0] = -exp(1.0 / (1.0 + p[i]));

        printf("p = %.1f\n\n", p[i]);
        
        // time steps
        data.N = 10;

        // steps loop
        for(double k = 0; k<=15; k++) {

            err = 0.0;
            erk_solve(&data, &process);

            if (k>0)
                printf("N = %6d\t EOC = %.2f\n", data.N, log(err_prev / err) / log(2));

            err_prev = err;
            data.N *= 2;
        }

        printf("\n\n");

    }

    free_vec(v);

}

