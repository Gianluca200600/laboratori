/*
** 8.3. Un integrale primo non polinomiale.
*/

#include "../library/header.h"


FILE *f;
LOTVOL_DATA lotvol_data;

void process(REAL t, REAL tau, VECTOR *approx) {

    fprintf(f, "%11.6f\t %16.11f\n", t, lotvol_data.beta*(approx->v[0]+approx->v[1])
                                       -lotvol_data.gamma*log(approx->v[0])
                                       -lotvol_data.alpha*log(approx->v[1]));

}

int main() {

    f = fopen("l_8_3.txt", "w");

    VECTOR *v = get_vec(2);
    v->v[0] = 0.0;
    v->v[1] = 3.3;

    REAL t0 = 0.0;
    REAL T = 100.0;

    lotvol_data = (LOTVOL_DATA) {
        1.1,
        0.4,
        0.1
    };

    // Gauss 3
    ODE_DATA_IRK data_ga3 = {
        &lotvol, 
        &lotvol_der, 
        v, 
        t0, 
        T, 
        100, 
        1e-13,
        15,
        get_rkm_ra(3),
        &lotvol_data
    };

    fprintf(f, "%11.6f\t %16.11f\n", t0, lotvol_data.beta*(v->v[0]+v->v[1])
                                       -lotvol_data.gamma*log(v->v[0])
                                       -lotvol_data.alpha*log(v->v[1]));

    irk_solve(&data_ga3, &process);
    
    fclose(f);
    free_vec(v);

    return 0;

}