/*
** 5.3. Brussellator.
*/

#include "../library/header.h"


FILE *f;

void process(REAL t, REAL tau, VECTOR *approx) {

    fprintf(f, "%11.6f\t %11.6f\t %11.6f\n", t, approx->v[0], approx->v[1]);

}

int main() {

    f = fopen("l_5_3.txt", "w");

    VECTOR *v = get_vec(2);
    v->v[0] = 2.0;
    v->v[1] = 2.0;

    REAL t0 = 0.0;
    REAL T = 100.0; 

    BRUSS_DATA bruss_data = {
        30.0,
        30.0
    };

    ODE_DATA_ERK data = {
        &bruss,
        v,
        t0,
        T,
        50000,
        &bruss_data,
        get_rkm_fe()
    };

    fprintf(f, "%11.6f\t %11.6f\t %11.6f\n", t0, v->v[0], v->v[1]);

    erk_solve(&data, &process);  
    
    fclose(f);
    free_vec(v);

    return 0;

}