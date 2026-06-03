/*
** 3.3. Oscillatore armonico
**
** Note:
** - IVP u′′ = −π^2*u  (u(0), u′(0)) = (π/4,0)
** - T = 10
** - N = 1000 passi uniformi
** - Metodo di Eulero Esplicito -> FE (forward)
** - Metodo di Eulero Modificato -> ME
*/

#include "../library/header.h"


FILE *f;

void process(REAL t, REAL tau, VECTOR *approx){

    fprintf(f, "%.5f; %.5f; %.5f\n", t, approx->v[0], approx->v[1]);

}

int main() {

    FILE *f_fe = fopen("l_3_3_fe.txt", "w");
    FILE *f_me = fopen("l_3_3_me.txt", "w");

    int dim = 2;
    REAL t0 = 0.0;
    REAL T = 10.0;
    VECTOR *v = get_vec(dim);
    v->v[0] = PI/4.0;
    v->v[1] = 0.0;

    HARM_OSC_DATA harm_osc_data = {PI};

    ODE_DATA *data = malloc(sizeof(*data));
    data->f = &harm_osc;
    data->N = 1000;
    data->t0 = t0;
    data->T = T;
    data->v = v;
    data->param = &harm_osc_data;
    
    fprintf(f_fe, "%.5f; %.5f; %.5f\n", t0, v->v[0], v->v[1]);
    fe_solve(data, &process);
    fclose(f_fe);

    fprintf(f_me, "%.5f; %.5f; %.5f\n", t0, v->v[0], v->v[1]);
    me_solve(data, &process);
    fclose(f_me);

    free_vec(v);
    free(data);

    return 0;

}