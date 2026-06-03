/*
** 2.3. Implementazione di solutori non lineari
*/

#include "../library/header.h"


// Target function f
void f(VECTOR *input, VECTOR *output, void *param) {

    REAL x1 = input->v[0];
    REAL x2 = input->v[1];

    output->v[0] = atan(x1) + x1*x2;
    output->v[1] = x2*(2.0 + cos(x2)) + x1*x1;

}

// Derivative of f
void df(VECTOR *input, MATRIX *output, void *param) {

    REAL x1 = input->v[0];
    REAL x2 = input->v[1];

    output->A[0][0] = 1.0/(1.0 + x1*x1) + x2;
    output->A[0][1] = x1;
    output->A[1][0] = 2* x1;
    output->A[1][1] = 2.0 + cos(x2) - x2*sin(x2);

}


int main() {

    // random seed
    srand(1);

    int dim = 2;

    // starting point
    REAL theta = (REAL)rand() / RAND_MAX;
    REAL rho = 0.2;

    VECTOR *sol = get_vec(dim);
    sol->v[0] = rho*cos(2*PI*theta);
    sol->v[1] = rho*sin(2*PI*theta);

    NONLIN_DATA *data = malloc(sizeof(NONLIN_DATA));
    data->F = &f;
    data->DF = &df;
    data->param = NULL;
    data->sol = sol;
    data->kmax = 100;
    data->tol = 1e-10;

    REAL res_nrm_newton = newton_solve(data, true);
    // REAL res_nrm_chord = chord_solve(data, true);

    printf("\nApproximate solution:\t");
    print_vec(sol);
    printf("\nFinal residual = %.2e\n", res_nrm_newton);

    free(data);

    return 0;

}