/*
** Concrete examples of IVPs for testing the ODE solvers
*/

#include "header.h"

// Exponential growth 
// IVP: f(t,z) = lambda*z, lambda>0
//      initial value: v
// Solution: u(t) = v*exp(lambda*t)

void exp_growth(REAL t, VECTOR *input, VECTOR *output, void *param) {

    EXP_GROWTH_DATA *data = (EXP_GROWTH_DATA*) param;
    REAL z = input->v[0];
    REAL lambda = data->lambda;

    output->v[0] = lambda * z;

}

void exp_growth_der(REAL t, VECTOR *input, MATRIX *output, void *param) {

    EXP_GROWTH_DATA *data = (EXP_GROWTH_DATA*) param;
    REAL z = input->v[0];
    REAL lambda = data->lambda;

    output->A[0][0] = lambda;

}

void exp_growth_solve(REAL t, VECTOR *input, VECTOR *output, void *param) {

    EXP_GROWTH_DATA *data = (EXP_GROWTH_DATA*) param;
    REAL v = input->v[0];
    REAL lambda = data->lambda;

    output->v[0] = v * exp(lambda * t);

}


// Logistic growth
// IVP: f(t,z) = lambda*(1-z)*z, lambda>0
//      initial value: v
// Solution: u(t) = 1/(1 + (1/v - 1)*exp(-lambda*t))

void log_growth(REAL t, VECTOR *input, VECTOR *output, void *param) {

    LOG_GROWTH_DATA *data = (LOG_GROWTH_DATA*) param;
    REAL z = input->v[0];
    REAL lambda = data->lambda;
    REAL M = data->M;

    output->v[0] = lambda * (1 - z / M) * z;

}

void log_growth_der(REAL t, VECTOR *input, MATRIX *output, void *param) {

    LOG_GROWTH_DATA *data = (LOG_GROWTH_DATA*) param;
    REAL z = input->v[0];
    REAL lambda = data->lambda;
    REAL M = data->M;

    output->A[0][0] = lambda * (1 - 2*z / M);

}

void log_growth_solve(REAL t, VECTOR *input, VECTOR *output, void *param) {

    LOG_GROWTH_DATA *data = (LOG_GROWTH_DATA*) param;
    REAL v = input->v[0];
    REAL lambda = data->lambda;
    REAL M = data->M;

    output->v[0] = v*M*exp(lambda*t) / (M + v*exp(lambda*t) - v);

}


// Harmonic oscillator
// IVP: f(t,z) = [ z_2; -omega^2*z_1 ]
//      initial value: v = [v_1; v_2]
// Solution: u_1(t) = v_1*cos(omega*t) + v_2/omega*sin(omega*t)
//           u_2(t) = -v_1*omega*sin(omega*t) + v_2*cos(omega*t)

void harm_osc(REAL t, VECTOR *input, VECTOR *output, void *param) {

    HARM_OSC_DATA *data = (HARM_OSC_DATA*) param;
    REAL z1 = input->v[0];
    REAL z2 = input->v[1];
    REAL omega = data->omega;

    output->v[0] = z2;
    output->v[1] = -omega * omega * z1;

}

void harm_osc_der(REAL t, VECTOR *input, MATRIX *output, void *param) {

    HARM_OSC_DATA *data = (HARM_OSC_DATA*) param;
    REAL omega = data->omega;

    output->A[0][0] = 0.0;
    output->A[0][1] = 1.0;

    output->A[1][0] = -omega * omega;
    output->A[1][1] = 0.0;

}

void harm_osc_solve(REAL t, VECTOR *input, VECTOR *output, void *param) {

    HARM_OSC_DATA *data = (HARM_OSC_DATA*) param;
    REAL v1 = input->v[0];
    REAL v2 = input->v[1];
    REAL omega = data->omega;

    output->v[0] = v1*cos(omega*t) + (v2 / omega)*sin(omega*t);
    output->v[1] = -v1*omega*sin(omega*t) + v2*cos(omega*t);

}