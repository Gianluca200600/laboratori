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


// Power function in t
// IVP: f(t,z) = -(1-t)^p*z
//      initial value: v
// Solution: u(t) = -exp{(1-t)^{p+1} / (p+1)}

void powert(REAL t, VECTOR *input, VECTOR *output, void *param) {

    POWERt_DATA *data = (POWERt_DATA*) param;
    output->v[0] = -pow(1.0-t, data->p) * input->v[0];

}

void powert_solve(REAL t, VECTOR *input, VECTOR *output, void *param) {

    POWERt_DATA *data = (POWERt_DATA*) param;
    output->v[0] = -exp(pow(1.0-t, 1.0+data->p) / (1.0+data->p));

}


// SIR model
// z = (z1, z2, z3) = (S, I, R)
//          ( -beta*z1*z2            )
// f(t,z) = ( (beta*z1 - gamma) * z2 )
//          ( gamma*z2               )

void sir(REAL t, VECTOR *input, VECTOR* output, void *param) {

    SIR_DATA *data = (SIR_DATA*) param;
    REAL beta = data->beta;
    REAL gamma = data->gamma;
    REAL z1 = input->v[0];
    REAL z2 = input->v[1];
    REAL z3 = input->v[2];

    output->v[0] = -beta*z1*z2;
    output->v[1] = (beta*z1 - gamma)*z2;
    output->v[2] = gamma*z2;

}

void sir_der(REAL t, VECTOR *input, MATRIX *output, void *param) {

    SIR_DATA *data = (SIR_DATA*) param;
    REAL beta = data->beta;
    REAL gamma = data->gamma;
    REAL z1 = input->v[0];
    REAL z2 = input->v[1];
    REAL z3 = input->v[2];
    
    output->A[0][0] = -beta*z2;
    output->A[0][1] = -beta*z1;
    output->A[0][2] = 0.0;

    output->A[1][0] = beta*z2;
    output->A[1][1] = beta*z1 - gamma;
    output->A[1][2] = 0.0;

    output->A[2][0] = 0.0;
    output->A[2][1] = gamma;
    output->A[2][2] = 0.0;

}


// 1D Heat equation (discretized)
// f(t,z) = -M^2*A*z + b where
// b = [1,...1]^T 
//     [ 2 -1  ..       ]
//     [-1  2  -1  ..   ]
// A = [   -1   2  -1 ..] 
//     [   ..  ..  ..   ]
//     [           -1  2]

void heat(REAL t, VECTOR *input, VECTOR *output, void *param) {

    HEAT_DATA *data = (HEAT_DATA*) param;
    REAL m = data->m;

    // components loop
    for(int i=0; i<data->m-1; i++) {

        // diag
        output->v[i] = 2.0*input->v[i];

        if (i>0)
            output->v[i] -= input->v[i-1];
        if (i<data->m-2)
            output->v[i] -= input->v[i+1];

        output->v[i] *= -m*m;
        output->v[i] += 1.0;

    }

}


// Brussellator model for autocatalytic reactions
// f(t,z) = ( A + z1^2*z2 - (B+1)*z1 )
//          ( B*z1 - z1^2*z2 )

void bruss(REAL t, VECTOR *input, VECTOR *output, void *param) {

    BRUSS_DATA *data = (BRUSS_DATA*) param;
    REAL A = data->A;
    REAL B = data->B;
    REAL z1 = input->v[0];
    REAL z2 = input->v[1];

    output->v[0] = A + z1*z1*z2 - (B + 1.0)*z1;
    output->v[1] = B*z1 - z1*z1*z2;

}

void bruss_der(REAL t, VECTOR *input, MATRIX *output, void *param) {

    BRUSS_DATA *data = (BRUSS_DATA*) param;
    REAL A = data->A;
    REAL B = data->B;
    REAL z1 = input->v[0];
    REAL z2 = input->v[1];

    output->A[0][0] = 2*z1*z2 - (B + 1.0);
    output->A[0][1] = z1*z1;
    output->A[1][0] = B - 2*z1*z2;
    output->A[1][1] = -z1*z1;

}


// 1D Allen-Cahn discretized equation
// f(t,z) = -(1/2h)*A*z - F(z)/(eps^2) where
//     [ 1 -1  ..       ]          [ z1^3 -z1 ]
//     [-1  2  -1  ..   ]          [   ....   ]
// A = [   -1   2  -1 ..]  and F = [   ....   ]
//     [   ..  ..  ..   ]          [   ....   ]
//     [           -1  1]          [ zd^3 -zd ]

void allencahn(REAL t, VECTOR *input, VECTOR *output, void* param) {

    ALLENCAHN_DATA *data = (ALLENCAHN_DATA*) param;
    
    int d = input->n + 1;

    for(int i=0; i<d-1; i++) {

        // Arix
        if(i==0 || i==d-2)
            output->v[i] = 1.0 * input->v[i];
        else    
            output->v[i] = 2.0 * input->v[i];

        if(i>0)
            output->v[i] -= input->v[i+1];
        if(i<d-2)
            output->v[i] -= input->v[i+1];
        
        output->v[i] *= -d / 2.0;
        output->v[i] -= (input->v[i] -1.0)*input->v[i] / (data->eps*data->eps);

    }

}

void allencahn_der(REAL t, VECTOR *input, MATRIX *output, void *param) {

    ALLENCAHN_DATA *data = (ALLENCAHN_DATA*) param;

    int d = input->n + 1;

    set_mat(0.0, output);

    //components loop
    for(int i=0; i<d-1; i++) {

        if(i==0 || i==d-2)
            output->A[i][i] = -d/2.0;
        else
            output->A[i][i] = -d/1.0;

        // Off-diagonal entries
        if(i>0)   output->A[i][i-1] = d/2.0;
        if(i<d-2) output->A[i][i+1] = d/2.0;

        // Add nonlinear contribution
        output->A[i][i] -= (3.0*input->v[i]*input->v[i]-1.0) / (data->eps*data->eps);

    }

}


// Prothero_robinson equation
// f(t,z) = lambda(z - sin(2*PI*t)) + 2*PI*cos(2*PI*t)
// Solution: u(t) = exp(lambda*t)v + sin(2*PI*t)

void prorob(REAL t, VECTOR *input, VECTOR *output, void *param) {

    PROROB_DATA *data = (PROROB_DATA*) param;
    REAL lambda = data->lambda;

    output->v[0] = lambda*(input->v[0] - sin(2.0*PI*t)) + 2.0*PI*cos(2.0*PI*t);

}

void prorob_der(REAL t, VECTOR *input, MATRIX *output, void *param) {

    PROROB_DATA *data = (PROROB_DATA*) param;
    REAL lambda = data->lambda;

    output->A[0][0] = lambda;

}

void prorob_solve(REAL t, VECTOR *input, VECTOR *output, void *param) {

    PROROB_DATA *data = (PROROB_DATA*) param;
    REAL lambda = data->lambda;

    output->v[0] = exp(lambda*t)*input->v[0] + sin(2.0*PI*t);

}


// Lotka-Volterra equations
// f(t,z) = ( alpha*z1 - beta*z1*z2 )
//          ( beta*z1*z2 - gamma*z2 )

void lotvol(REAL t, VECTOR *input, VECTOR *output, void *param) {

    LOTVOL_DATA *data = (LOTVOL_DATA*) param;
    REAL alpha = data->alpha;
    REAL beta = data->beta;
    REAL gamma = data->gamma;

    output->v[0] = input->v[0]*(alpha - beta*input->v[1]);
    output->v[1] = (beta*input->v[0] - gamma)*input->v[1];

}

void lotvol_der(REAL t, VECTOR *input, MATRIX *output, void *param) {

    LOTVOL_DATA *data = (LOTVOL_DATA*) param;
    REAL alpha = data->alpha;
    REAL beta = data->beta;
    REAL gamma = data->gamma;

    output->A[0][0] = alpha - beta*input->v[1];
    output->A[0][1] = -beta*input->v[0];

    output->A[1][0] = beta*input->v[1];
    output->A[1][1] = beta*input->v[0] - gamma;

}