/*
** Library
** Concrete examples of IVPs for testing the ODE solvers
*/

// Exponential growth u'=lambda*u  lambda>0
typedef struct{
    REAL lambda;
}EXP_GROWTH_DATA;

void exp_growth(REAL t, VECTOR *input, VECTOR *output, void *param);
void exp_growth_der(REAL t, VECTOR *input, MATRIX *output, void *param);
void exp_growth_solve(REAL t, VECTOR *input, VECTOR *output, void *param);


// Logistic growth u′ = λ*(1 − u)u
typedef struct{
    REAL lambda;
    REAL M;
}LOG_GROWTH_DATA;

void log_growth(REAL t, VECTOR *input, VECTOR *output, void *param);
void log_growth_der(REAL t, VECTOR *input, MATRIX *output, void *param);
void log_growth_solve(REAL t, VECTOR *input, VECTOR *output, void *param);


// Harmonic oscillator u′′ = −ω^2*u
typedef struct{
    REAL omega;
}HARM_OSC_DATA;

void harm_osc(REAL t, VECTOR *input, VECTOR *output, void *param);
void harm_osc_der(REAL t, VECTOR *input, MATRIX *output, void *param);
void harm_osc_solve(REAL t, VECTOR *input, VECTOR *output, void *param);
