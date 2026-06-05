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


// Power function 
typedef struct{
    REAL p;
}POWERt_DATA;

void powert(REAL t, VECTOR *input, VECTOR *output, void *param);
void powert_solve(REAL t, VECTOR *input, VECTOR *output, void *param);


// SIR model
typedef struct {
    REAL beta;
    REAL gamma;
}SIR_DATA;

void sir(REAL t, VECTOR *input, VECTOR* output, void *param);
void sir_der(REAL t, VECTOR *input, MATRIX *output, void *param);


// 1D Heat equation (discretized)
typedef struct {
    int m;
}HEAT_DATA;

void heat(REAL t, VECTOR *input, VECTOR *output, void *param);


// Brussellator model for autocatalytic reactions
typedef struct{
    REAL A;
    REAL B;
}BRUSS_DATA;

void bruss(REAL t, VECTOR *input, VECTOR *output, void *param);
void bruss_der(REAL t, VECTOR *input, MATRIX *output, void *param);


// 1D Allen-Cahn equation (discretized)
typedef struct {
    REAL eps;
}ALLENCAHN_DATA;

void allencahn(REAL t, VECTOR *input, VECTOR *output, void* param);
void allencahn_der(REAL t, VECTOR *input, MATRIX *output, void *param);