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


// Prothero_robinson equation
typedef struct {
    REAL lambda;
}PROROB_DATA;

void prorob(REAL t, VECTOR *input, VECTOR *output, void *param);
void prorob_der(REAL t, VECTOR *input, MATRIX *output, void *param);
void prorob_solve(REAL t, VECTOR *input, VECTOR *output, void *param);


// Lotka-Volterra equations
typedef struct {
    REAL alpha;
    REAL beta;
    REAL gamma;
}LOTVOL_DATA;

void lotvol(REAL t, VECTOR *input, VECTOR *output, void *param);
void lotvol_der(REAL t, VECTOR *input, MATRIX *output, void *param);


// Oregonator
void oregon(REAL t, VECTOR *input, VECTOR *output, void *param);
void oregon_der(REAL t, VECTOR *input, MATRIX *output, void *param);


// Hodgkin-Huxley
typedef struct{
    REAL gNa;
    REAL gK;
    REAL gL;
    REAL vNa;
    REAL vK;
    REAL vL;
}HOD_HUX_DATA;

void hodhux(REAL t, VECTOR *input, VECTOR *output, void *param);
void hodhux_der(REAL t, VECTOR *input, MATRIX *output, void *param);