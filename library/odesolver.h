/*
** Library
** Solution of ordinary differential equations
*/

// Approximation
typedef void (*APPROX)(
    REAL t,
    REAL tau,
    VECTOR *approx
);

typedef struct {
    FUNt f;
    VECTOR *v; // initial value
    REAL t0; 
    REAL T;
    int N;
    void *param;
}ODE_DATA;

typedef struct {
    FUNt f;
    VECTOR *v; // initial value
    REAL t0;
    REAL T;
    int N;
    void *param;
    RKM *mtd; // Runge-Kutta method
}ODE_DATA_ERK;


// get_uniform_grid -> Generate uniform grid
VECTOR *get_uniform_grid(REAL t0, REAL T, int N);

// max_error_grid -> Max error at grid points
REAL max_error_grid(FUNt sol, MATRIX *approx, VECTOR *grid, int dim, VECTOR *v, void *param);

// fe_solve -> Forward Euler method
void fe_solve(ODE_DATA *data, APPROX process);

// me_solve -> Modified Euler method
void me_solve(ODE_DATA *data, APPROX process);

// erk_solve -> Explicit Runge-Kutta method
void erk_solve(ODE_DATA_ERK *data, APPROX process);