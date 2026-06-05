/*
** Library
** Functions for (non)linear solvers
*/


// Definitions  
typedef struct{
    int dim; 
    int *perm;
}PERM;

typedef struct{
    FUN F; 
    DFUN DF; 
    VECTOR *sol; 
    int kmax; 
    REAL tol;
    void *param;
}NONLIN_DATA;


// Functions
PERM*  get_perm(int dim);
void  set_perm(PERM *P);
void print_perm(PERM *P);
void  free_perm(PERM *P);

// LU factorization with partial pivoting
void lup(MATRIX *M, PERM *P);

// Gauss method for linear systems
void lup_solve(MATRIX *LU, PERM *P, VECTOR *rhs, VECTOR *sol);

// Newton's method
REAL newton_solve(NONLIN_DATA *data, bool info);

// Chord method
REAL chord_solve(NONLIN_DATA *data, bool info);
