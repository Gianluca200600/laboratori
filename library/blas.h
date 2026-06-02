/*
** Basic Linear Algebra Subprograms (BLAS).
*/

// FUNCTIONS

// Vector
VECTOR *get_vec(const int n);
void set_vec(REAL a, VECTOR *v);
void print_vec(VECTOR *v);
void free_vec(VECTOR *v);
void prod_scal_vec(REAL a, VECTOR *v, VECTOR *out); // out = a*v
void lin_comb_vec(REAL a, VECTOR *v, REAL b, VECTOR *w, VECTOR *out); // out = a*v + b*w
REAL prod_vec_vec(VECTOR *v, VECTOR *w); // v*w (scalar product)
REAL nrm2_vec(VECTOR *v); // Euclidean norm
REAL nrminf_vec(VECTOR *v); // Infinity norm


// Matrix
MATRIX *get_mat(int n_row, int n_col);
void set_mat(REAL a, MATRIX *A);
void print_mat(MATRIX *A);
void free_mat(MATRIX *A);
void id_mat(MATRIX *A);
void prod_scal_mat(REAL a, MATRIX *A, MATRIX *out); // out = a*A
void lin_comb_mat(REAL a, MATRIX *A, REAL b, MATRIX *B, MATRIX *out); // out = a*A + b*B


// Operations
void prod_mat_vec(MATRIX *A, VECTOR *v, VECTOR *out); // out = A*v
void prod_vec_mat(VECTOR *v, MATRIX *A, VECTOR *out); // out = v*A