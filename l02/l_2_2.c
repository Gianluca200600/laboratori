/*
** 2.2. Metodo di Gauß per sistemi lineari
*/

#include "../library/header.h"

int main() {

    // RNG 
	srand(1);
	
	// Memory allocation
	for(int n=2; n<=64; n++) {
		
		int d = n*n;
		
		MATRIX *A = get_mat(d, d);
		VECTOR *exact_sol = get_vec(d);
		VECTOR *approx_sol = get_vec(d);
		VECTOR *b = get_vec(d);
		PERM *p = get_perm(d);
		clock_t time;
		
		// Build A
		set_mat(0.0, A);

		for(int i=0; i<d; i++)
			A->A[i][i] = 4.0; // diagonal
		
		for(int i=1; i<d; i++) {
			A->A[i][i-1] = -1.0; // sub-diagonal
			A->A[i-1][i] = -1.0; // super-diagonal
		}

		for(int i=0; i<d-n; i++) {
			A->A[i][i+n] = 0.0; // n-th super-diagonal
			A->A[i+n][i] = 0.0; // n-th sub-diagonal
		}

		//Build exact solution
		for(int i=0; i<d; i++)
			exact_sol->v[i] = 2.0*(REAL)rand()/RAND_MAX - 1.0; // random in [-1, 1]

		//Compute right-hand side
		prod_mat_vec(A, exact_sol, b);

		// Solve system
		time = -clock();
		lup(A, p);
		lup_solve(A, p, b, approx_sol);
		time += clock();

		// Error
		lin_comb_vec(-1.0, exact_sol, 1.0, approx_sol, approx_sol); // approx_sol = -exact_sol + approx_sol
		REAL err = nrm2_vec(approx_sol);


		// Print
		printf("n = %4d, dim = %4d, time = %6.3f sec, error = %.2e \n", n, d, (REAL) time/CLOCKS_PER_SEC, err);
		

		// Free memory
		free_mat(A);
		free_vec(exact_sol);
		free_vec(approx_sol);
		free_vec(b);
		free(p);
		
	}
	
	return 0;

}