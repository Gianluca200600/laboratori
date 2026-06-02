/*
** Library
** Functions for Runge-Kutta methods
*/

// Functions

// Stability of RK methods
void stability_rkm(COMPLEX *z, COMPLEX *R, RKM *mtd);

// Generic RK method
RKM *get_rkm(int s);

// Explicit RK methods
RKM *get_rkm_fe();  // Forward Euler
RKM *get_rkm_me();  // Modified Euler
RKM *get_rkm_ie();  // Improved Euler (Heun)
RKM *get_rkm_rl();  // Ralston (Heun)
RKM *get_rkm_rk3(); // 3rd-order explicit Runge-Kutta
RKM *get_rkm_rk4(); // 4th-order explicit Runge-Kutta