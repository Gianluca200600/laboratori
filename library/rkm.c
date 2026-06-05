/*
** Functions for Runge-Kutta methods
*/

#include "header.h"


// Stability of RK methods
void stability_rkm(COMPLEX *z, COMPLEX *R, RKM *mtd) {

    int s = mtd->b->n;

    MATRIX *A = get_mat(2*s, 2*s);
    VECTOR *b = get_vec(2*s);
    VECTOR *x = get_vec(2*s);
    PERM *p = get_perm(2*s);

    for(int i=0; i<s; i++) {

        for(int j=0; j<s; j++) {

            A->A[i][j] = - z->Re * mtd->A->A[i][j];
            A->A[i][j+s] = z->Im * mtd->A->A[i][j];
            A->A[i+s][j] = - z->Im * mtd->A->A[i][j];
            A->A[i+s][j+s] = z->Re * mtd->A->A[i][j];

        }

        A->A[i][i] += 1.0;
        A->A[i+s][i+s] += 1.0; 

        b->v[i] = z->Re;
        b->v[i+s] = z->Im;

    }

    // Solve A*x=b
    lup(A, p);
    lup_solve(A, p, b, x);

    // Stability: Real part
    R->Re = 1.0;

    for(int l=0; l<s; l++)
        R->Re += mtd->b->v[l] * x->v[l];

    // Stability: Imaginary part
    R->Im = 0.0;

    for(int k=0; k<s; k++)
        R->Im += mtd->b->v[k] * x->v[k+s];
    
    free_mat(A);
    free_vec(b);
    free_vec(x);
    free(p);

}


// Generic RK method with s stages, Butcher table 
RKM *get_rkm(int s) {

    RKM *mtd = malloc(sizeof(RKM));

    mtd->A = get_mat(s, s);
    mtd->b = get_vec(s);
    mtd->c = get_vec(s);

    return mtd;

}


// Explicit RK methods

// Forward Euler
RKM *get_rkm_fe() {

    RKM *mtd = get_rkm(1);

    mtd->A->A[0][0] = 0.0;
    mtd->b->v[0] = 1.0;
    mtd->c->v[0] = 0.0;

    return mtd;

}


// Modified Euler (midpoint)
RKM *get_rkm_me() {

    RKM *mtd = get_rkm(2);

    mtd->A->A[0][0] = 0.0;
    mtd->A->A[0][1] = 0.0;
    mtd->A->A[1][0] = 0.5;
    mtd->A->A[1][1] = 0.0;

    mtd->b->v[0] = 0.0;
    mtd->b->v[1] = 1.0;

    mtd->c->v[0] = 0.0;
    mtd->c->v[1] = 0.5;

    return mtd;

}


// Heun: Improved Euler
RKM *get_rkm_ie() {

    RKM *mtd = get_rkm(2);

    mtd->A->A[0][0] = 0.0;
    mtd->A->A[0][1] = 0.0;
    mtd->A->A[1][0] = 1.0;
    mtd->A->A[1][1] = 0.0;

    mtd->b->v[0] = 0.0;
    mtd->b->v[1] = 1.0;

    mtd->c->v[0] = 0.5;
    mtd->c->v[1] = 0.5;

    return mtd;

}


// Heun: Ralston
RKM *get_rkm_rl() {

    RKM *mtd = get_rkm(2);

    mtd->A->A[0][0] = 0.0;
    mtd->A->A[0][1] = 0.0;
    mtd->A->A[1][0] = 2.0 / 3.0;
    mtd->A->A[1][1] = 0.0;

    mtd->b->v[0] = 0.0;
    mtd->b->v[1] = 2.0 / 3.0;

    mtd->c->v[0] = 0.25;
    mtd->c->v[1] = 0.75;

    return mtd;

}


// 3rd-order explicit Runge-Kutta
RKM *get_rkm_rk3() {

    RKM *mtd = get_rkm(3);

    mtd->A->A[0][0] = 0.0;
    mtd->A->A[0][1] = 0.0;
    mtd->A->A[0][2] = 0.0;
    mtd->A->A[1][0] = 0.5;
    mtd->A->A[1][1] = 0.0;
    mtd->A->A[1][2] = 0.0;
    mtd->A->A[2][0] = -1.0;
    mtd->A->A[2][1] = 2.0;
    mtd->A->A[2][2] = 0.0;

    mtd->b->v[0] = 1.0 / 6.0;
    mtd->b->v[1] = 2.0 / 3.0;
    mtd->b->v[2] = 1.0 / 6.0;

    mtd->c->v[0] = 0.0;
    mtd->c->v[1] = 0.5;
    mtd->c->v[2] = 1.0;

    return mtd;

}


// 4th-order explicit Runge-Kutta
RKM *get_rkm_rk4() {

    RKM *mtd = get_rkm(4);

    mtd->A->A[0][0] = 0.0;
    mtd->A->A[0][1] = 0.0;
    mtd->A->A[0][2] = 0.0;
    mtd->A->A[0][3] = 0.0;

    mtd->A->A[1][0] = 0.5;
    mtd->A->A[1][1] = 0.0;
    mtd->A->A[1][2] = 0.0;
    mtd->A->A[1][3] = 0.0;

    mtd->A->A[2][0] = 0.0;
    mtd->A->A[2][1] = 0.5;
    mtd->A->A[2][2] = 0.0;
    mtd->A->A[2][3] = 0.0;

    mtd->A->A[3][0] = 0.0;
    mtd->A->A[3][1] = 0.0;
    mtd->A->A[3][2] = 1.0;
    mtd->A->A[3][3] = 0.0;


    mtd->b->v[0] = 1.0 / 6.0;
    mtd->b->v[1] = 1.0 / 3.0;
    mtd->b->v[2] = 1.0 / 3.0;
    mtd->b->v[3] = 1.0 / 6.0;

    mtd->c->v[0] = 0.0;
    mtd->c->v[1] = 0.5;
    mtd->c->v[2] = 0.5;
    mtd->c->v[3] = 1.0;

    return mtd;

}


// Radau
RKM *get_rkm_ra1() {

    RKM *mtd = get_rkm(1);

    mtd->A->A[0][0] = 1.0;
    mtd->b->v[0] = 1.0;
    mtd->c->v[0] = 1.0;

    return mtd;

}

RKM *get_rkm_ra2() {

    RKM *mtd = get_rkm(2);

    mtd->A->A[0][0] = 5.0/12.0; 
    mtd->A->A[0][1] = -1.0/12.0;
    mtd->A->A[1][0] = 3.0/4.0;
    mtd->A->A[1][1] = 1.0/4.0;

    // Fill b
    mtd->b->v[0] = 3.0/4.0;
    mtd->b->v[1] = 1.0/4.0;

    // Vector c
    mtd->c->v[0] = 1.0/3.0;
    mtd->c->v[1] = 1.0;

    return mtd;

}

RKM *get_rkm_ra3() {

    REAL s6 = sqrt(6.0);
    RKM *mtd = get_rkm(3);

    mtd->A->A[0][0] = (+88.0 - 7.0*s6) / 360.0;
    mtd->A->A[0][1] = (+296.0 - 169.0*s6) / 1800.0;
    mtd->A->A[0][2] = (-2.0 + 3.0*s6) / 225.0;

    mtd->A->A[1][0] = (+296.0 + 169.0*s6) / 1800.0;
    mtd->A->A[1][1] = (+88.0 + 7.0*s6) / 360.0;
    mtd->A->A[1][2] = (-2.0 - 3.0*s6) / 225.0;

    mtd->A->A[2][0] = (+16.0 - 1.0*s6) / 36.0;
    mtd->A->A[2][1] = (+16.0 + 1.0*s6) / 36.0;
    mtd->A->A[2][2] = (+1.0) / 9.0;

    // Fill b
    mtd->b->v[0] = (16.0 - s6) / 36.0;
    mtd->b->v[1] = (16.0 + s6) / 36.0;
    mtd->b->v[2] = 1.0 / 9.0;

    // Fill c
    mtd->c->v[0] = (4.0 - s6) / 10.0;
    mtd->c->v[1] = (4.0 + s6) / 10.0;
    mtd->c->v[2] = 1.0;

    return mtd;

}

RKM *get_rkm_ra(int s) {

    switch(s) {

        case 1: return get_rkm_ra1();
        case 2: return get_rkm_ra2();
        case 3: return get_rkm_ra3();

        default: printf("\nget_rkm_ra: no Radau method with %d stadia is implemented\n", s);
        exit(1);

    }

}


// Gauss methods
RKM *get_rkm_ga1() {

    RKM *mtd = get_rkm(1);

    mtd->A->A[0][0] = 1.0 / 2.0;
    mtd->b->v[0] = 1.0;
    mtd->c->v[0] = 1.0 / 2.0;

    return mtd;

}

RKM *get_rkm_ga2(){

    // Get tableau
    RKM *mtd = get_rkm(2);

    // Fill A
    mtd->A->A[0][0] = 1.0 / 4.0;  
    mtd->A->A[0][1] = 1.0 / 4.0 - sqrt(3.0) / 6.0;
    mtd->A->A[1][0] = 1.0 / 4.0 + sqrt(3.0) / 6.0;
    mtd->A->A[1][1] = 1.0 / 4.0;

    // Fill b
    mtd->b->v[0] = 1.0 / 2.0;
    mtd->b->v[1] = 1.0 / 2.0;

    // Vector c
    mtd->c->v[0] = 1.0 / 2.0 - sqrt(3.0) / 6.0;
    mtd->c->v[1] = 1.0 / 2.0 + sqrt(3.0) / 6.0;

    return mtd;
}

RKM *get_rkm_ga3(){

    REAL s15 = sqrt(15);

    // Get tableau
    RKM *mtd = get_rkm(3);

    // Fill A
    mtd->A->A[0][0] = 5.0 / 36.0;
    mtd->A->A[0][1] = 8.0 / 36.0 - s15 / 15.0;
    mtd->A->A[0][2] = 5.0 / 36.0 - s15 / 30.0;

    mtd->A->A[1][0] = 5.0 / 36.0 + s15 / 24.0;
    mtd->A->A[1][1] = 8.0 / 36.0;
    mtd->A->A[1][2] = 5.0 / 36.0 - s15 / 24.0;

    mtd->A->A[2][0] = 5.0 / 36.0 + s15 / 30.0;
    mtd->A->A[2][1] = 8.0 / 36.0 + s15 / 15.0;
    mtd->A->A[2][2] = 5.0 / 36.0;

    // Fill b
    mtd->b->v[0] = 5.0 / 18.0;
    mtd->b->v[1] = 8.0 / 18.0;
    mtd->b->v[2] = 5.0 / 18.0;

    // Fill c
    mtd->c->v[0] = 1.0 / 2.0 - s15 / 10.0;
    mtd->c->v[1] = 1.0 / 2.0;
    mtd->c->v[2] = 1.0 / 2.0 + s15 / 10.0;

    return mtd;
}

RKM *get_rkm_ga(int s){

    switch(s){

        case 1: return get_rkm_ga1();
        case 2: return get_rkm_ga2();
        case 3: return get_rkm_ga3();

        default: printf("\nget_rkm_ga: no Gauss method with %d stadia is implemented\n", s);
        exit(1);
    }
    
 }