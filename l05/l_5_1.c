/*
** 5.1. Dominio di stabilità
*/

#include "../library/header.h"


int main() {

    FILE *f = fopen("l_5_1.txt", "w");

    COMPLEX z, R;

    // grid loop
    for(z.Re=-5.0; z.Re<2.0; z.Re+=0.05) {

        for(z.Im=-3.0; z.Im<3.0; z.Im+=0.05) {

            stability_rkm(&z, &R, get_rkm_rk4());

            fprintf(f, "%.6f\t %.6f\t %.6f\t %.6f\n", z.Re, z.Im, R.Re, R.Im);

        }

    }

    fclose(f);

    return 0;

}