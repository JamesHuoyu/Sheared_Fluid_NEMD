#define _USE_MATH_DEFINES  // 启用数学常量
#include <math.h>
#include "md.h"
#include <stdio.h>
#include <stdlib.h>

c_real epsilon_11 = 1.0;
c_real epsilon_12 = 1.5;
c_real epsilon_22 = 0.5;
c_real sigma_11 = 1.0;
c_real simga_12 = 0.8;
c_real sigma_22 = 0.88;
void LJ_potential_force(real r, real *f, int *type1, int *type2){
    real r2 = r * r;
    real r6 = r2 * r2 * r2;
    real r12 = r6 * r6;

    real epsilon;
    real sigma;
    if (*type1 == 0 && *type2 == 0){
        epsilon = epsilon_11;
        sigma = sigma_11;
    } else if (*type1 == 1 && *type2 == 1){
        epsilon = epsilon_22;
        sigma = sigma_22;
    } else if ((*type1 == 0 && *type2 == 1) || (*type1 == 1 && *type2 == 0)){
        epsilon = epsilon_12;
        sigma = simga_12;
    } else {
        printf("Error: Invalid particle types\n");
        return;
    }
    sigma_2 = sigma * sigma;
    sigma_6 = sigma_2 * sigma_2 * sigma_2;
    sigma_12 = sigma_6 * sigma_6;
    if (r < CUTOFF){
        *f = 24 * epsilon * (2 * sigma_12 / r12 - sigma_6 / r6) / r;
    } else{
        *f = 0.0;
    }
}
void calculate_forces(Particle *particles, int nparticle){
    for(int i = 0; i < nparticle; i++){
        particles[i].fx = 0.0;
        particles[i].fy = 0.0;
        particles[i].fz = 0.0;
    }
    for(int i = 0; i < nparticle; i++){
        for(int j = i + 1; j < nparticle; j++){
            real dx = particles[j].x - particles[i].x;
            real dy = particles[j].y - particles[i].y;
            real dz = particles[j].z - particles[i].z;
            real r = sqrt(dx * dx + dy * dy + dz * dz);
            real f;
            int type1 = particles[i].type;
            int type2 = particles[j].type;
            LJ_potential_force(r, &f, &type1, &type2);
            particles[i].fx += f * dx / r;
            particles[i].fy += f * dy / r;
            particles[i].fz += f * dz / r;
            particles[j].fx -= f * dx / r;
            particles[j].fy -= f * dy / r;
            particles[j].fz -= f * dz / r;
        }
    }
}