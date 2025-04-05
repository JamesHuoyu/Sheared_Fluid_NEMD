#define _USE_MATH_DEFINES  // 启用数学常量
#include <math.h>
#include "md.h"
#include <stdio.h>
#include <stdlib.h>

static const real epislon[2][2] = {
    {1.0, 1.5},
    {1.5, 0.5}
};

static const real sigma[2][2] = {
    {1.0, 0.8},
    {0.8, 0.88}
};

static c_real sr2 = 1.0 / (CUTOFF * CUTOFF);
static c_real sr6 = sr2 * sr2 * sr2;
static c_real sr12 = sr6 * sr6;
static c_real pot_cut = sr12 - sr6;

real LJ_cut_force(real r, int type1, int type2){
    if (r < 1.33){
        return 0.0;
    }
    real r2 = r * r;
    real r6 = r2 * r2 * r2;
    real r12 = r6 * r6;

    real epsilon = epislon[type1][type2];
    real sigma_1 = sigma[type1][type2];
    real sigma_2 = sigma_1 * sigma_1;
    real sigma_6 = sigma_2 * sigma_2 * sigma_2;
    real sigma_12 = sigma_6 * sigma_6;
    real force = 24.0 * epsilon / r * (2.0 * sigma_12 / r12 - sigma_6 / r6);

    return (r <= CUTOFF) ? force : 0.0;
}

void compute_forces(Particle *particles, int nparticle, real box_size ,real cutoff, ForceFunc func){
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
            // 处理周期性边界条件
            if (dx > 0.5 * box_size) dx -= box_size;
            if (dx < -0.5 * box_size) dx += box_size;
            if (dy > 0.5 * box_size) dy -= box_size;
            if (dy < -0.5 * box_size) dy += box_size;
            if (dz > 0.5 * box_size) dz -= box_size;
            if (dz < -0.5 * box_size) dz += box_size;
            // 计算距离
            real r = sqrt(dx * dx + dy * dy + dz * dz);
            real f = func(r, particles[i].type, particles[j].type);
            if (f == 0.0) continue; // 跳过零力
            particles[i].fx += f * dx / r;
            particles[i].fy += f * dy / r;
            particles[i].fz += f * dz / r;
            particles[j].fx -= f * dx / r;
            particles[j].fy -= f * dy / r;
            particles[j].fz -= f * dz / r;
        }
    }
}