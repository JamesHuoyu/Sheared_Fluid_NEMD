#define _USE_MATH_DEFINES  // 启用数学常量
#include <math.h>
#include "md.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef NEAREST_DISTANCE_SQUARE
#define NEAREST_DISTANCE_SQUARE 1.77
#endif

static const real epsilon[2][2] = {
    {1.0, 1.0},
    {1.0, 1.0}
};

static const real sigma[2][2] = {
    {1.0, 1.0},
    {1.0, 1.0}
};

static c_real sr2 = 1.0 / (CUTOFF * CUTOFF);
static c_real sr6 = sr2 * sr2 * sr2;
static c_real sr12 = sr6 * sr6;
static c_real pot_cut = sr12 - sr6;

ForceEnergyPair LJ_cut_force_energy(real r, int type1, int type2){
    ForceEnergyPair result = {0.0, 0.0};
    real rc = 1.0 / (r * r);
    bool overlap = rc > NEAREST_DISTANCE_SQUARE;
    if(overlap){
        printf("粒子重叠，距离: %.4f, 最近距离: %.4f\n", r, sqrt(1.0 / NEAREST_DISTANCE_SQUARE));
    }
    // if(r < nearest_distance){
    //     result.force = 0.0;
    //     result.potential = 0.0;
    //     return result;
    // }

    real r2 = r * r;
    real r6 = r2 * r2 * r2;
    real r12 = r6 * r6;

    real eps = epsilon[type1][type2];
    real sig = sigma[type1][type2];
    real sig2 = sig * sig;
    real sig6 = sig2 * sig2 * sig2;
    real sig12 = sig6 * sig6;

    // 计算力
    result.force = 24.0 * eps / r * (2.0 * sig12 / r12 - sig6 / r6);
    // 计算势能
    result.potential = 4.0 * eps * (sig12 / r12 - sig6 / r6);

    if (r > CUTOFF) {
        result.force = 0.0;
        result.potential = 0.0;
    }

    return result;
}

real calculate_shift_term(int type1, int type2){
    real eps = epsilon[type1][type2];
    real sig = sigma[type1][type2];
    real sig2 = sig * sig;
    real sig6 = sig2 * sig2 * sig2;
    real sig12 = sig6 * sig6;

    return 4.0 * eps * (sig12 * sr12 - sig6 *sr6);
}

void compute_forces_and_energy(Particle *particles, int nparticle, real box_size ,real cutoff, EnergyComponents *energy){
    energy->potential_cut = 0.0;
    energy->potential_shifted = 0.0;
    // 初始化力
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
            while(dx > box_size / 2.0) dx -= box_size;
            while(dx < -box_size / 2.0) dx += box_size;
            while(dy > box_size / 2.0) dy -= box_size;
            while(dy < -box_size / 2.0) dy += box_size;
            while(dz > box_size / 2.0) dz -= box_size;
            while(dz < -box_size / 2.0) dz += box_size;
            // 计算距离
            real r = sqrt(dx * dx + dy * dy + dz * dz);
            ForceEnergyPair pair = LJ_cut_force_energy(r, particles[i].type, particles[j].type);
            real f = pair.force;
            if (f == 0.0) continue; // 跳过零力
            particles[i].fx += f * dx / r;
            particles[i].fy += f * dy / r;
            particles[i].fz += f * dz / r;
            particles[j].fx -= f * dx / r;
            particles[j].fy -= f * dy / r;
            particles[j].fz -= f * dz / r;
            if (r <= CUTOFF){
                energy->potential_cut += pair.potential;

                real shift = calculate_shift_term(particles[i].type, particles[j].type);
                energy->potential_shifted += (pair.potential - shift);
            }
        }
    }
}

real calculate_tail_correction(Particle *particles, int nparticle, real density, real box_size){
    // 粒子类型分类计算
    int type_count[2] = {0, 0};

    for(int i = 0; i < nparticle; i++){
        if(particles[i].type < 2){
            type_count[particles[i].type]++;
        }
    }

    real tail_corr = 0.0;

    for(int type1 = 0; type1 < 2; type1++){
        for(int type2 = 0; type2 < 2; type2 ++){
            real n1 = type_count[type1];
            real n2 = type_count[type2];
            real pairs;

            if (type1 == type2) {
                pairs = n1 * (n1 - 1) / 2.0; // 相同类型的粒子对数目
            } else {
                pairs = n1 * n2; // 不同类型的粒子对数目
            }

            if (pairs > 0) {
                real eps = epsilon[type1][type2];
                real sig = sigma[type1][type2];
                real sig3 = sig * sig * sig;
                real sig6 = sig3 * sig3;
                real sig12 = sig6 * sig6;

                real rc = CUTOFF;
                real rc3 = rc * rc * rc;
                real rc9 = rc3 * rc3 * rc3;

                real term1 = sig12 / (9.0 * rc9);
                real term2 = sig6 / (3.0 * rc3);
                real correction = (8.0 / 3.0) * M_PI * pairs * density * eps * (term1 - term2);
                tail_corr += correction;
            }
        }
    }
    return tail_corr;
}