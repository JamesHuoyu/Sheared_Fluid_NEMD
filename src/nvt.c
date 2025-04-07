#include "md.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// 周期性边界条件
void apply_periodic_boundary(Particle *p, int n, real box_size){
    for(int i=0; i<n; i++){
        if(p[i].x < -box_size/2.0) p[i].x += box_size;
        if(p[i].x > box_size/2.0) p[i].x -= box_size;
        if(p[i].y < -box_size/2.0) p[i].y += box_size;
        if(p[i].y > box_size/2.0) p[i].y -= box_size;
        if(p[i].z < -box_size/2.0) p[i].z += box_size;
        if(p[i].z > box_size/2.0) p[i].z -= box_size;
    }
}

// 计算系统动能和温度
void compute_thermo(const Particle *p, int n, real *ke, real *temp) {
    *ke = 0.0;
    for(int i=0; i<n; i++) {
        *ke += p[i].vx*p[i].vx + p[i].vy*p[i].vy + p[i].vz*p[i].vz;
    }
    *ke *= 0.5; // 动能 = 0.5*m*v² (假设质量m=1)
    *temp = (*ke) * 2.0 / (3.0 * n); // 温度公式 T = (2/3)(KE/N)
}

// 移除系统总动量
void remove_momentum(Particle *p, int n) {
    real px = 0.0, py = 0.0, pz = 0.0;
    
    // 计算总动量
    for(int i=0; i<n; i++) {
        px += p[i].vx;
        py += p[i].vy;
        pz += p[i].vz;
    }
    
    // 计算平均动量并扣除
    px /= n; py /= n; pz /= n;
    for(int i=0; i<n; i++) {
        p[i].vx -= px;
        p[i].vy -= py;
        p[i].vz -= pz;
    }
}

// 速度放缩温度调节
void velocity_rescale(Particle *p, int n, real target_temp) {
    real ke, current_temp;
    compute_thermo(p, n, &ke, &current_temp);
    
    if(current_temp <= 1e-8) return; // 避免除以零
    
    real lambda = sqrt(target_temp / current_temp);
    for(int i=0; i<n; i++) {
        p[i].vx *= lambda;
        p[i].vy *= lambda;
        p[i].vz *= lambda;
    }
}

// 速度Verlet积分第一步
void verlet_step1(Particle *p, int n, real dt) {
    for(int i=0; i<n; i++) {        
        p[i].vx += 0.5 * p[i].fx * dt;
        p[i].vy += 0.5 * p[i].fy * dt;
        p[i].vz += 0.5 * p[i].fz * dt;
    }
}

// 速度Verlet积分第二步
void verlet_step2(Particle *p, int n, real dt) {
    for(int i=0; i<n; i++) {
        p[i].x += p[i].vx * dt;
        p[i].y += p[i].vy * dt;
        p[i].z += p[i].vz * dt;
    }
}

// 主模拟循环
void nvt_simulation(Particle *particles, const SimulationParams *params) {
    EnergyComponents energy;

    FILE *thermo = fopen("thermo.log", "w");
    fprintf(thermo, "Step Temperature KineticEnergy PotentialEnergy PotentialCut ShiftTerm TailCorrection TotalEnergy\n");
    // 计算长程修正项
    energy.tail_corr = calculate_tail_correction(particles, params->nparticle, params->density, params->box_size);

    for(int nblock=0; nblock < params->nblocks; nblock++){
        real ke, temp;
        real ke_sum = 0.0, temp_sum = 0.0, pot_sum = 0.0, total_sum = 0.0;
        real pot_cut_sum = 0.0, shift_sum = 0.0;
        for(int step=0; step<params->nsteps; step++) {
            // 1. 计算力（此处需要具体实现）
            compute_forces_and_energy(particles, params->nparticle, params->box_size, CUTOFF, &energy);
            // 2. 积分运动方程
            verlet_step1(particles, params->nparticle, params->dt);
            verlet_step2(particles, params->nparticle, params->dt);
            // 3. 周期性边界条件
            apply_periodic_boundary(particles, params->nparticle, params->box_size);
            // 4. 计算力（此处需要具体实现）
            compute_forces_and_energy(particles, params->nparticle, params->box_size, CUTOFF, &energy);
            // 5. 积分运动方程
            verlet_step1(particles, params->nparticle, params->dt);
            
            // 周期性温度调节
            if(step % params->thermo_freq == 0) {
                remove_momentum(particles, params->nparticle);
                velocity_rescale(particles, params->nparticle, params->target_temp);
            }
            compute_thermo(particles, params->nparticle, &energy.kinetic, &temp);
            real pot_total = energy.potential_cut + energy.shift_term + energy.tail_corr;
            energy.total = energy.kinetic + pot_total;

            ke_sum += energy.kinetic;
            temp_sum += temp;
            pot_sum += pot_total;
            pot_cut_sum += energy.potential_cut;
            shift_sum += energy.shift_term;
            total_sum += energy.total;
        }
        // 5. 输出热力学量
        real ke_avg = ke_sum / params->nsteps;
        real temp_avg = temp_sum / params->nsteps;
        real pot_avg = pot_sum / params->nsteps;
        real pot_cut_avg = pot_cut_sum / params->nsteps;
        real shift_avg = shift_sum / params->nsteps;
        real total_avg = total_sum / params->nsteps;
        fprintf(thermo, "%d %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n", nblock, temp_avg, ke_avg, pot_avg, pot_cut_avg, shift_avg, energy.tail_corr, total_avg);
    }
    fclose(thermo);
}

