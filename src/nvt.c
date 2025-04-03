#include "md.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// 计算系统动能和温度
void compute_thermo(Particle *p, int n, real *ke, real *temp) {
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
        p[i].x += p[i].vx * dt + 0.5 * p[i].fx * dt * dt;
        p[i].y += p[i].vy * dt + 0.5 * p[i].fy * dt * dt;
        p[i].z += p[i].vz * dt + 0.5 * p[i].fz * dt * dt;
        
        p[i].vx += 0.5 * p[i].fx * dt;
        p[i].vy += 0.5 * p[i].fy * dt;
        p[i].vz += 0.5 * p[i].fz * dt;
    }
}

// 速度Verlet积分第二步
void verlet_step2(Particle *p, int n, real dt) {
    for(int i=0; i<n; i++) {
        p[i].vx += 0.5 * p[i].fx * dt;
        p[i].vy += 0.5 * p[i].fy * dt;
        p[i].vz += 0.5 * p[i].fz * dt;
    }
}

// 主模拟循环
void nvt_simulation(Particle *particles, NVT_Params params) {
    FILE *thermo = fopen("thermo.log", "w");
    fprintf(thermo, "Step Temperature KineticEnergy\n");
    
    for(int step=0; step<params.nsteps; step++) {
        // 1. 计算力（此处需要具体实现）
        // compute_forces(particles, params.nparticle);
        
        // 2. 积分运动方程
        verlet_step1(particles, params.nparticle, params.dt);
        // 重新计算力（此处省略）
        // compute_forces(particles, params.nparticle);
        verlet_step2(particles, params.nparticle, params.dt);
        
        // 3. 周期性温度调节
        if(step % params.thermo_freq == 0) {
            remove_momentum(particles, params.nparticle);
            velocity_rescale(particles, params.nparticle, params.target_temp);
        }
        
        // 4. 输出热力学量
        real ke, temp;
        compute_thermo(particles, params.nparticle, &ke, &temp);
        fprintf(thermo, "%d %.4f %.4f\n", step, temp, ke);
    }
    
    fclose(thermo);
}

int main() {
    // 初始化参数
    NVT_Params params = {
        .nparticle = 2916,     // 粒子总数
        .target_temp = 1.0,    // 目标温度
        .dt = 0.001,           // 时间步长
        .nsteps = 10000,       // 总步数
        .thermo_freq = 100     // 每100步调节一次温度
    };
    
    // 分配内存并初始化粒子（此处需与初始化程序对接）
    Particle *particles = malloc(params.nparticle * sizeof(Particle));
    // initialize_from_file(particles, "init.config"); 
    
    // 运行NVT模拟
    nvt_simulation(particles, params);
    
    free(particles);
    return 0;
}
