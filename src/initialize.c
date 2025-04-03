#define _USE_MATH_DEFINES  // 启用数学常量
#include <math.h>
#include "md.h"
#include <stdio.h>
#include <stdlib.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


void xorshift128p_init(xorshift128p_state *state, uint64_t seed) {
    state->s[0] = seed;
    state->s[1] = seed ^ 0x6ac744b51297a5c3ULL;
}

real xorshift128p_double(xorshift128p_state *state){
    uint64_t x = state->s[0];
    uint64_t y = state->s[1];
    state->s[0] = y;
    x ^= x << 23;
    state->s[1] = x ^ y ^ (x >> 17) ^ (y >> 26);
    return (real)(state->s[1] + y) / (real)UINT64_MAX;
}

void gaussian_random(xorshift128p_state *rng, real *v1, real *v2) {
    real u1 = xorshift128p_double(rng);
    real u2 = xorshift128p_double(rng);
    real radius = sqrt(-2.0 * log(u1));
    real theta = 2.0 * M_PI * u2;
    *v1 = radius * cos(theta);
    *v2 = radius * sin(theta);
}

void initialize_velocities(Particle *particles, int n, real target_temp){
    xorshift128p_state rng;
    xorshift128p_init(&rng, time(NULL));

    real total_momentum[3] = {0.0, 0.0, 0.0};
    for(int i=0; i < n; i += 2){
        real v1x, v1y, v1z, v2x, v2y, v2z;
        gaussian_random(&rng, &v1x, &v2x);
        gaussian_random(&rng, &v1y, &v2y);
        gaussian_random(&rng, &v1z, &v2z);

        particles[i].vx = v1x;
        particles[i].vy = v1y;
        particles[i].vz = v1z;

        if(i + 1 < n){
            particles[i + 1].vx = v2x;
            particles[i + 1].vy = v2y;
            particles[i + 1].vz = v2z;
        }

       for(int j = i; j < i + 2 && j < n; j++){
            total_momentum[0] += particles[j].vx;
            total_momentum[1] += particles[j].vy;
            total_momentum[2] += particles[j].vz;
        }
    }

    real avg_momentum[3] = {total_momentum[0] / n, total_momentum[1] / n, total_momentum[2] / n};

    for (int i = 0; i < n; i++) {
        particles[i].vx -= avg_momentum[0];
        particles[i].vy -= avg_momentum[1];
        particles[i].vz -= avg_momentum[2];
    }

    real kinetic_energy = 0.0;
    for (int i = 0; i < n; i++) {
        kinetic_energy += 0.5 * (particles[i].vx * particles[i].vx +
                                  particles[i].vy * particles[i].vy +
                                  particles[i].vz * particles[i].vz);
    }

    real scale = sqrt((3.0 * target_temp) / (2.0 * kinetic_energy / n));
    for (int i = 0; i < n; i++) {
        particles[i].vx *= scale;
        particles[i].vy *= scale;
        particles[i].vz *= scale;
    }
}

void validate_system(Particle *particles, int n, real target_temp){
    real momentum[3] = {0.0, 0.0, 0.0};
    real kinetic_energy = 0.0;

    for(int i=0; i < n; i ++){
        momentum[0] += particles[i].vx;
        momentum[1] += particles[i].vy;
        momentum[2] += particles[i].vz;
        kinetic_energy += 0.5 * (particles[i].vx * particles[i].vx +
                                  particles[i].vy * particles[i].vy +
                                  particles[i].vz * particles[i].vz);
    }
    printf( "系统验证:\n");
    printf("总动量: (%.4f, %.4f, %.4f)\n", momentum[0], momentum[1], momentum[2]);
    printf("动能: %.4f\n", kinetic_energy);
}

void initialize_fcc(Particle *particles, int n, real density, real *box) {
    int nc = (int)cbrt(n / 4.0); // FCC单元数每维
    *box = cbrt(n / density);    // 计算盒子尺寸
    real cell_size = *box / nc; // 每个FCC单元的大小
    
    // FCC基元内四个原子的相对位置
    real bases[4][3] = {
        {0.25, 0.25, 0.25}, {0.25, 0.75, 0.75},
        {0.75, 0.25, 0.75}, {0.75, 0.75, 0.25}
    };
    
    int index = 0;
    for (int x = 0; x < nc; x++) {
        for (int y = 0; y < nc; y++) {
            for (int z = 0; z < nc; z++) {
                for (int b = 0; b < 4; b++) {
                    // 计算原子在盒子内的坐标（中心在原点）
                    particles[index].x = (x + bases[b][0]) * cell_size - *box/2.0;
                    particles[index].y = (y + bases[b][1]) * cell_size - *box/2.0;
                    particles[index].z = (z + bases[b][2]) * cell_size - *box/2.0;
                    index++;
                }
            }
        }
    }
    
    // 分配粒子类型（8:2比例）
    int type_a = (int)(n * 0.8);
    for (int i = 0; i < n; i++) {
        particles[i].type = (i < type_a) ? 0 : 1;
    }
}

// int main() {
//     const int n = 2916;        // 总粒子数
//     c_real density = 1.2; // 约化密度
//     real box;                 // 盒子尺寸
    
//     // 分配内存
//     Particle *particles = (Particle*)malloc(n * sizeof(Particle));
//     if (!particles) {
//         fprintf(stderr, "内存分配失败\n");
//         return 1;
//     }
    
//     // 初始化FCC结构和类型
//     initialize_fcc(particles, n, density, &box);
    
//     // 输出验证信息
//     printf("系统初始化完成:\n");
//     printf("粒子总数: %d\n", n);
//     printf("盒子尺寸: %.4f\n", box);
//     printf("类型0粒子: %d\n", (int)(n * 0.8));
//     printf("类型1粒子: %d\n", n - (int)(n * 0.8));
    
//     // 此处可添加文件输出或其他初始化逻辑
//     initialize_velocities(particles, n, 1.0); // 初始化速度
//     validate_system(particles, n, 1.0); // 验证系统

//     free(particles);
//     return 0;
// }