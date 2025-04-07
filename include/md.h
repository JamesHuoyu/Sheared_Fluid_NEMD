#ifndef MD_H
#define MD_H

#define _USE_MATH_DEFINES  // 启用数学常量
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

/*-----------------------------
  基本类型和常量定义
-----------------------------*/
typedef double real;
typedef const double c_real;       // 统一浮点类型
#define CUTOFF 2.5         // 截断半径
#define BOLTZMANN 1.0      // 玻尔兹曼常数（约化单位）

/*-----------------------------
  粒子系统相关结构体
-----------------------------*/
typedef struct {
    real x, y, z;    // 位置分量
    real vx, vy, vz; // 速度分量
    real fx, fy, fz; // 力分量
    int type;        // 粒子类型
} Particle;

typedef struct {
    uint64_t s[2];   // Xorshift128+状态
} xorshift128p_state;

/*-----------------------------
  力场相关结构体
-----------------------------*/
typedef struct {
  real force;
  real potential;
} ForceEnergyPair;

/*--------------------------
函数能量和压力结构体
---------------------------*/
typedef struct{
  real kinetic; // 动能
  real potential_cut; // 截断势能
  real shift_term; // cut-shift项
  real tail_corr; // 长程修正项
  real total; // 总能量
} EnergyComponents;
/*-------------------------------*/
typedef struct{
  real virial; // 压力张量
  real pressure; // 压强
} PressureComponents;

/*-----------------------------
  系统参数结构体
-----------------------------*/
// NVT系综参数
typedef struct {
    int nparticle;     // 总粒子数
    real target_temp;  // 目标温度
    real dt;           // 时间步长
    real density;      // 初始密度
    real box_size;     // 模拟盒子尺寸
    int nsteps;        // 每块模拟步数
    int nblocks;      // 块数
    int thermo_freq;   // 温度调节频率
} SimulationParams;

/*-----------------------------
  函数声明
-----------------------------*/
// 初始化模块
void initialize_fcc(Particle *p, int n, real density, real *box);
void initialize_velocities(Particle *p, int n, real temp);

// 随机数生成
void xorshift128p_init(xorshift128p_state *state, uint64_t seed);
real xorshift128p_double(xorshift128p_state *state);
void gaussian_random(xorshift128p_state *rng, real *v1, real *v2);

// 热力学计算
void compute_thermo(const Particle *p, int n, real *ke, real *temp);
void remove_momentum(Particle *p, int n);

// 温度控制
void velocity_rescale(Particle *p, int n, real target_temp);

// 运动方程积分
void verlet_step1(Particle *p, int n, real dt);
void verlet_step2(Particle *p, int n, real dt);

// 力场计算（不同势函数的接口）
typedef ForceEnergyPair (*ForceFunc)(real r2, int type1, int type2);
void compute_forces_and_energy(Particle *p, int nparticle, real box_size, real cutoff, EnergyComponents *energy);
void compute_forces(Particle *p, int n, real box_size, real cutoff, ForceFunc func);
real calculate_tail_correction(Particle *p, int nparticle, real density, real box_size);

// 主模拟循环
void nvt_simulation(Particle *p, const SimulationParams *params);
void run_md_simulation(Particle *p, const SimulationParams *params);

// 输入输出
void save_config(const Particle *p, int n, const char *filename);
void load_config(Particle **p, int *n, const char *filename);

#endif // MD_H