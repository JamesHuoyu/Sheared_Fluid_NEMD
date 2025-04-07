#include "md.h"

int main(){
    SimulationParams params = {
        .nparticle = 256,
        .target_temp = 1.0,
        .dt = 0.005,
        .density = 0.75,
        .nsteps = 1000,
        .nblocks = 1,
        .thermo_freq = 10,
    };
    
    Particle *particles = malloc(params.nparticle * sizeof(Particle));

    initialize_fcc(particles, params.nparticle, params.density, &params.box_size);
    printf("系统初始化完成:\n");
    printf("粒子总数: %d\n", params.nparticle);
    printf("盒子尺寸: %.4f\n", params.box_size);
    printf("类型0粒子: %d\n", (int)(params.nparticle * 0.8));
    printf("类型1粒子: %d\n", params.nparticle - (int)(params.nparticle * 0.8));
    initialize_velocities(particles, params.nparticle, params.target_temp);
    printf("速度初始化完成:\n");
    nvt_simulation(particles, &params);
    printf("NVT模拟完成:\n");
}