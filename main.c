#include "md.h"

int main(){
    SimulationParams params = {
        .nparticle = 2916,
        .target_temp = 1.0,
        .dt = 0.001,
        .density = 1.2,
        .nsteps = 100000,
        .thermo_freq = 100,
    };
    
    Particle *particles = malloc(params.nparticle * sizeof(Particle));

    initialize_fcc(particles, params.nparticle, params.density, &params.box_size);
    printf("系统初始化完成:\n");
    printf("粒子总数: %d\n", params.nparticle);
    printf("盒子尺寸: %.4f\n", params.box_size);
    printf("类型0粒子: %d\n", (int)(params.nparticle * 0.8));
    printf("类型1粒子: %d\n", params.nparticle - (int)(params.nparticle * 0.8));
    initialize_velocities(particles, params.nparticle, params.target_temp);
}