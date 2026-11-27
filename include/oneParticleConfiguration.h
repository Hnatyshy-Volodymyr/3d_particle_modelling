#ifndef ONE_PARTICLE_CONFIGURATION_H
#define ONE_PARTICLE_CONFIGURATION_H

#include "particle.h"

#define TIME_STEP 0.01
#define TOTAL_STEPS 10000

#define PARTICLE_TYPE 2

Particle generateOneParticle(void);

void generateTargetCoordinates(double* x, double* y, double* z);
int moveParticleToTarget(Particle* particle, double targetX, double targetY, double targetZ, double stepSize);

#endif