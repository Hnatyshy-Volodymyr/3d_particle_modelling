#include <stdlib.h>
#include <math.h>

#include "../include/oneParticleConfiguration.h"

static double randomCoordinate(void)
{
    return L * rand() / ((double)RAND_MAX + 1.0);
}

Particle generateOneParticle(void)
{
    Particle particle;

    particle.x = randomCoordinate();
    particle.y = randomCoordinate();
    particle.z = randomCoordinate();

    particle.type = PARTICLE_TYPE;

    if (particle.type == 1) {
        particle.radius = RADIUS_1;
    }
    else {
        particle.radius = RADIUS_2;
    }

    return particle;
}

void generateTargetCoordinates(double* targetX, double* targetY, double* targetZ)
{
    *targetX = randomCoordinate();
    *targetY = randomCoordinate();
    *targetZ = randomCoordinate();
}

// Move particle towards target
// "Move particle... MOVE! *A lash of the whip* MOVE"

int moveParticleToTarget(Particle* particle, double targetX, double targetY, double targetZ, double stepSize)
{
    double dx = targetX - particle->x;
    double dy = targetY - particle->y;
    double dz = targetZ - particle->z;

    double distance = sqrt(dx * dx + dy * dy + dz * dz);

    if (distance <= stepSize) {
        particle->x = targetX;
        particle->y = targetY;
        particle->z = targetZ;

        return 1;
    }

    particle->x += dx / distance * stepSize;
    particle->y += dy / distance * stepSize;
    particle->z += dz / distance * stepSize;

    return 0;
}