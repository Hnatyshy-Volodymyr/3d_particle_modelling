#ifndef RANDOM_CONFIGURATION_H
#define RANDOM_CONFIGURATION_H

#include "spatialDomain.h"

// M-1 Random Config

double randomCoordinate(void);
Particle generateParticle(int type);

int particlesOverlap(const Particle* particleA, const Particle* particleB);

int canPlaceParticle(SpatialDomain* domain, Particle* candidateParticle, int* pbcUsedCounter);

int generateConfiguration(SpatialDomain* domain, int* attempts, int* pbcUsedCounter);

#endif