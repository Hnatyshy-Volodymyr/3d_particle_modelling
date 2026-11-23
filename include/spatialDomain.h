#ifndef SPATIAL_DOMAIN_H
#define SPATIAL_DOMAIN_H

#include "particle.h"

// Domain
typedef struct {
    Particle particles[N];
    int head[NCELLS][NCELLS][NCELLS];
    int next[N];
    int particleCount;
} SpatialDomain;

int initializeDomain(SpatialDomain* domain, int particleCount);
void initializeCellLinkedList(SpatialDomain* domain);
void insertParticle(SpatialDomain* domain, int particleIndex);

int getCellIndex(double coordinate);
int translateCellParticle(int index, int pbcEnabled);

#endif