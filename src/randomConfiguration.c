#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../include/randomConfiguration.h"

// Random position inside the box
double randomCoordinate(void)
{
    return L * rand() / ((double)RAND_MAX + 1.0);
}

Particle generateParticle(int type)
{
    Particle p;

    p.x = randomCoordinate();
    p.y = randomCoordinate();
    p.z = randomCoordinate();

    p.type = type;

    if (type == 1) {
        p.radius = RADIUS_1;
    }
    else {
        p.radius = RADIUS_2;
    }

    return p;
}

// Check the distance between particle centers with using periodic boundary conditions
int particlesOverlap(const Particle* particleA, const Particle* particleB)
{
    double dx = particleA->x - particleB->x;
    double dy = particleA->y - particleB->y;
    double dz = particleA->z - particleB->z;

    double r2 = dx * dx + dy * dy + dz * dz;

    double minDistance = particleA->radius + particleB->radius;

    return r2 < minDistance * minDistance;
}

// Check particles in the same and nearby cells
int canPlaceParticle(SpatialDomain* domain, Particle* candidateParticle, int* pbcUsedCounter)
{
    int candidateIX = getCellIndex(candidateParticle->x);
    int candidateJY = getCellIndex(candidateParticle->y);
    int candidateKZ = getCellIndex(candidateParticle->z);

    //ceil() - rounding number up
    int overlapCheckRange = (int)ceil((candidateParticle->radius + MAX_RADIUS) / (L / NCELLS));


    //NeightborCell for candidateCell.
    for (int neighborCellIX = candidateIX - overlapCheckRange; neighborCellIX <= candidateIX + overlapCheckRange; neighborCellIX++) {
        for (int neighborCellJY = candidateJY - overlapCheckRange; neighborCellJY <= candidateJY + overlapCheckRange; neighborCellJY++) {
            for (int neighborCellKZ = candidateKZ - overlapCheckRange; neighborCellKZ <= candidateKZ + overlapCheckRange; neighborCellKZ++) {

                int originalNeighborCellIX = translateCellParticle(neighborCellIX, PBC_FLAG_X);
                int originalNeighborCellJY = translateCellParticle(neighborCellJY, PBC_FLAG_Y);
                int originalNeighborCellKZ = translateCellParticle(neighborCellKZ, PBC_FLAG_Z);

                if (originalNeighborCellIX == -1 || originalNeighborCellJY == -1 || originalNeighborCellKZ == -1) {
                    continue;
                }

                for (int i = domain->head[originalNeighborCellIX][originalNeighborCellJY][originalNeighborCellKZ]; i != -1; i = domain->next[i]) {

                    Particle periodicParticleImage = domain->particles[i];

                    bool pbcUsedFlag = (neighborCellIX != originalNeighborCellIX || neighborCellJY != originalNeighborCellJY || neighborCellKZ != originalNeighborCellKZ);

                    // Creating particle image
                    periodicParticleImage.x += (neighborCellIX - originalNeighborCellIX) * (L / NCELLS);
                    periodicParticleImage.y += (neighborCellJY - originalNeighborCellJY) * (L / NCELLS);
                    periodicParticleImage.z += (neighborCellKZ - originalNeighborCellKZ) * (L / NCELLS);

                    if (particlesOverlap(candidateParticle, &periodicParticleImage)) {
                        if (pbcUsedFlag) {
                            (*pbcUsedCounter)++;

                            /*Counter for the number of particle deviations
                            caused by the periodic application of boundary conditions*/
                            if (*pbcUsedCounter % 10000 == 0 && *pbcUsedCounter <= 50000) {
                                printf("PBC overlap rejection count>=: %d\n", *pbcUsedCounter);
                            }
                            else if (*pbcUsedCounter % 100000 == 0 && *pbcUsedCounter <= 500000) {
                                printf("PBC overlap rejection count>=: %d\n", *pbcUsedCounter);
                            }
                            else if (*pbcUsedCounter % 1000000 == 0 && *pbcUsedCounter <= 5000000) {
                                printf("PBC overlap rejection count>=: %d\n", *pbcUsedCounter);
                            }
                        }

                        return 0;
                    }
                }
            }
        }
    }

    return 1;
}

// Generate a configuration without overlapping particles
int generateConfiguration(SpatialDomain* domain, int* attempts, int* pbcUsedCounter)
{
    int count = 0;

    /*int type_1Remaining = type_1Number;
    int type_2Remaining = type_2Number;*/

    *attempts = 0;

    while (count < domain->particleCount && *attempts < MAX_ATTEMPTS) {
        int type;

        type = (count < type_1Number) ? 1 : 2;

        /*Randomly selects which type to generate.
        If one type of particle has already been generated, we move on to the next.

        P.S: I came up with the task myself.
        I created the problem myself.
        *AAAAAAAAH!*
        I solved the problem myself...”
        */

        /*if (type_1Remaining > 0 && type_2Remaining > 0) {
            type = rand() % 2 + 1;
        }
        else if (type_1Remaining > 0) {
            type = 1;
        }
        else {
            type = 2;
        }*/

        Particle candidateParticle = generateParticle(type);
        (*attempts)++;

        if (!canPlaceParticle(domain, &candidateParticle, pbcUsedCounter)) {
            continue;
        }

        domain->particles[count] = candidateParticle;
        insertParticle(domain, count);

        /*if (type == 1) {
            type_1Remaining--;
        }
        else {
            type_2Remaining--;
        }*/

        count++;

        if (count % 100 == 0 || count == domain->particleCount) {
            printf("Placed particles: %d / %d\n", count, domain->particleCount);
        }
    }

    return count;
}