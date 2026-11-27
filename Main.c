// H.V.V.
/* Cell Linked List has been used to reduce the number of checks and optimise particle generation. */
//v0.10-24.11

#define _CRT_SECURE_NO_WARNINGS

#define CONFIGURATION_MODE 2

#define HARD_SPHERE_CONFIGURATION 1
#define ONE_PARTICLE_CONFIGURATION 2

//#include <windows.h>
#include <stdio.h>
#include <stdlib.h>

#include "include/particle.h"
#include "include/spatialDomain.h"
#include "include/randomConfiguration.h"
#include "include/output.h"
#include "include/oneParticleConfiguration.h"

void runHardSphereConfiguration(void)
{
    SpatialDomain domain;
    int attempts = 0;
    int pbcRejectCounter = 0;

    initializeDomain(&domain, N);

    printf("Hard-sphere configuration generator\n");
    printf("Target particle count: %d\n", N);
    printf("Box dimensions: %.1f x %.1f x %.1f Angstroms\n", L, L, L);

    printf("Type 1: %d particles, diameter %.1f Angstroms\n", type_1Number, DIAMETER_1);
    printf("Type 2: %d particles, diameter %.1f Angstroms\n\n", type_2Number, DIAMETER_2);

    int count = generateConfiguration(&domain, &attempts, &pbcRejectCounter);

    int type1Count = 0;
    int type2Count = 0;

    for (int i = 0; i < count; i++) {
        if (domain.particles[i].type == 1) {
            type1Count++;
        }
        else if (domain.particles[i].type == 2) {
            type2Count++;
        }
    }

    if (count != N) {
        printf("\nWARNING: Could place only %d of %d particles.\n", count, N);
        printf("Maximum number of attempts reached: %d\n", attempts);
        printf("Saving the generated partial configuration.\n");
    }

    if (!saveConfiguration(domain.particles, count)) {
        return;
    }

    if (count == N) {
        printf("\nGeneration completed successfully.\n");
    }
    else {
        printf("\nGeneration stopped. Partial configuration saved.\n");
    }

    printf("Particles placed: %d / %d\n", count, N);
    printf("Type 1 particles: %d / %d\n", type1Count, type_1Number);
    printf("Type 2 particles: %d / %d\n", type2Count, type_2Number);
    printf("Number of attempts: %d\n", attempts);
    printf("Spatial cells: %d\n", NCELLS * NCELLS * NCELLS);

    if (OUTPUT_FORMAT == OUTPUT_FORMAT_XYZ) {
        printf("Output file: %sparticles.xyz\n", OUTPUT_DIRECTORY);
    }
    else if (OUTPUT_FORMAT == OUTPUT_FORMAT_PDB) {
        printf("Output file: %sparticles.pdb\n", OUTPUT_DIRECTORY);
    }
}

void runOneParticleSimulation(void)
{
    /*We do not run the Domain initialization 
    because it includes the initialization of the Cell - Linked List, 
    which we do not need for a single particle.*/

    SpatialDomain domain;

    domain.particleCount = 1;

    domain.particles[0] = generateOneParticle();

    Particle* particle = &domain.particles[0];

    printf("One particle simulation\n");
    printf("Box dimensions: %.1f x %.1f x %.1f Angstroms\n", L, L, L);
    printf("Particle type: %d\n", particle->type);
    printf("Particle radius: %.2f Angstroms\n", particle->radius);
    printf("Initial position: (%.4f, %.4f, %.4f)\n",
        particle->x, particle->y, particle->z);

    if (!saveConfiguration(domain.particles, domain.particleCount)) {
        printf("Failed to save configuration.\n");
        return;
    }

    printf("Configuration saved successfully.\n");
}

int main(void)
{
    srand(SEED);

    if (CONFIGURATION_MODE == HARD_SPHERE_CONFIGURATION) {
        runHardSphereConfiguration();
    }
    else if (CONFIGURATION_MODE == ONE_PARTICLE_CONFIGURATION) {
        runOneParticleSimulation();
    }

    return 0;
}