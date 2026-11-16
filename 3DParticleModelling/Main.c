// H.V.V.
/* Cell Linked List has been used to reduce the number of checks and optimise particle generation. */

#define _CRT_SECURE_NO_WARNINGS

//#include <windows.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define N1 50
#define N2 450

#define N (N1 + N2)

#define L 10.0

#define DIAMETER_1 2.0
#define DIAMETER_2 1.0

#define RADIUS_1 (DIAMETER_1 / 2.0)
#define RADIUS_2 (DIAMETER_2 / 2.0)

#define MAX_RADIUS RADIUS_1

#define NCELLS 10
#define MAX_ATTEMPTS 10000000

typedef struct {
    double x, y, z;
    double radius;
    int type;
} Particle;

// Find the cell for a given coordinate
int getCellIndex(double coordinate)
{
    int index = (int)(coordinate / (L / NCELLS));

    if (index < 0)
        index = 0;
    if (index >= NCELLS)
        index = NCELLS - 1;

    return index;
}

/* If we need to access a cell outside the box, we jump to the opposite side of the box 
and retrieve the cells and the particles within them from there*/
int translateCellParticle(int index)
{
    if (index < 0)
        return index + NCELLS;

    if (index >= NCELLS)
        return index - NCELLS;

    return index;
}

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
int particlesOverlap(const Particle* a, const Particle* b)
{
    double dx = a->x - b->x;
    double dy = a->y - b->y;
    double dz = a->z - b->z;

    double r2 = dx * dx + dy * dy + dz * dz;

    double minDistance = a->radius + b->radius;

    return r2 < minDistance * minDistance;
}

// Check particles in the same and nearby cells
int canPlaceParticle(
    const Particle particles[],
    const Particle* candidate,
    const int head[NCELLS][NCELLS][NCELLS],
    const int next[N],
    int* pbc_used_counter)
{
    int ix = getCellIndex(candidate->x);
    int jy = getCellIndex(candidate->y);
    int kz = getCellIndex(candidate->z);

    //ceil() - rounding number up
    int overlapCheckRange = (int)ceil((candidate->radius + MAX_RADIUS) / (L / NCELLS));

    for (int x = ix - overlapCheckRange; x <= ix + overlapCheckRange; x++) {
        for (int y = jy - overlapCheckRange; y <= jy + overlapCheckRange; y++) {
            for (int z = kz - overlapCheckRange; z <= kz + overlapCheckRange; z++) {

                //bx - tranclate x cordinate
                int bx = translateCellParticle(x);
                int by = translateCellParticle(y);
                int bz = translateCellParticle(z);

                // Each cell has a list of particle indices
                for (int i = head[bx][by][bz]; i != -1; i = next[i]) {

                    Particle image = particles[i];

                    bool pbcUsedFlag = (x != bx || y != by || z != bz);

                    //int pbc_used_count = 0;

                    // Creating particle image
                    image.x += (x - bx) * (L / NCELLS);
                    image.y += (y - by) * (L / NCELLS);
                    image.z += (z - bz) * (L / NCELLS);

                    if (particlesOverlap(candidate, &image)) {
                        if (pbcUsedFlag) {
                            (*pbc_used_counter)++;

                            /*Counter for the number of particle deviations 
                            caused by the periodic application of boundary conditions*/
                            if (*pbc_used_counter % 10000 == 0) {
                                printf("PBC overlap rejection count: %d\n", *pbc_used_counter);
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

// Add the particle to the list for its cell
void insertParticle(
    int particle_index,
    const Particle* particle,
    int head[NCELLS][NCELLS][NCELLS],
    int next[N])
{
    int cx = getCellIndex(particle->x);
    int cy = getCellIndex(particle->y);
    int cz = getCellIndex(particle->z);

    next[particle_index] = head[cx][cy][cz];
    head[cx][cy][cz] = particle_index;
}

// Generate a configuration without overlapping particles
int generate_configuration(
    Particle particles[],
    int head[NCELLS][NCELLS][NCELLS],
    int next[N],
    int* attempts,
    int* pbc_used_counter)
{
    int count = 0;
    //Number of particles still to create
    //int typeN1Remaining = N1;
    //int typeN2Remaining = N2;

    *attempts = 0;

    while (count < N && *attempts < MAX_ATTEMPTS) {
        int type;

        type = (count < N1) ? 1 : 2;

        /*Randomly selects which type to generate.
        If one type of particle has already been generated, we move on to the next.

        P.S: I came up with the task myself. 
        I created the problem myself. 
        *AAAAAAAAH!* 
        I solved the problem myself...”
        */

        /*
        if (typeN1Remaining == 0) {
            type = 2;
        }
        else if (typeN2Remaining == 0) {
            type = 1;
        }
        else {
            type = (rand() % (typeN1Remaining + typeN2Remaining) < typeN1Remaining) ? 1 : 2;
        }
        */

        Particle candidate = generateParticle(type);
        (*attempts)++;
        
        if (!canPlaceParticle(particles, &candidate, head, next, pbc_used_counter)) {
            continue;
        }

        particles[count] = candidate;
        insertParticle(count, &candidate, head, next);

        /*
        if (type == 1) {
            typeN1Remaining--;
        }
        else {
            typeN2Remaining--;
        }
        */

        count++;


        if (count % 100 == 0 || count == N) {
            printf("Placed particles: %d / %d\n", count, N);
        }
    }

    return count;
}

/*
int save_pdb(const Particle particles[], int count)
{
    FILE* file = fopen("particles.pdb", "w");

    if (file == NULL) {
        perror("Error opening particles.pdb");
        return 0;
    }

    // Box dimensions and angles
    fprintf(file,
        "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
        L, L, L, 90.0, 90.0, 90.0);

    for (int i = 0; i < count; i++) {
        fprintf(file,
            "ATOM  %5d  C   UNK A%4d    "
            "%8.3f%8.3f%8.3f%6.2f%6.2f          C  \n",
            i + 1,
            i + 1,
            particles[i].x,
            particles[i].y,
            particles[i].z,
            1.00,
            0.00);
    }

    fprintf(file, "TER\n");
    fprintf(file, "END\n");

    if (fclose(file) != 0) {
        perror("Error closing particles.pdb");
        return 0;
    }

    return 1;
}
*/

// Save particle configuration in Extended XYZ format for OVITO
int saveXYZ(const Particle particles[], int count)
{
    FILE* file = fopen("particles.xyz", "w");

    if (file == NULL) {
        perror("Error opening particles.xyz");
        return 0;
    }

    fprintf(file, "%d\n", count);

    fprintf(file,
        "Lattice=\"%.6f 0 0 0 %.6f 0 0 0 %.6f\" "
        "Properties=species:S:1:pos:R:3:Radius:R:1 "
        "pbc=\"T T T\"\n",
        L, L, L);

    for (int i = 0; i < count; i++) {

        const char* species = (particles[i].type == 1) ? "N" : "O";

        fprintf(file, "%s %.6f %.6f %.6f %.6f\n",
            species,
            particles[i].x,
            particles[i].y,
            particles[i].z,
            particles[i].radius);
    }

    if (fclose(file) != 0) {
        perror("Error closing particles.xyz");
        return 0;
    }

    return 1;
}

int main(void)
{
    srand((unsigned int)time(NULL));

    Particle particles[N];
    int head[NCELLS][NCELLS][NCELLS];
    int next[N];
    int attempts = 0;
    int pbc_used_counter = 0;

    // -1 means that the cell is empty
    for (int x = 0; x < NCELLS; x++) {
        for (int y = 0; y < NCELLS; y++) {
            for (int z = 0; z < NCELLS; z++) {
                head[x][y][z] = -1;
            }
        }
    }

    printf("Hard-sphere configuration generator\n");
    printf("Target particle count: %d\n", N);
    printf("Box dimensions: %.1f x %.1f x %.1f Angstroms\n",L, L, L);

    printf("Type 1: %d particles, diameter %.1f Angstroms\n",N1, DIAMETER_1);
    printf("Type 2: %d particles, diameter %.1f Angstroms\n\n",N2, DIAMETER_2);

    int count = generate_configuration(
        particles, head, next, &attempts, &pbc_used_counter);

    int type1Count = 0;
    int type2Count = 0;

    for (int i = 0; i < count; i++) {
        if (particles[i].type == 1) {
            type1Count++;
        }
        else if (particles[i].type == 2) {
            type2Count++;
        }
    }

    if (count != N) {
        printf("\nWARNING: Could place only %d of %d particles.\n", count, N);
        printf("Maximum number of attempts reached: %d\n", attempts);
        printf("Saving the generated partial configuration.\n");
    }

    /*
    if (!save_pdb(particles, count)) {
        return 1;
    }
    */

    if (!saveXYZ(particles, count)) {
        return 1;
    }

    if (count == N) {
        printf("\nGeneration completed successfully.\n");
    }
    else {
        printf("\nGeneration stopped. Partial configuration saved.\n");
    }

    printf("Particles placed: %d / %d\n", count, N);
    printf("Type 1 particles: %d / %d\n", type1Count, N1);
    printf("Type 2 particles: %d / %d\n", type2Count, N2);
    printf("Number of attempts: %d\n", attempts);
    printf("Spatial cells: %d\n", NCELLS * NCELLS * NCELLS);
    //printf("Output file: particles.pdb\n");
    printf("Output file: particles.xyz\n");

    return 0;
}