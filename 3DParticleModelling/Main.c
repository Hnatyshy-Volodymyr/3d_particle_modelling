// H.V.V.
/* Cell Linked List has been used to reduce the number of checks and optimise particle generation. */

#define _CRT_SECURE_NO_WARNINGS

//#include <windows.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define type_1Quantity 50
#define type_2Quantity 450

#define N (type_1Quantity + type_2Quantity)

#define L 10.0

#define DIAMETER_1 2.0
#define DIAMETER_2 1.0

#define RADIUS_1 (DIAMETER_1 / 2.0)
#define RADIUS_2 (DIAMETER_2 / 2.0)

#define MAX_RADIUS RADIUS_1

#define NCELLS 10
#define MAX_ATTEMPTS 10000000

//PBC: 1 - axis on, 0 - axis off.
#define PBC_FLAG_X 1
#define PBC_FLAG_Y 1
#define PBC_FLAG_Z 1

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
and retrieve the cells and the particles within them from there
If PBC does not work for axis, we do not transpose the cells.
*/
int translateCellParticle(int index, int pbcEnabled)
{
    if (pbcEnabled == 0) {
        if (index < 0 || index >= NCELLS)
            return -1;

        return index;
    }

    index = index % NCELLS;

    if (index < 0) { index += NCELLS; }

    if (index >= NCELLS) { index -= NCELLS; }

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
int canPlaceParticle(
    const Particle particles[],
    const Particle* candidateParticle,
    const int head[NCELLS][NCELLS][NCELLS],
    const int next[N],
    int* pbcUsedCounter)
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

                for (int i = head[originalNeighborCellIX][originalNeighborCellJY][originalNeighborCellKZ]; i != -1; i = next[i]) {

                    Particle periodicParticleImage = particles[i];

                    bool pbcUsedFlag = (neighborCellIX != originalNeighborCellIX || neighborCellJY != originalNeighborCellJY || neighborCellKZ != originalNeighborCellKZ);

                    //int pbc_used_count = 0;

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
                            } else if (*pbcUsedCounter % 100000 == 0 && *pbcUsedCounter <= 500000) {
                                printf("PBC overlap rejection count>=: %d\n", *pbcUsedCounter);
                            } else if (*pbcUsedCounter % 1000000 == 0 && *pbcUsedCounter <= 5000000) {
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

// Add the particle to the list for its cell
void insertParticle(
    int particleIndex,
    const Particle* candidateParticle,
    int head[NCELLS][NCELLS][NCELLS],
    int next[N])
{
    int candidateIX = getCellIndex(candidateParticle->x);
    int candidateJY = getCellIndex(candidateParticle->y);
    int candidateKZ = getCellIndex(candidateParticle->z);

    next[particleIndex] = head[candidateIX][candidateJY][candidateKZ];
    head[candidateIX][candidateJY][candidateKZ] = particleIndex;
}

// Generate a configuration without overlapping particles
int generateConfiguration(
    Particle particles[],
    int head[NCELLS][NCELLS][NCELLS],
    int next[N],
    int* attempts,
    int* pbcUsedCounter)
{
    int count = 0;
    //Number of particles still to create
    //int typeN1Remaining = N1;
    //int typeN2Remaining = N2;

    *attempts = 0;

    while (count < N && *attempts < MAX_ATTEMPTS) {
        int type;

        type = (count < type_1Quantity) ? 1 : 2;

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

        Particle candidateParticle = generateParticle(type);
        (*attempts)++;

        if (!canPlaceParticle(particles, &candidateParticle, head, next, pbcUsedCounter)) {
            continue;
        }

        particles[count] = candidateParticle;
        insertParticle(count, &candidateParticle, head, next);

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

    fprintf(file,"Properties=species:S:1:pos:R:3:Radius:R:1\n");

    for (int i = 0; i < count; i++) {

        const char* species = (particles[i].type == 1) ? "N" : "O";

        //Centred coordinate for .XYZ file.
        //- L/2.0 - Half of box

        fprintf(file, "%s %.6f %.6f %.6f %.6f\n",
            species,
            particles[i].x - L / 2.0,
            particles[i].y - L / 2.0,
            particles[i].z - L / 2.0,
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
    int pbcRejectCounter = 0;

    // -1 means that the cell is empty
    for (int iX = 0; iX < NCELLS; iX++) {
        for (int jY = 0; jY < NCELLS; jY++) {
            for (int kZ = 0; kZ < NCELLS; kZ++) {
                head[iX][jY][kZ] = -1;
            }
        }
    }

    printf("Hard-sphere configuration generator\n");
    printf("Target particle count: %d\n", N);
    printf("Box dimensions: %.1f x %.1f x %.1f Angstroms\n", L, L, L);

    printf("Type 1: %d particles, diameter %.1f Angstroms\n", type_1Quantity, DIAMETER_1);
    printf("Type 2: %d particles, diameter %.1f Angstroms\n\n", type_2Quantity, DIAMETER_2);

    int count = generateConfiguration(particles, head, next, &attempts, &pbcRejectCounter);

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
    printf("Type 1 particles: %d / %d\n", type1Count, type_1Quantity);
    printf("Type 2 particles: %d / %d\n", type2Count, type_2Quantity);
    printf("Number of attempts: %d\n", attempts);
    printf("Spatial cells: %d\n", NCELLS * NCELLS * NCELLS);
    //printf("Output file: particles.pdb\n");
    printf("Output file: particles.xyz\n");

    return 0;
}