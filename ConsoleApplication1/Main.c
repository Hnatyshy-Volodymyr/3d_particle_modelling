// H.V.V.
/* Cell Linked List has been used to reduce the number of checks and optimise particle generation. */

#define _CRT_SECURE_NO_WARNINGS

#include <windows.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define N 512
#define L 10.0
#define SIGMA 1.0
#define NCELLS 10
#define MAX_ATTEMPTS 100000

typedef struct {
    double x, y, z;
} Particle;

// Find the cell for a given coordinate
int get_cell_index(double coordinate)
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
int translate_cell_particle(int index)
{
    if (index < 0)
        return index + NCELLS;

    if (index >= NCELLS)
        return index - NCELLS;

    return index;
}

// Random position inside the box
double random_coordinate(void)
{
    return L * rand() / ((double)RAND_MAX + 1.0);
}

Particle generate_particle(void)
{
    Particle p;

    p.x = random_coordinate();
    p.y = random_coordinate();
    p.z = random_coordinate();

    return p;
}

// Check the distance between particle centers with using periodic boundary conditions
int particles_overlap(const Particle* a, const Particle* b)
{
    double dx = a->x - b->x;
    double dy = a->y - b->y;
    double dz = a->z - b->z;

    double r2 = dx * dx + dy * dy + dz * dz;

    return r2 < SIGMA * SIGMA;
}

// Check particles in the same and nearby cells
int can_place_particle(
    const Particle particles[],
    const Particle* candidate,
    const int head[NCELLS][NCELLS][NCELLS],
    const int next[N],
    int count,
    int* pbc_used_counter)
{
    int ix = get_cell_index(candidate->x);
    int jy = get_cell_index(candidate->y);
    int kz = get_cell_index(candidate->z);

    for (int x = ix - 1; x <= ix + 1; x++) {
        for (int y = jy - 1; y <= jy + 1; y++) {
            for (int z = kz - 1; z <= kz + 1; z++) {

                //bx - broadcast x cordinate
                int bx = translate_cell_particle(x);
                int by = translate_cell_particle(y);
                int bz = translate_cell_particle(z);

                // Each cell has a list of particle indices
                for (int i = head[bx][by][bz]; i != -1; i = next[i]) {

                    Particle image = particles[i];

                    bool pbc_used_flag = (x != bx || y != by || z != bz);

                    //int pbc_used_count = 0;

                    // Creating particle image
                    image.x += (x - bx) * (L / NCELLS);
                    image.y += (y - by) * (L / NCELLS);
                    image.z += (z - bz) * (L / NCELLS);

                    if (particles_overlap(candidate, &image)) {
                        if (pbc_used_flag) {
                            (*pbc_used_counter)++;

                            /*Counter for the number of particle deviations 
                            caused by the periodic application of boundary conditions*/
                            if (*pbc_used_counter % 100 == 0) {
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
void insert_particle(
    int particle_index,
    const Particle* particle,
    int head[NCELLS][NCELLS][NCELLS],
    int next[N])
{
    int cx = get_cell_index(particle->x);
    int cy = get_cell_index(particle->y);
    int cz = get_cell_index(particle->z);

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
    *attempts = 0;

    while (count < N && *attempts < MAX_ATTEMPTS) {
        Particle candidate = generate_particle();
        (*attempts)++;

        if (!can_place_particle(particles, &candidate, head, next,
            count, pbc_used_counter)) {
            continue;
        }

        particles[count] = candidate;
        insert_particle(count, &candidate, head, next);
        count++;

        if (count % 100 == 0) {
            printf("Placed particles: %d / %d\n", count, N);
        } else if(count % N == 0){
            printf("Placed particles: %d / %d\n", count, N);
        }
    }

    return count;
}

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
    printf("Box dimensions: %.1f x %.1f x %.1f Angstroms\n",
        L, L, L);
    printf("Particle diameter: %.1f Angstroms\n\n", SIGMA);

    int count = generate_configuration(
        particles, head, next, &attempts, &pbc_used_counter);

    if (count != N) {
        printf("Could place only %d of %d particles.\n", count, N);
        printf("Number of attempts: %d\n", attempts);
        printf("Try increasing the box size or the attempt limit.\n");
        return 1;
    }

    if (!save_pdb(particles, count)) {
        return 1;
    }

    printf("\nGeneration completed successfully.\n");
    printf("Particles placed: %d\n", count);
    printf("Number of attempts: %d\n", attempts);
    printf("Spatial cells: %d\n", NCELLS * NCELLS * NCELLS);
    printf("Output file: particles.pdb\n");

    return 0;
}