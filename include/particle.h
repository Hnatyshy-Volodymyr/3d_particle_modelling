#ifndef PARTICLE_H
#define PARTICLE_H

#define SEED 12345

#define type_1Number 50
#define type_2Number 450

#define N (type_1Number + type_2Number)

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

#endif