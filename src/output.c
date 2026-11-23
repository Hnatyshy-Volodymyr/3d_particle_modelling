#include <stdio.h>

#include "../include/output.h"

//XYZ format for OVITO
int saveXYZ(const Particle particles[], int count)
{
    FILE* file = fopen(OUTPUT_DIRECTORY "particles.xyz", "w");

    if (file == NULL) {
        perror("Error opening particles.xyz");
        return 0;
    }

    fprintf(file, "%d\n", count);

    fprintf(file, "Origin=\"-%.6f -%.6f -%.6f\" "
        "Lattice=\"%.6f 0 0 0 %.6f 0 0 0 %.6f\" "
        "Properties=species:S:1:pos:R:3:Radius:R:1\n",
        L / 2.0, L / 2.0, L / 2.0, L, L, L);

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

//PDB format
int savePDB(const Particle particles[], int count)
{
    FILE* file = fopen(OUTPUT_DIRECTORY "particles.pdb", "w");

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

//Save config with chosen format
int saveConfiguration(const Particle particles[], int count)
{
    if (OUTPUT_FORMAT == OUTPUT_FORMAT_XYZ) {
        return saveXYZ(particles, count);
    }
    else if (OUTPUT_FORMAT == OUTPUT_FORMAT_PDB) {
        return savePDB(particles, count);
    }

    return 0;
}