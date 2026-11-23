#ifndef OUTPUT_H
#define OUTPUT_H

#include "particle.h"
//Output directory
#define OUTPUT_DIRECTORY "C:\\Users\\1\\Desktop\\PhD\\3DParticleModelling\\3DParticleModelling\\output\\"

//Output format: 1 - XYZ, 2 - PDB
#define OUTPUT_FORMAT_XYZ 1
#define OUTPUT_FORMAT_PDB 2

//If we need PBD, write:
//OUTPUT_FORMAT OUTPUT_FORMAT_PDB
//If we need XYZ, write:
//OUTPUT_FORMAT OUTPUT_FORMAT_XYZ

#define OUTPUT_FORMAT OUTPUT_FORMAT_XYZ

int saveXYZ(const Particle particles[], int count);
int savePDB(const Particle particles[], int count);
int saveConfiguration(const Particle particles[], int count);

#endif