#include "../include/spatialDomain.h"

// Initialize spatial domain
int initializeDomain(SpatialDomain* domain, int particleCount)
{
    domain->particleCount = particleCount;

    initializeCellLinkedList(domain);

    return 1;
}

// Initialize Cell Linked List
void initializeCellLinkedList(SpatialDomain* domain)
{
    // -1 means that the cell is empty
    for (int iX = 0; iX < NCELLS; iX++) {
        for (int jY = 0; jY < NCELLS; jY++) {
            for (int kZ = 0; kZ < NCELLS; kZ++) {
                domain->head[iX][jY][kZ] = -1;
            }
        }
    }
}

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

// Add the particle to the list for its cell
void insertParticle(SpatialDomain* domain, int particleIndex)
{
    Particle* candidateParticle = &domain->particles[particleIndex];

    int candidateIX = getCellIndex(candidateParticle->x);
    int candidateJY = getCellIndex(candidateParticle->y);
    int candidateKZ = getCellIndex(candidateParticle->z);

    domain->next[particleIndex] = domain->head[candidateIX][candidateJY][candidateKZ];
    domain->head[candidateIX][candidateJY][candidateKZ] = particleIndex;
}