/**
 * SDCSIM - strand displacement simulator on SIMD||DNA architecture
 * @file molecule.h
 * @brief Contain headers for molecule class
 * @author Lukáš Plevač <xpleva07@vut.cz>
 */

#pragma once

class Molecule;

#define MAX_GPU_STRANDS    1000
#define MAX_GPU_STRAND_LEN 255
#define GPU_DEPTH          4
#define GPU_BACK_DEPTH     2

#define GPU_MEM_SIZE       (MAX_GPU_STRANDS * MAX_GPU_STRAND_LEN * GPU_DEPTH + 1)
#define GPU_BACK_MEM_SIZE  (MAX_GPU_STRANDS * MAX_GPU_STRAND_LEN * GPU_BACK_DEPTH + 1)

// definition of depth
#define GPU_PARTNER_STRAND 0
#define GPU_PARTNER_ATOM   1 
#define GPU_DOMAIN         2
#define GPU_NUCLEOTIDES    3

// on first position is alwais length
#define GPU_ELEMENT(y, x, z)    (2 + x + y * MAX_GPU_STRANDS + z * MAX_GPU_STRAND_LEN * MAX_GPU_STRANDS)
#define GET_GPU_COUNT_COL(x, z) GPU_ELEMENT(x, -1, z)
#define GET_GPU_COUNT_ROW()     0

#include <vector>
#include "strand.h"

/**
 * @brief CLass reprezenting DNA Molecule in system strand displacement
 */
class Molecule {
    public:
        /**
         * COnstructor of class
         * @param mainStrand bottom strand of molecule in SIMD||DNA architecture
         */
        Molecule(Strand *mainStrand);
        ~Molecule();

        /**
         * Add strand do molecule and bind from base to base
         * @param strand strand to bind
         * @param bindFrom position on masinstrand to bind from
         * @return index of strand in molecule
         */
        unsigned addStrand(Strand *strand, int bindFrom);

        /**
         * Mark strand in molecule as deleted
         * @param index index of strand
         */
        void     remStrand(unsigned index);

        /**
         * pair domian of strand with main what is not binded and can bind
         * @param index index of strand
         */
        void     pairUncomplete(unsigned index);

        /**
         * Get pointer to strand by index
         * @param undex of strand 
         */
        Strand*  getStrand(unsigned index);

        /**
         * Force bind strand with main strand
         * @param index index of strand to bind
         * @param from  index in incoming strand from what we binding
         * @param to    index in incoming strand to what we binding
         */
        void     donePairStrands(unsigned index, unsigned from, unsigned to);

        /**
         * Get Number of strands in molecule
         * @return unsigned count
         */
        unsigned size();

        /**
         * Delete marked strands 
         */
        bool     finishDelete();

        /**
         * Do strochastic simulation in molecule
         * @param nucleotides    map of nucleotides reprezentation
         * @param temp           temperature of solution with molecule
         * @param instruction    vector of incoming molecules / strands
         * @param strands_count  quntity in incoming molecules / strads
         * @param time           time of simulation
         */
        void     simulate(std::map<DOMAIN_DT, Nucleotides*> &nucleotides, float temp, std::vector<Molecule*> &instruction, unsigned strands_count, unsigned time);

        /**
         * Deterministickly done binding of all strads in molecule
         */
        void     deterministicBind();

        void initGpu();
        void freeGpu();
        void updateSelf();
        void updateGpu();

    private:
        std::vector<Strand*>* strands; //< @brief strads in molecule

        /* GPU */
        int* Gmem = NULL; // 3D ARRAY 
        
};
