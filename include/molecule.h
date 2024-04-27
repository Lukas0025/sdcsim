/**
 * SDCSIM - strand displacement simulator on SIMD||DNA architecture
 * @file molecule.h
 * @brief Contain headers for molecule class
 * @author Lukáš Plevač <xpleva07@vut.cz>
 */

#pragma once

class Molecule;

#include <vector>
#include <map>
#include "gpu.h"
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

        /**
         * Init GPU for simulation
         * @param nucleotides nucleotides for simulation
         */
        void initGpu(std::map<DOMAIN_DT, Nucleotides*> &nucleotides);

        /**
          * Free gpu resources (deinit)
          */
        void freeGpu();

        /**
          * Copy data from GPU to HOST
          */
        void updateSelf();

        /**
          * Copy data from HOST to GPU
          */
        void updateGpu();

        /**
         * Do strochastic simulation in molecule on GPU acelerator
         * @param nucleotides    map of nucleotides reprezentation
         * @param temp           temperature of solution with molecule
         * @param instruction    vector of incoming molecules / strands
         * @param strands_count  quntity in incoming molecules / strads
         * @param time           time of simulation
         */
        void simulateGpu(std::map<DOMAIN_DT, Nucleotides*> &nucleotides, float temp, std::vector<Molecule*> &instruction, unsigned strands_count, unsigned sim_time);

    private:
        std::vector<Strand*>* strands; //< @brief strads in molecule

        /* GPU */
        int* Gmem  = NULL; // 3D ARRAY 
        bool onGPU = false;
        
};
