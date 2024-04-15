#pragma once

class Molecule;

#include <vector>
#include "strand.h"

class Molecule {
    public:
        Molecule(Strand *mainStrand);
        ~Molecule();

        /**
         * Bind from base to base
         * @param strand strand to bind
         */
        unsigned addStrand(Strand *strand, int bindFrom);
        void     remStrand(unsigned index);
        void     pairUncomplete(unsigned index);
        Strand*  getStrand(unsigned index);

        void     donePairStrands(unsigned index, unsigned from, unsigned to);
        unsigned size();
        bool     finishDelete();

        void     simulate(std::map<DOMAIN_DT, Nucleotides*> &nucleotides, float temp, std::vector<Molecule*> &instruction, unsigned strands_count, unsigned time);
        void     deterministicBind();

    private:
        std::vector<Strand*>* strands;
        std::vector<std::vector<NUCLEOTIDE_DT>>* nucleotides;
        std::vector<std::vector<int>>*           nucleotidesBindings;
        std::vector<Strand*>*                    nucleotides2strand;
        std::vector<std::vector<unsigned>>*      nucleotides2domain;

};
