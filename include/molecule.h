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
        void     unPair(unsigned index);
        unsigned size();

    private:
        std::vector<Strand*>* strands;

};