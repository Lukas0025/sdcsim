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

        void     createNucleotidesLevel(std::map<DOMAIN_DT, Nucleotides*> &nucleotides);
        void     addNucleotidesFreeStrand(Strand* strand, unsigned desity, std::map<DOMAIN_DT, Nucleotides*> &nucleotides);
        void     updateDomainLevel();

    private:
        std::vector<Strand*>* strands;
        std::vector<std::vector<NUCLEOTIDE_DT>>* nucleotides;
        std::vector<std::vector<int>>*           nucleotidesBindings;
        std::vector<Strand*>*                    nucleotides2domain;

};
