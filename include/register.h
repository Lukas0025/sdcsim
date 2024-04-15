#pragma once

#include "molecule.h"

class Register {
    public:
        Register(Molecule *init);

        Molecule* get();

        void applyInstruction(std::vector<Molecule*> instruction, int time);
        void enableNucleotidesLevel(std::map<DOMAIN_DT, Nucleotides*> &nucleotides, float temp, int strands);

    private:
        void doAllBinding(Molecule* mol);
        void removeReplaced(unsigned oldSize);
        void removeUnstable(int to);
        void removeUnbinded(Molecule* mol);
        void unbindOfMultiple();
        void bindOfMultiple();
        void elution(std::map<DOMAIN_DT, Nucleotides*> &nucleotides);

        Molecule *reg; // molecule reprezenting register

        bool nucleotidesSim;
        float temp;
        int strands_count;
        std::map<DOMAIN_DT, Nucleotides*>* nucleotides;
};