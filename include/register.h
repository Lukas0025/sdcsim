#pragma once

#include "molecule.h"

class Register {
    public:
        Register(Molecule *init);

        Molecule* get();

        void applyInstruction(std::vector<Molecule*> instruction);
        void enableNucleotidesLevel(std::map<DOMAIN_DT, Nucleotides*> &nucleotides, float temp);

    private:
        void doAllBinding(Molecule* mol);
        void removeReplaced(unsigned oldSize);
        void removeUnstable(int to);
        void removeUnbinded(Molecule* mol);
        void unbindOfMultiple();
        void bindOfMultiple();

        Molecule *reg; // molecule reprezenting register

        bool nucleotidesSim;
        float temp;
        std::map<DOMAIN_DT, Nucleotides*>* nucleotides;
};