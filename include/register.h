#pragma once

#include "molecule.h"

class Register {
    public:
        Register(Molecule *init);

        Molecule* get();

        void applyInstruction(std::vector<Molecule*> instruction);

    private:
        void doAllBinding(Molecule* mol);
        void removeReplaced(unsigned oldSize);
        void removeUnstable();
        void removeUnbinded(Molecule* mol);

        Molecule *reg; // molecule reprezenting register
};