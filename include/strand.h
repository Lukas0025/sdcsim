#pragma once

class Strand;

#include <vector>
#include "domain.h"
#include <cstddef>

typedef struct atom_ {
    Domain         domain;
    struct atom_*  partner;
    Strand*        strand;
    int            partnersCount;
} Atom;

extern Atom* multiAtom;

class Strand {
    public:
        Strand();
        ~Strand();

        Strand* copy();

        unsigned addDomain(Domain domain);
        static void pairDomain(Atom* atom1, Atom* atom2);
        static void halfPairDomain(Atom* atom1, Atom* atom2);
        static void unpairDomain(Atom* atom1);
        Atom* getAtom(unsigned index);
        unsigned length();
        void complementaryLast();
        void del();
        bool isDeleted();

    private:
        Atom createAtom(Domain domain);

        std::vector<Atom>* atoms;
        
        bool readOnly;
        bool deleted;

};
