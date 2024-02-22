#pragma once

class Strand;

#include <vector>
#include "domain.h"

typedef struct atom_ {
    Domain         domain;
    struct atom_*  partner;
    Strand*        strand;
} Atom;

class Strand {
    public:
        Strand();
        ~Strand();

        Strand* copy();

        unsigned addDomain(Domain domain);
        static void pairDomain(Atom* atom1, Atom* atom2);
        static void halfPairDomain(Atom* atom1, Atom* atom2);
        static void halfUnpairDomain(Atom* atom1, Atom* atom2);
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
