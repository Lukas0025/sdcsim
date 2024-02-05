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

        unsigned addDomain(Domain domain);
        static void pairDomain(Atom* atom1, Atom* atom2);
        void alignmentPrint();
        Atom* getAtom(unsigned index);
        unsigned length();

    private:
        Atom createAtom(Domain domain);

        std::vector<Atom>* atoms;
        bool readOnly;

};
