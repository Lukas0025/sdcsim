class Molecule;

#include <vector>
#include "strand.h"

class Molecule {
    public:
        Molecule();
        ~Molecule();

        unsigned addStrand(Strand strand);
        Strand   remStrand(unsigned index);

        void     pairStrands();

    private:


};
