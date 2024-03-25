#include <cstddef>
#include "molecule.h"
#include <cstdlib>
#include <iostream>

#define RANDOM_RANGE(TO)       (rand() % (TO))
#define RANDOM_FLOAT()         ((float)rand()/(float)RAND_MAX)

Molecule::Molecule(Strand *strand) {
    this->strands = new std::vector<Strand*>;    
    this->strands->push_back(strand);
}

Molecule::~Molecule() {
    delete this->strands;
}

inline bool isPossibleToBind(Atom* iAtom, Atom* jAtom) {

    bool status = jAtom->partner == NULL && iAtom->partner == NULL && iAtom->strand != jAtom->strand;

    auto strand = iAtom->strand;
    int  diffN  = (int)(iAtom->offsetN) - jAtom->offsetN;

    // now loop over chain and check if other is binded
    for (unsigned i = 0; i < strand->length(); i++) {
        auto atom = strand->getAtom(i);
        
        status = status && (atom->partner == NULL || atom->partner->strand != jAtom->strand || (int)(atom->offsetN - atom->partner->offsetN) == diffN);
    }

    return status;

}

void Molecule::simulate(std::map<DOMAIN_DT, Nucleotides*> &nucleotides, float temp, std::vector<Molecule*> &instruction) {
    // set random seed
    srand((unsigned)time(NULL));

    unsigned time = 100000;

    while (time > 0) {

        //select i and j strands
        auto iStrand = this->getStrand(RANDOM_RANGE(this->size()));

        // allow select from incomming
        auto jSIndex = RANDOM_RANGE(this->size() + instruction.size());

        // use incomming
        if (jSIndex >= this->size()) {
            auto newStrand = instruction.at(jSIndex - this->size())->getStrand(0)->copy();
            jSIndex        = this->size();
            this->strands->push_back(newStrand);
        }

        auto jStrand = this->getStrand(jSIndex);

        //select i and j atom
        auto iAtomId = RANDOM_RANGE(iStrand->length());
        auto iAtom   = iStrand->getAtom(iAtomId);

        auto jAtomId = RANDOM_RANGE(jStrand->length());
        auto jAtom   = jStrand->getAtom(jAtomId);

        if (RANDOM_FLOAT() < 0.5) { // renaturation
            if (nucleotides[iAtom->domain.get()]->ReDeNaturationP(nucleotides[jAtom->domain.get()], temp) > RANDOM_FLOAT()) {
                if (isPossibleToBind(iAtom, jAtom)) {
                    Strand::pairDomain(iAtom, jAtom);

                    //try bind other other parts if possible - this is created by energy of this creastion of pair
                    for (unsigned i = 0; i <= std::min(iAtomId, jAtomId); i++) {
                        iAtom   = iStrand->getAtom(iAtomId - i);
                        jAtom   = jStrand->getAtom(jAtomId - i);

                        if (!(nucleotides[iAtom->domain.get()]->ReDeNaturationP(nucleotides[jAtom->domain.get()], temp) > RANDOM_FLOAT())) break;

                        if (iAtom->partner != NULL) {
                            iAtom->partner->partner = NULL;
                            iAtom->partner          = NULL;
                        }

                        if (jAtom->partner != NULL) {
                            jAtom->partner->partner = NULL;
                            jAtom->partner          = NULL;
                        }

                        Strand::pairDomain(iAtom, jAtom);   
                    }

                    /*for (unsigned i = iAtomId; i < iStrand->length(); i++) {
                        
                    }*/
                }
            }
        } else if (iAtom->partner != NULL && nucleotides[iAtom->domain.get()]->ReDeNaturationP(nucleotides[iAtom->partner->domain.get()], temp) < RANDOM_FLOAT()) {
            iAtom->partner->partner = NULL;
            iAtom->partner          = NULL;
        }

        time--;
    }
}

unsigned Molecule::addStrand(Strand *strand, int bindFrom) {
    unsigned mainBindFrom   = bindFrom;
    unsigned incomeBindFrom = 0;

    if (bindFrom < 0) {
        incomeBindFrom = -bindFrom;
        mainBindFrom   = 0;
    }

    //first add starnd to array
    unsigned id = this->strands->size();
    this->strands->push_back(strand);

    auto mainStrand = this->getStrand(0);

    //ok now do half pairing
    for (unsigned i = mainBindFrom; i < std::min(strand->length() + mainBindFrom, mainStrand->length()); i++) {
        auto bottom = mainStrand->getAtom(i);
        auto top    = strand->getAtom(i - mainBindFrom + incomeBindFrom);

        if (bottom->domain == ~top->domain) {
            //we can pair this
            Strand::halfPairDomain(top, bottom);
        }
    }

    return id;
}

void Molecule::remStrand(unsigned index) {

    if (index == 0) return; // main strand is nundeleteble

    //do unlink
    for (unsigned i = 0; i < this->strands->at(index)->length(); i++) {
        Strand::unpairDomain(this->strands->at(index)->getAtom(i));
    }

    // set deleted flag
    this->strands->at(index)->del();
}

bool Molecule::finishDelete() {
    bool someDelete = false;

    for (int i = this->strands->size() - 1; i >= 0; i--) {
        if (this->strands->at(i)->isDeleted()) {
            //get lastest index
            this->strands->at(i) = this->strands->at(this->strands->size() - 1);

            //and remove last
            this->strands->pop_back();

            someDelete = true;
        }
    }

    return someDelete;
}

Strand* Molecule::getStrand(unsigned index) {
    return this->strands->at(index);
}

void Molecule::donePairStrands(unsigned index, unsigned from, unsigned to) {
    auto strand = this->getStrand(index);

    for (unsigned i = from; i <= to; i++) {
        auto top    = strand->getAtom(i);
        
        if (top->partner != NULL) {
            Strand::pairDomain(top, top->partner);
        }
    }
}

void Molecule::pairUncomplete(unsigned index) {

    auto strand = this->strands->at(index);
    bool binded = false;

    for (int i = strand->length() - 1; i >= 0; i--) {
        auto top    = strand->getAtom(i);
        
        if (top->partner != NULL && top->partner->partner == top) {
            binded = true;
        }

        if (top->partner != NULL && top->partner->partner == NULL && binded) {
            Strand::pairDomain(top, top->partner);
        }
    }

    binded = false;

    for (int i = 0; i < strand->length(); i++) {
        auto top    = strand->getAtom(i);
        
        if (top->partner != NULL && top->partner->partner == top) {
            binded = true;
        }

        if (top->partner != NULL && top->partner->partner == NULL && binded) {
            Strand::pairDomain(top, top->partner);
        }
    }
}

unsigned Molecule::size() {
    return this->strands->size();
}