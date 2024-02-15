#include <cstddef>
#include "molecule.h"

Molecule::Molecule(Strand *strand) {
    this->strands = new std::vector<Strand*>;
    this->strands->push_back(strand);
}

Molecule::~Molecule() {
    delete this->strands;
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
    Strand* tmp = this->strands->at(index);

    //get lastest index
    this->strands[index] = this->strands[this->strands->size() - 1];

    //and remove last
    this->strands->pop_back();
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

    for (unsigned i = 0; i < strand->length(); i++) {
        auto top    = strand->getAtom(i);
            
        if (top->partner != NULL && top->partner->partner == NULL) {
            Strand::pairDomain(top, top->partner);
        }
    }
}

void Molecule::unPair(unsigned index) {
    auto strand = this->strands->at(index);

    for (unsigned i = 0; i < strand->length(); i++) {
        auto top    = strand->getAtom(i);
            
        Strand::halfUnpairDomain(top, top->partner);
    }
}

unsigned Molecule::size() {
    return this->strands->size();
}