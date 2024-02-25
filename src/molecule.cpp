#include <cstddef>
#include "molecule.h"
#include <iostream>

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

    if (index == 0) return; // main strand is nundeleteble

    //do unlink
    for (unsigned i = 0; i < this->strands->at(index)->length(); i++) {
        Strand::unpairDomain(this->strands->at(index)->getAtom(i));
    }

    // set deleted flag
    this->strands->at(index)->del();
}

void Molecule::finishDelete() {
    for (int i = this->strands->size() - 1; i >= 0; i--) {
        if (this->strands->at(i)->isDeleted()) {
            //get lastest index
            this->strands->at(i) = this->strands->at(this->strands->size() - 1);

            //and remove last
            this->strands->pop_back();
        }
    }
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