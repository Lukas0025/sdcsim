#include <cstddef>
#include "molecule.h"
#include <iostream>

#define NO_PAIR -1

Molecule::Molecule(Strand *strand) {
    this->strands = new std::vector<Strand*>;
    
    this->nucleotides         = new std::vector<std::vector<NUCLEOTIDE_DT>>;
    this->nucleotidesBindings = new std::vector<std::vector<int>>;
    this->nucleotides2domain  = new std::vector<Strand*>;
    
    this->strands->push_back(strand);
}

Molecule::~Molecule() {
    delete this->strands;
}

void Molecule::createNucleotidesLevel(std::map<DOMAIN_DT, Nucleotides*> &nucleotides) { // convert current reprezetation to nucleodies level
    // Create map 
    // for performance
    std::map<Strand*, int> strandsMap;
    for (unsigned i = 0; i < this->size(); i++) {
        strandsMap[this->getStrand(i)] = i;
    }

    for (unsigned i = 0; i < this->size(); i++) {
        auto strand   = this->getStrand(i);
        
        std::vector<NUCLEOTIDE_DT> strandN;
        std::vector<int>           bindings;

        for (unsigned j = 0; j < strand->length(); j++) {
            auto atom = strand->getAtom(j);

            if (!nucleotides.count(NORMALIZE_DOMAIN(atom->domain.get()))) std::cout << "error no domain for " << NORMALIZE_DOMAIN(atom->domain.get()) << "\n";

            auto nuc = nucleotides[NORMALIZE_DOMAIN(atom->domain.get())];
            
            bool complementary = IS_COMPLEMENTARY(atom->domain.get());
            unsigned partner   = HAVE_PARTNER(atom) ? strandsMap[atom->partner->strand]  : NO_PAIR;
            unsigned pos       = HAVE_PARTNER(atom) ? atom->partner->offsetN    : j;

            for (unsigned k = 0; k < nuc->length(); k++) {
                strandN.push_back(complementary ? ~nuc->get(k) : nuc->get(k));
                bindings.push_back(partner);
                bindings.push_back(pos + k);
            }
        }

        this->nucleotides->push_back(strandN);
        this->nucleotidesBindings->push_back(bindings);
        this->nucleotides2domain->push_back(strand);
    }
}

void Molecule::addNucleotidesFreeStrand(Strand* strand, unsigned density, std::map<DOMAIN_DT, Nucleotides*> &nucleotides) {
    for (unsigned i = 0; i < density; i++) {
        
        std::vector<NUCLEOTIDE_DT> strandN;
        std::vector<int>           bindings;

        for (unsigned j = 0; j < strand->length(); j++) {
            auto atom = strand->getAtom(j);

            if (!nucleotides.count(NORMALIZE_DOMAIN(atom->domain.get()))) std::cout << "error no domain for " << NORMALIZE_DOMAIN(atom->domain.get()) << "\n";

            auto nuc           = nucleotides[NORMALIZE_DOMAIN(atom->domain.get())];
            bool complementary = IS_COMPLEMENTARY(atom->domain.get());

            for (unsigned k = 0; k < nuc->length(); k++) {
                strandN.push_back(complementary ? ~nuc->get(k) : nuc->get(k));
                bindings.push_back(NO_PAIR);
                bindings.push_back(j + k);
            }
        }

        this->nucleotides->push_back(strandN);
        this->nucleotidesBindings->push_back(bindings);
        this->nucleotides2domain->push_back(strand);
    }
}

void Molecule::updateDomainLevel() {      // convert current reprezentation to domain level

}

/*void Molecule::simulate() {

}*/

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