#include "strand.h"
#include "error.h"
#include <iostream>
#include <stdexcept>

/* Definice zarážky */
Atom multiAtomStack = {0, NULL, NULL, 0};
Atom *multiAtom     = &multiAtomStack;

Strand::Strand(std::map<DOMAIN_DT, Nucleotides*> *nucleotides) {
    this->atoms       = new std::vector<Atom>;
    this->deleted     = false;
    this->readOnly    = false;
    this->nucleotides = nucleotides;
}

Strand::~Strand() {
    delete this->atoms;
}

unsigned Strand::length() {
    return this->atoms->size();
}

Atom* Strand::getAtom(unsigned index) {
    this->readOnly = true;

    if (index >= this->atoms->size()) error::inAppError("Indexing out of range in strand"); 

    return &((*this->atoms)[index]);
}

Strand* Strand::copy() {
    auto str = new Strand(this->nucleotides);

    for (auto &atom : *this->atoms) {
        str->addDomain(atom.domain);   
    }

    return str;
}

void Strand::del() {
    this->deleted = true;
}

bool Strand::isDeleted() {
    return this->deleted;
}

unsigned Strand::addDomain(Domain d) {
    if (this->readOnly) {
        throw std::runtime_error("Try to add domain to read only strand");
    }

    auto atom = this->createAtom(d);

    if (this->atoms->size() > 0 && this->nucleotides != NULL) {
        auto last = (*this->atoms)[this->atoms->size() - 1];

        if (this->nucleotides->count(NORMALIZE_DOMAIN(last.domain.get()))) {
            atom.offsetN = last.offsetN + ((*(this->nucleotides))[NORMALIZE_DOMAIN(last.domain.get())])->length();
        }
    }

    this->atoms->push_back(atom);

    return this->atoms->size() - 1;
}

void Strand::pairDomain(Atom* atom1, Atom* atom2) {
    //do link a1 <-> a2
    atom1->partner = atom2;
    atom2->partner = atom1;
}

void Strand::halfPairDomain(Atom* atom1, Atom* atom2) {
    //do link a1 -> a2
    atom1->partner        = atom2;
    
    atom1->partnersCount += 1;
    atom2->partnersCount += 1;
}

void Strand::unpairDomain(Atom* atom) {
    if (atom->partner == NULL)           return;
    if (atom->partner->partner == atom)  atom->partner->partner = NULL;

    atom->partner->partnersCount -= 1;
    atom->partnersCount          -= 1;

    atom->partner = NULL;
}

void Strand::complementaryLast() {
    this->atoms->at(this->length() - 1).domain = ~(this->atoms->at(this->length() - 1).domain);
}

Atom Strand::createAtom(Domain d) {
    Atom a = {d, NULL, this, 0, 0};

    return a;
}