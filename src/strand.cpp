#include "strand.h"
#include <iostream>
#include <stdexcept>

Strand::Strand() {
    this->atoms = new std::vector<Atom>;
    this->readOnly = false;
}

Strand::~Strand() {
    delete this->atoms;
}

unsigned Strand::length() {
    return this->atoms->size();
}

Atom* Strand::getAtom(unsigned index) {
    this->readOnly = true;

    return &((*this->atoms)[index]);
}

unsigned Strand::addDomain(Domain d) {
    if (this->readOnly) {
        throw std::runtime_error("Try to add domain to read only strand");
    }

    auto atom = this->createAtom(d);

    this->atoms->push_back(atom);

    return this->atoms->size() - 1;
}

void Strand::pairDomain(Atom* atom1, Atom* atom2) {
    //do cross linking
    atom1->partner = atom2;
    atom2->partner = atom1;
}

void Strand::alignmentPrint() {
    for (auto atom = this->atoms->begin(); atom < this->atoms->end(); atom++) {
        std::cout << atom->domain.get();
        std::cout << " ";
    }
}

Atom Strand::createAtom(Domain d) {
    Atom a = {d, NULL, this};

    return a;
}