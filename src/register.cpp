#include "register.h"
#include <cmath>
#include <iostream>

Register::Register(Molecule *init) {
    this->reg = init;
}

Molecule* Register::get() {
    return this->reg;
}

void Register::applyInstruction(std::vector<Molecule*> instruction) {
    
    for (unsigned i = 0; i < 1; i++) {
        auto preAddRegSize = this->reg->size();

        for (auto &inst : instruction) {
            this->doAllBinding(inst);
        }

        this->removeReplaced(preAddRegSize);

        this->removeUnstable();

        for (unsigned i = 0; i < this->reg->size(); i++) {
            this->reg->pairUncomplete(i);
        }
    }
    
}

void Register::doAllBinding(Molecule* mol) {

    if (mol->size() == 0) return; // nothing to bind

    auto strand     = mol->getStrand(0);
    auto mainStrand = this->reg->getStrand(0);

    // for all possible bindigs
    for (int i = 1 - (int)strand->length(); i < (int)mainStrand->length(); i++) {

        unsigned alignScore = 0;
        unsigned compareEnd = std::min(strand->length(), mainStrand->length() - i);

        //try aligment it and check score
        for (int baseID = i < 0 ? std::abs(i) : 0; baseID < compareEnd; baseID++) {
            
            if (~(strand->getAtom(baseID)->domain) == mainStrand->getAtom(i + baseID)->domain) {
                alignScore += 1;
            } else {
                alignScore  = 0;
            }

            // this definitly bind
            if (alignScore >= 2) {
                auto newStrand = strand->copy();
                this->reg->addStrand(newStrand, i);
                break;
            }

        }
    }
}

void Register::removeReplaced(unsigned oldSize) {

}

void Register::removeUnstable() {

}