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
        for (auto &inst : instruction) {
            this->removeUnbinded(inst);
        }

        auto preAddRegSize = this->reg->size();

        for (auto &inst : instruction) {
            this->doAllBinding(inst);
        }

        this->removeReplaced(preAddRegSize);

        this->removeUnstable(preAddRegSize - 1);

        for (unsigned i = 0; i < this->reg->size(); i++) {
            this->reg->pairUncomplete(i);
        }

        this->removeUnstable(this->reg->size() - 1);

        for (unsigned i = 0; i < this->reg->size(); i++) {
            this->reg->pairUncomplete(i);
        }
        
        this->reg->finishDelete();
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
    //bind all incomming
    for (unsigned i = oldSize; i < this->reg->size(); i++) {
        auto workStrand = this->reg->getStrand(i);

        // find binding spot for strand
        for (unsigned bsi = 0; bsi < workStrand->length(); bsi++) {

            if (workStrand->getAtom(bsi)->partner != NULL &&
                workStrand->getAtom(bsi)->partner->partner == NULL) { //apply unbind force

                // bind on binding spot
                workStrand->getAtom(bsi)->partner->partner = workStrand->getAtom(bsi);
                
                break;

            }

        }
    }

    //apply unbind force
    for (unsigned i = oldSize; i < this->reg->size(); i++) {
        auto workStrand = this->reg->getStrand(i);

        // find binding spot for strand
        for (unsigned bsi = 0; bsi < workStrand->length(); bsi++) {

            if (workStrand->getAtom(bsi)->partner != NULL &&
                workStrand->getAtom(bsi)->partner->partner == workStrand->getAtom(bsi)) {

                    for (unsigned j = bsi + 1; j < workStrand->length() - bsi; j++) {
                        if (workStrand->getAtom(j)->partner != NULL) workStrand->getAtom(j)->partner->partner = NULL; // do unbind
                        else break; // end of binding
                    }

                    for (int j = bsi - 1; j >= 0; j--) {
                        if (workStrand->getAtom(j)->partner != NULL) workStrand->getAtom(j)->partner->partner = NULL; // do unbind
                        else break; // end of binding
                    }
                }
        }
    }
}

void Register::removeUnstable(int to) {
    for (int i = to; i >= 1; i--) {
        auto mainStrand     = this->reg->getStrand(i);
        unsigned bindScore  = 0;

        // calculate bind score
        for (int i = 0; i < mainStrand->length(); i++) {
            if (mainStrand->getAtom(i)->partner != NULL &&
                mainStrand->getAtom(i)->partner->partner == mainStrand->getAtom(i)) {
                bindScore++;
            } else {
                bindScore = 0;
            }

            if (bindScore >= 2) break;
        }

        //this do unbind?
        if (bindScore < 2) {
            this->reg->remStrand(i);
        }
    }
}

void Register::removeUnbinded(Molecule* mol) {
    if (mol->size() == 0) return; // nothing to bind

    auto strand = mol->getStrand(0);

    for (int j = this->reg->size() - 1; j >= 1; j--) {
        auto mainStrand     = this->reg->getStrand(j);
        unsigned bindScore  = 0;

        // calculate bind score
        for (int i = 0; i < mainStrand->length(); i++) {
            if (mainStrand->getAtom(i)->partner != NULL &&
                mainStrand->getAtom(i)->partner->partner == mainStrand->getAtom(i)) {
                bindScore++;
            }
        }


        // for all possible bindigs
        for (int i = 1 - (int)strand->length(); i < (int)mainStrand->length(); i++) {

            unsigned alignScore   = 0;
            unsigned pullingScore = 0;
            unsigned compareEnd   = std::min(strand->length(), mainStrand->length() - i);
            bool     canPull      = false;

            //try aligment it and check score
            for (int baseID = i < 0 ? std::abs(i) : 0; baseID < compareEnd; baseID++) {
                
                if (~(strand->getAtom(baseID)->domain) == mainStrand->getAtom(i + baseID)->domain) {

                    if (mainStrand->getAtom(i + baseID)->partner == NULL || mainStrand->getAtom(i + baseID)->partner->partner != mainStrand->getAtom(i + baseID)) {
                        canPull = true;
                    } else {
                        pullingScore += 1;
                    }

                    alignScore += 1;

                }
            }

            //this do unbind?
            if (canPull && alignScore >= 2 && (bindScore - pullingScore) <= 1) {
                this->reg->remStrand(j);
                break;
            }
        }

    }



}