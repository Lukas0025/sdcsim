/**
 * SDCSIM - strand displacement simulator on SIMD||DNA architecture
 * @file molecule.cpp
 * @brief Contain implementation of molecule class
 * @author Lukáš Plevač <xpleva07@vut.cz>
 */

#include <cstddef>
#include "molecule.h"
#include "error.h"
#include <cstdlib>
#include <iostream>
#include <random>

#define RANDOM_RANGE(TO)       (mt() % (TO))
#define RANDOM_FLOAT()         ((float)mt()/(float)mt.max())

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

void Molecule::simulate(std::map<DOMAIN_DT, Nucleotides*> &nucleotides, float temp, std::vector<Molecule*> &instruction, unsigned strands_count, unsigned sim_time) {
    // set random seed
    std::mt19937 mt(time(nullptr)); 

    // insert instructions strands
    for (unsigned i = 0; i < strands_count; i++) {
        for (const auto &s : instruction) {
            this->strands->push_back(s->getStrand(0)->copy());
        }
    }

    for (int i = 0; i < sim_time; i++) {
        //select i and j strands
        auto iStrand = this->getStrand(RANDOM_RANGE(this->size()));
        auto jStrand = this->getStrand(RANDOM_RANGE(this->size()));

        //select i and j atom
        auto iAtomId = RANDOM_RANGE(iStrand->length());
        auto iAtom   = iStrand->getAtom(iAtomId);

        auto jAtomId = RANDOM_RANGE(jStrand->length());
        auto jAtom   = jStrand->getAtom(jAtomId);
        
        if (RANDOM_FLOAT() < 0.5) { // renaturation
            if (nucleotides[iAtom->domain.get()]->ReNaturationP(nucleotides[jAtom->domain.get()], temp) > RANDOM_FLOAT()) {
                if (isPossibleToBind(iAtom, jAtom)) {
                    Strand::pairDomain(iAtom, jAtom);
                }
            }
        } else if (iAtom->partner != NULL && nucleotides[iAtom->domain.get()]->DeNaturationP(nucleotides[iAtom->partner->domain.get()], temp) > RANDOM_FLOAT()) {
            iAtom->partner->partner = NULL;
            iAtom->partner          = NULL;
        }

        if (i % 100 == 99) {
            this->deterministicBind();
        }
    }

    this->deterministicBind(); //do correct unbind
}

void Molecule::deterministicBind() {

    //#pragma omp parallel for
    for (int main_i = 0; main_i < this->size(); main_i++) {
        auto mainStrand = this->getStrand(main_i);

        for (int i = 0; i < this->size(); i++) {
            auto strand = this->getStrand(i);
            int  strandBS = -1; 
            int  mainBS   = -1;

            //find bind spot
            for (int j = 0; j < strand->length(); j++) {
                auto atom = strand->getAtom(j);
                if (atom->partner != NULL && atom->partner->strand == mainStrand) {
                    strandBS = j;
                    mainBS   = atom->partner - atom->partner->strand->getAtom(0);
                }
            }

            if (strandBS == -1) continue;

            for (int j = strandBS; j < strand->length(); j++) {
                auto k = mainBS + (j - strandBS);

                if (k >= mainStrand->length()) break;
                if (strand->getAtom(j)->partner != NULL) continue;

                if (mainStrand->getAtom(k)->domain == ~strand->getAtom(j)->domain) {
                    if (mainStrand->getAtom(k)->partner == NULL) { // do bind
                        Strand::pairDomain(mainStrand->getAtom(k), strand->getAtom(j));
                    } else if (mainStrand->getAtom(k)->partner->strand != strand) { // do ubnind of existing but not bind new and remeber it for second loop
                        mainStrand->getAtom(k)->partnersCount = 100;
                        mainStrand->getAtom(k)->partner->partnersCount = 100;
                    }
                }
            }


            for (int j = strandBS; j >= 0; j--) {
                auto k = mainBS - (strandBS - j);

                if (k < 0) break;
                if (strand->getAtom(j)->partner != NULL) continue;

                if (mainStrand->getAtom(k)->domain == ~strand->getAtom(j)->domain) {
                    if (mainStrand->getAtom(k)->partner == NULL) { // do bind
                        Strand::pairDomain(mainStrand->getAtom(k), strand->getAtom(j));
                    } else if (mainStrand->getAtom(k)->partner->strand != strand) { // do ubnind of existing but not bind new
                        mainStrand->getAtom(k)->partnersCount = 100;
                        mainStrand->getAtom(k)->partner->partnersCount = 100;
                    }
                }
            }
        }
    }

    for (int main_i = 1; main_i < this->size(); main_i++) {
        auto mainStrand = this->getStrand(main_i);

        int bind_size = 0;
        bool unbindforce = false;
            
        for (int i = 0; i < mainStrand->length(); i++) {
            if (mainStrand->getAtom(i)->partnersCount != 100 &&
                mainStrand->getAtom(i)->partner != NULL &&
                mainStrand->getAtom(i)->partner->partner == mainStrand->getAtom(i) &&
                mainStrand->getAtom(i)->partner->strand == this->getStrand(0)) bind_size += 1;

            if (mainStrand->getAtom(i)->partnersCount == 100 &&
                mainStrand->getAtom(i)->partner != NULL) unbindforce = true;

            mainStrand->getAtom(i)->partnersCount = 0;
        }

        if (bind_size < 1 || (bind_size <= 1 && unbindforce)) {
            for (int i = 0; i < mainStrand->length(); i++) {
                Strand::unpairDomain(mainStrand->getAtom(i));
            }
        }
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
    for (unsigned i = mainBindFrom; i < std::min(strand->length() + mainBindFrom - incomeBindFrom, mainStrand->length()); i++) {
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
    if (index >= this->strands->size()) error::inAppError("Indexing out of range in molecule"); 

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