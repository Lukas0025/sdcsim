/**
 * SDCSIM - strand displacement simulator on SIMD||DNA architecture
 * @file molecule_gpu.cpp
 * @brief Contain gpu implementation of molecule class
 * @author Lukáš Plevač <xpleva07@vut.cz>
 */

#include <cstddef>
#include "molecule.h"
#include "error.h"
#include <cstdlib>
#include <iostream>
#include <random>

void Molecule::initGpu() {
    this->Gmem = new int[GPU_MEM_SIZE];

    #pragma acc enter data copyin(this)
    #pragma acc enter data create(Gmem[0:GPU_MEM_SIZE])
}

void Molecule::freeGpu() {
    #pragma acc exit data delete(Gmem)
    #pragma acc exit data delete(this)

    delete Gmem;
}

void Molecule::updateSelf() {
    #pragma acc update self(Gmem[0:GPU_BACK_MEM_SIZE])

    std::map<int, Strand*> strandPointers;

    for (int main_i = 0; main_i < this->size(); main_i++) {
        strandPointers[main_i] = this->getStrand(main_i);
    }

    for (int main_i = 0; main_i < this->size(); main_i++) {
        auto mainStrand = this->getStrand(main_i);

        for (int atom_i = 0; atom_i < mainStrand->length(); atom_i++) {
            auto atom = mainStrand->getAtom(atom_i);
            auto partnerStrand = Gmem[GPU_ELEMENT(atom_i, main_i, GPU_PARTNER_STRAND)] == -1 ? 
                NULL : strandPointers[Gmem[GPU_ELEMENT(atom_i, main_i, GPU_PARTNER_STRAND)]];

            atom->partner = NULL;

            if (partnerStrand != NULL)
                atom->partner = partnerStrand->getAtom(Gmem[GPU_ELEMENT(atom_i, main_i, GPU_PARTNER_ATOM)]);
        }
    }
}

void Molecule::updateGpu() {
    Gmem[GET_GPU_COUNT_ROW()] = this->size();

    std::map<Strand*, int> strandIndexes;

    //cretae simple map
    for (int main_i = 0; main_i < this->size(); main_i++) {
        strandIndexes[this->getStrand(main_i)] = main_i;
    }

    for (int main_i = 0; main_i < this->size(); main_i++) {
        auto mainStrand = this->getStrand(main_i);    

        Gmem[GET_GPU_COUNT_COL(main_i, GPU_PARTNER_STRAND)] = mainStrand->length();
        Gmem[GET_GPU_COUNT_COL(main_i, GPU_PARTNER_ATOM  )] = mainStrand->length();
        Gmem[GET_GPU_COUNT_COL(main_i, GPU_DOMAIN        )] = mainStrand->length();

        for (int atom_i = 0; atom_i < mainStrand->length(); atom_i++) {
            auto atom = mainStrand->getAtom(atom_i);
            Gmem[GPU_ELEMENT(atom_i, main_i, GPU_PARTNER_STRAND)] = atom->partner == NULL ? -1 : strandIndexes[atom->partner->strand];
            Gmem[GPU_ELEMENT(atom_i, main_i, GPU_PARTNER_ATOM  )] = atom->partner == NULL ? -1 : atom->partner - atom->partner->strand->getAtom(0);
            Gmem[GPU_ELEMENT(atom_i, main_i, GPU_DOMAIN        )] = atom->domain.get();
        }
    }

    #pragma acc update device(Gmem[0:GPU_MEM_SIZE])
}