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

#include <openacc.h>
#include "openacc_curand.h"

#define TO_GPU_DOMAIN(X)  (DOMAIN_DT)(IS_COMPLEMENTARY(X) ? LAST_STRAND_INDEX - NORMALIZE_DOMAIN(X) : X)
#define ON_GPU_IS_COMP(X) ((LAST_STRAND_INDEX >> 1) < X)
#define ON_GPU_COMP(X)    (LAST_STRAND_INDEX - X)

#define randomGangSize   4
#define randomVectorSize 128
#define loopRandomIter   2

void Molecule::initGpu(std::map<DOMAIN_DT, Nucleotides*> &nucleotides) {
    this->Gmem  = new int[GPU_MEM_SIZE];
    this->onGPU = true;

    //domain repres
    //div 2 we dont need cmplementary part
    for (DOMAIN_DT i = 0; i < nucleotides.size() >> 1; i++) {
        auto nuc = nucleotides[i];

        Gmem[GET_GPU_COUNT_COL(i,                       GPU_NUCLEOTIDES)] = nuc->length();
        Gmem[GET_GPU_COUNT_COL(LAST_STRAND_INDEX - i,   GPU_NUCLEOTIDES)] = nuc->length(); //complementary

        for (unsigned j = 0; j < nuc->length(); j++) {
            Gmem[GPU_ELEMENT(j, i, GPU_NUCLEOTIDES)]                     =  nuc->get(j);  //Add standart
            Gmem[GPU_ELEMENT(j, LAST_STRAND_INDEX - i, GPU_NUCLEOTIDES)] = ~nuc->get(j); //Add complementary
        }
    }


    #pragma acc enter data copyin(this)
    #pragma acc enter data create(Gmem[0:GPU_MEM_SIZE])

    #pragma acc parallel loop present(Gmem)
    for (int i = 0; i < MAX_GPU_STRAND_LEN; i++) {

        #pragma acc loop
        for (int j = 0; j < MAX_GPU_STRANDS; j++) {
            Gmem[GPU_ELEMENT(j, i, GPU_NOTE)] = 0;
        }

    }
}

void Molecule::freeGpu() {
    #pragma acc exit data delete(Gmem)
    #pragma acc exit data delete(this)

    this->onGPU = false;

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
            Gmem[GPU_ELEMENT(atom_i, main_i, GPU_DOMAIN        )] = TO_GPU_DOMAIN(atom->domain.get());
        }
    }

    #pragma acc update device(Gmem[0:GPU_FRONT_MEM_SIZE])
}

inline bool GPUisPossibleToBind(int iAtom, int iStrand, int jAtom, int jStrand, int* Gmem) {
    bool status = Gmem[GPU_ELEMENT(iAtom, iStrand, GPU_PARTNER_STRAND)] == -1 && 
                  Gmem[GPU_ELEMENT(jAtom, jStrand, GPU_PARTNER_STRAND)] == -1 &&
                  jStrand != iStrand;

    return status;
}


inline void GPUdeterministicBind(int* Gmem) {    
    #pragma acc parallel loop present(Gmem)
    for (int main_i = 0; main_i < Gmem[GET_GPU_COUNT_ROW()]; main_i++) {

        #pragma acc loop
        for (int strand_i = 0; strand_i < Gmem[GET_GPU_COUNT_ROW()]; strand_i++) {
            auto main_len = Gmem[GET_GPU_COUNT_COL(main_i, GPU_PARTNER_STRAND)];
            auto strand_len = Gmem[GET_GPU_COUNT_COL(strand_i, GPU_PARTNER_STRAND)];

            int  strandBS = -1; 
            int  mainBS   = -1;

            //find bind spot
            #pragma acc loop seq
            for (int j = 0; j < strand_len; j++) {
                if (
                    Gmem[GPU_ELEMENT(j, strand_i, GPU_PARTNER_STRAND)] == main_i
                ) {
                    strandBS = j;
                    mainBS   = Gmem[GPU_ELEMENT(j, strand_i, GPU_PARTNER_ATOM)];
                }
            }

            if (strandBS == -1) continue;

            #pragma acc loop seq
            for (int j = strandBS; j < strand_len; j++) {
                auto k = mainBS + (j - strandBS);

                if (k < main_len && Gmem[GPU_ELEMENT(j, strand_i, GPU_PARTNER_STRAND)] == -1) {
                    if (Gmem[GPU_ELEMENT(k, main_i, GPU_DOMAIN)] == ON_GPU_COMP(Gmem[GPU_ELEMENT(j, strand_i, GPU_DOMAIN)])) {
                        if (Gmem[GPU_ELEMENT(k, main_i, GPU_PARTNER_STRAND)] == -1) { // do bind
                            
                            Gmem[GPU_ELEMENT(k, main_i, GPU_PARTNER_STRAND)]   = strand_i;
                            Gmem[GPU_ELEMENT(k, main_i, GPU_PARTNER_ATOM)]     = j;

                            Gmem[GPU_ELEMENT(j, strand_i, GPU_PARTNER_STRAND)] = main_i;
                            Gmem[GPU_ELEMENT(j, strand_i, GPU_PARTNER_ATOM)]   = k;
                        
                        } else if (GPU_ELEMENT(k, main_i, GPU_PARTNER_STRAND) != strand_i) { // do ubnind of existing but not bind new and remeber it for second loop
                            auto partnerStrand = Gmem[GPU_ELEMENT(k, main_i, GPU_PARTNER_STRAND)];
                            auto partnerAtom   = Gmem[GPU_ELEMENT(k, main_i, GPU_PARTNER_ATOM)];

                            Gmem[GPU_ELEMENT(k, main_i, GPU_NOTE)]                  = 100;
                            Gmem[GPU_ELEMENT(partnerAtom, partnerStrand, GPU_NOTE)] = 100;
                        }
                    }
                }
            }

            #pragma acc loop seq
            for (int j = strandBS; j >= 0; j--) {
                auto k = mainBS - (strandBS - j);

                if (k >= 0 && Gmem[GPU_ELEMENT(j, strand_i, GPU_PARTNER_STRAND)] == -1) {
                    if (Gmem[GPU_ELEMENT(k, main_i, GPU_DOMAIN)] == ON_GPU_COMP(Gmem[GPU_ELEMENT(j, strand_i, GPU_DOMAIN)])) {
                        if (Gmem[GPU_ELEMENT(k, main_i, GPU_PARTNER_STRAND)] == -1) { // do bind
                            
                            Gmem[GPU_ELEMENT(k, main_i, GPU_PARTNER_STRAND)]   = strand_i;
                            Gmem[GPU_ELEMENT(k, main_i, GPU_PARTNER_ATOM)]     = j;

                            Gmem[GPU_ELEMENT(j, strand_i, GPU_PARTNER_STRAND)] = main_i;
                            Gmem[GPU_ELEMENT(j, strand_i, GPU_PARTNER_ATOM)]   = k;

                        } else if (GPU_ELEMENT(k, main_i, GPU_PARTNER_STRAND) != strand_i) { // do ubnind of existing but not bind new
                            auto partnerStrand = Gmem[GPU_ELEMENT(k, main_i, GPU_PARTNER_STRAND)];
                            auto partnerAtom   = Gmem[GPU_ELEMENT(k, main_i, GPU_PARTNER_ATOM)];

                            Gmem[GPU_ELEMENT(k, main_i, GPU_NOTE)]                  = 100;
                            Gmem[GPU_ELEMENT(partnerAtom, partnerStrand, GPU_NOTE)] = 100;
                        }
                    }
                }
            }
        }
    }

    #pragma acc parallel loop present(Gmem)
    for (int main_i = 1; main_i < Gmem[GET_GPU_COUNT_ROW()]; main_i++) {
        auto main_len = Gmem[GET_GPU_COUNT_COL(main_i, GPU_PARTNER_STRAND)];

        int bind_size = 0;
        bool unbindforce = false;
        
        #pragma acc loop seq
        for (int atom_i = 0; atom_i < main_len; atom_i++) {
            bind_size += Gmem[GPU_ELEMENT(atom_i, main_i, GPU_NOTE)] != 100 && Gmem[GPU_ELEMENT(atom_i, main_i, GPU_PARTNER_STRAND)] ==  0;
            unbindforce = unbindforce || Gmem[GPU_ELEMENT(atom_i, main_i, GPU_NOTE)] == 100 && Gmem[GPU_ELEMENT(atom_i, main_i, GPU_PARTNER_STRAND)] != -1;

            Gmem[GPU_ELEMENT(atom_i, main_i, GPU_NOTE)] = 0;
        }

        #pragma acc loop seq
        for (int atom_i = 0; atom_i < main_len; atom_i++) {
            // i partner
            int iPStrand = Gmem[GPU_ELEMENT(atom_i, main_i, GPU_PARTNER_STRAND)];
            int iPAtomId = Gmem[GPU_ELEMENT(atom_i, main_i, GPU_PARTNER_ATOM)];

            if ((bind_size < 1 || (bind_size <= 1 && unbindforce)) && iPStrand != -1) {
                Gmem[GPU_ELEMENT(atom_i, main_i, GPU_PARTNER_STRAND)]       = -1;
                Gmem[GPU_ELEMENT(atom_i, main_i, GPU_PARTNER_ATOM)]         = -1;

                Gmem[GPU_ELEMENT(iPAtomId, iPStrand,   GPU_PARTNER_STRAND)] = -1;
                Gmem[GPU_ELEMENT(iPAtomId, iPStrand,   GPU_PARTNER_ATOM)]   = -1;
            }
        }
    }
}

#pragma acc routine seq nohost
inline int randomRange(curandState_t* state, int rangeI) {
    return ((int)(curand_uniform(state) * (float)rangeI));
}

void Molecule::simulateGpu(std::map<DOMAIN_DT, Nucleotides*> &nucleotides, float temp, std::vector<Molecule*> &instruction, unsigned strands_count, unsigned sim_time) {
    // insert instructions strands
    for (unsigned i = 0; i < strands_count; i++) {
        for (const auto &s : instruction) {
            this->strands->push_back(s->getStrand(0)->copy());
        }
    }   

    //update GPU DATA
    this->updateGpu();

    curandState_t states[randomGangSize*randomVectorSize];
    #pragma acc enter data create(states[0:randomGangSize*randomVectorSize])

    //init states
    uint64_t seed = ((uint64_t)time(NULL)) * 3817290472947109740298748;

    #pragma acc parallel loop present(states) copyin(seed) num_gangs(randomGangSize) vector_length(randomVectorSize)
    for (unsigned i = 0; i < randomGangSize; i++) {

        #pragma acc loop vector
        for (unsigned j = 0; j < randomVectorSize; j++) {
            curand_init(seed * i, j , 0, &states[__pgi_gangidx() * randomVectorSize + __pgi_vectoridx()]);
        }
    }

    for (int iter = 0; iter < sim_time >> 11; iter++) {
        #pragma acc parallel loop present(Gmem, states) num_gangs(randomGangSize) vector_length(randomVectorSize)
        for (unsigned i = 0; i < randomGangSize; i++) {

            #pragma acc loop vector independent
            for (unsigned threat = 0; threat < randomVectorSize; threat++) {

                curandState_t state = states[__pgi_gangidx() * randomVectorSize + __pgi_vectoridx()];

                #pragma acc loop seq
                for (unsigned loop = 0; loop < loopRandomIter; loop++) {

                    // select i and j strands
                    int len     = Gmem[GET_GPU_COUNT_ROW()];

                    int iStrand = randomRange(&state, len);
                    int jStrand = randomRange(&state, len);

                    // select i and j atom
                    int iAtomId = randomRange(&state, Gmem[GET_GPU_COUNT_COL(iStrand, GPU_PARTNER_STRAND)]);
                    int jAtomId = randomRange(&state, Gmem[GET_GPU_COUNT_COL(jStrand, GPU_PARTNER_STRAND)]);

                    // domain
                    int iDomain        = Gmem[GPU_ELEMENT(iAtomId, iStrand, GPU_DOMAIN)];
                    int jDomain        = Gmem[GPU_ELEMENT(jAtomId, jStrand, GPU_DOMAIN)];

                    // i partner
                    int iPStrand       = Gmem[GPU_ELEMENT(iAtomId, iStrand, GPU_PARTNER_STRAND)];
                    int iPAtomId       = Gmem[GPU_ELEMENT(iAtomId, iStrand, GPU_PARTNER_ATOM)];

                    if (curand_uniform(&state) < 0.5) { // renaturation
                        if (
                            Nucleotides::GpuReNaturationP(iDomain, jDomain, Gmem, temp) > curand_uniform(&state) &&
                            GPUisPossibleToBind(iAtomId, iStrand, jAtomId, jStrand, Gmem)
                        ) {
                            Gmem[GPU_ELEMENT(iAtomId, iStrand, GPU_PARTNER_STRAND)] = jStrand;
                            Gmem[GPU_ELEMENT(jAtomId, jStrand, GPU_PARTNER_STRAND)] = iStrand;

                            Gmem[GPU_ELEMENT(iAtomId, iStrand, GPU_PARTNER_ATOM  )] = jAtomId;
                            Gmem[GPU_ELEMENT(jAtomId, jStrand, GPU_PARTNER_ATOM  )] = iAtomId;
                        }
                    } else if (
                        iPStrand != -1 &&
                        Nucleotides::GpuDeNaturationP(iDomain, Gmem[GPU_ELEMENT(iPAtomId, iPStrand, GPU_DOMAIN)], Gmem, temp) > curand_uniform(&state)
                    ) { // deneturartion
                        Gmem[GPU_ELEMENT(iAtomId, iStrand, GPU_PARTNER_STRAND)]     = -1;
                        Gmem[GPU_ELEMENT(iAtomId, iStrand, GPU_PARTNER_ATOM)]       = -1;

                        Gmem[GPU_ELEMENT(iPAtomId, iPStrand,   GPU_PARTNER_STRAND)] = -1;
                        Gmem[GPU_ELEMENT(iPAtomId, iPStrand,   GPU_PARTNER_ATOM)]   = -1;
                    }
                }

                //#pragma acc atomic write
                states[__pgi_gangidx() * randomVectorSize + __pgi_vectoridx()] = state;
            }
        }

        GPUdeterministicBind(Gmem);
    }

    #pragma acc exit data delete(states[0:randomGangSize*randomVectorSize])

    //UPDATE CPU data
    this->updateSelf();
}