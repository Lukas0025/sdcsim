/**
 * SDCSIM - strand displacement simulator on SIMD||DNA architecture
 * @file nucleotides.h
 * @brief Contain headers for nucleotides class
 * @author Lukáš Plevač <xpleva07@vut.cz>
 */

#pragma once

#include <vector>
#include <string>
#include <math.h>
#include <cstdint>
#include "gpu.h"

#define NUCLEOTIDE_DT uint8_t

// purines
#define ADENINE   0b01
#define GUANINE   0b10

// pyrmides
#define CYTOSINE  0xFD
#define THYMINE   0xFE

#define boltzmann_constant 1.987204259e-3

// ~pyrmides == purines <=> pyrmides == ~purines

/**
 * class reprezenting domain nucleotides reprezentation
 */
class Nucleotides {
    public:
        Nucleotides();
        ~Nucleotides();

        /**
         * Add nucleotide by string reprezentation (append)
         * @param na char reprezenting nucleotide
         * @param complementry bool say if compute of complementary before add (if true and you add A you in real adding T)
         */
        void addFromStr(char na, bool complementary = false);

        /**
         * Add nucleotide to domain reprezenation (append)
         * @param na NUCLEOTIDE_DT to add
         */
        void add(NUCLEOTIDE_DT na);

        /**
         * length of nucleotides strand reprezenting domain
         */
        unsigned length();

        /**
         * Get nucleotide on position
         * @param pos position in domain
         */
        NUCLEOTIDE_DT get(unsigned pos);

        /**
         * Compute bind power between to nucleotides
         * @param n1 nucleotide 
         * @param n2 nucleotide
         * @return energy in kcal
         */
        static inline float bindPower(NUCLEOTIDE_DT n1, NUCLEOTIDE_DT n2) {
            //Klump and Ackermann 1971 data
            if ((n1 == GUANINE && n2 == CYTOSINE) || (n2 == GUANINE && n1 == CYTOSINE)) {
                return -9.0; // 9kcal
            } else if ((n1 == ADENINE && n2 == THYMINE) || (n2 == ADENINE && n1 == THYMINE)) {
                return -7.2; // 7.2kcal
            } else {
                return -5.4; // 5.4kcal
            }
        }

        /**
         * Get string of nucleotides strand
         * @return string reprezentation
         */
        std::string getStr();

        /**
         * Compute complement of whole strand
         */
        Nucleotides operator~() {
            auto nuc = Nucleotides();

            for (const auto& val : *(this->strand)) {
                nuc.add(~val);
            }

            return nuc;
        }

        /**
         * Compute of probability of renaturating between to domains in Nucleotides reprezentation
         * @param partner partner to bind with
         * @param temp    temperature of solution
         * @return float probability 0..1
         */
        float ReNaturationP(Nucleotides* partner, float temp);

        /**
         * Compute of probability of denaturatin between to domains in Nucleotides reprezentation
         * @param partner partner to denaturate with
         * @param temp    temperature of solution
         * @return float probability 0..1
         */
        float DeNaturationP(Nucleotides* partner, float temp);

        /**
         * DeltaH between to domains in Nucleotides reprezentation
         * @param partner partner to compute
         * @return float
         */
        float deltaH(Nucleotides* partner);
        
        /**
         * DeltaS between to domains in Nucleotides reprezentation
         * @param partner partner to compute
         * @return float
         */
        float deltaS(Nucleotides* partner);

        /**
        * GPU Implementations
        */

        static inline float GpuReNaturationP(int d1, int d2, int* gpuData, float temp) {
            temp += 273.15; // convert to kelvin
            
            return std::min(1., 
                exp(
                    -(
                        (Nucleotides::GPUDeltaH(d1, d2, gpuData) + temp * Nucleotides::GPUDeltaS(d1, d2, gpuData)) /
                        (boltzmann_constant * temp)
                    )
                )
            );
        }

        static inline float GpuDeNaturationP(int d1, int d2, int* gpuData, float temp) {
            temp += 273.15; // convert to kelvin
            
            return std::min(1., 
                exp(
                    -(
                        ((-Nucleotides::GPUDeltaH(d1, d2, gpuData)) - temp * Nucleotides::GPUDeltaS(d1, d2, gpuData)) /
                        (boltzmann_constant * temp)
                    )
                )
            );
        }

        static inline float GPUDeltaS(int d1, int d2, int* gpuData) {
            return 0.023 * std::min(gpuData[GET_GPU_COUNT_COL(d1, GPU_NUCLEOTIDES)], gpuData[GET_GPU_COUNT_COL(d2, GPU_NUCLEOTIDES)]);
        }

        static inline float GPUDeltaH(int d1, int d2, int* gpuData) {
            float sum = 0;
            auto  len = std::min(gpuData[GET_GPU_COUNT_COL(d1, GPU_NUCLEOTIDES)], gpuData[GET_GPU_COUNT_COL(d2, GPU_NUCLEOTIDES)]);

            #pragma acc loop seq
            for (unsigned i = 0; i < len; i++) {
                sum += Nucleotides::bindPower(gpuData[GPU_ELEMENT(i, d1, GPU_NUCLEOTIDES)], gpuData[GPU_ELEMENT(i, d2, GPU_NUCLEOTIDES)]);
            }

            return sum;
        }

    private:
        std::vector<NUCLEOTIDE_DT>* strand;

};
