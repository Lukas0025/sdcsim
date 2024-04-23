/**
 * SDCSIM - strand displacement simulator on SIMD||DNA architecture
 * @file nucleotides.h
 * @brief Contain headers for nucleotides class
 * @author Lukáš Plevač <xpleva07@vut.cz>
 */

#pragma once

#include <vector>
#include <string>
#include <cstdint>

#define NUCLEOTIDE_DT uint8_t

// purines
#define ADENINE   0b01
#define GUANINE   0b10

// pyrmides
#define CYTOSINE  0xFD
#define THYMINE   0xFE

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
         * @param na1 nucleotide 
         * @param na2 nucleotide
         * @return energy in kcal
         */
        static float bindPower(NUCLEOTIDE_DT na1, NUCLEOTIDE_DT na2);

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

    private:
        std::vector<NUCLEOTIDE_DT>* strand;

};
