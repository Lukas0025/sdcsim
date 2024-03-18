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

class Nucleotides {
    public:
        Nucleotides();
        ~Nucleotides();

        void addFromStr(char na);
        void add(NUCLEOTIDE_DT na);

        static float bindPower(NUCLEOTIDE_DT na1, NUCLEOTIDE_DT na2);

        std::string getStr();

    private:
        std::vector<NUCLEOTIDE_DT>* strand;

};