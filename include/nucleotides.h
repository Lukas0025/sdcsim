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

        void addFromStr(char na, bool complementary = false);
        void add(NUCLEOTIDE_DT na);

        unsigned length();
        NUCLEOTIDE_DT get(unsigned pos);

        static float bindPower(NUCLEOTIDE_DT na1, NUCLEOTIDE_DT na2);

        std::string getStr();

        Nucleotides operator~() {
            auto nuc = Nucleotides();

            for (const auto& val : *(this->strand)) {
                nuc.add(~val);
            }

            return nuc;
        }

        float ReNaturationP(Nucleotides* partner, float temp);
        float DeNaturationP(Nucleotides* partner, float temp);
        float deltaH(Nucleotides* partner);
        float deltaS(Nucleotides* partner);

    private:
        std::vector<NUCLEOTIDE_DT>* strand;

};
