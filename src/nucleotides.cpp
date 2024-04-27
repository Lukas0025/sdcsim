/**
 * SDCSIM - strand displacement simulator on SIMD||DNA architecture
 * @file nucleotides.cpp
 * @brief Contain implementation of Nucleotides class
 * @author Lukáš Plevač <xpleva07@vut.cz>
 */

#include <cstddef>
#include "nucleotides.h"
#include <iostream>
#include <math.h>       /* exp */
#include "error.h"

Nucleotides::Nucleotides() {
    this->strand = new std::vector<NUCLEOTIDE_DT>;
}

Nucleotides::~Nucleotides() {
    delete this->strand;
}

void Nucleotides::add(NUCLEOTIDE_DT na) {
    this->strand->push_back(na);
}

unsigned Nucleotides::length() {
    return this->strand->size();
}

NUCLEOTIDE_DT Nucleotides::get(unsigned pos) {
    return this->strand->at(pos);
}

void Nucleotides::addFromStr(char na, bool complementary) {
    if        (na == 'A' || na == 'a') {
        return this->add(complementary ? ~ADENINE  : ADENINE);
    } else if (na == 'C' || na == 'c') {
        return this->add(complementary ? ~CYTOSINE : CYTOSINE);
    } else if (na == 'G' || na == 'g') {
        return this->add(complementary ? ~GUANINE  : GUANINE);
    } else if (na == 'T' || na == 't' || na == 'U' || na == 'u') {
        return this->add(complementary ? ~THYMINE  : THYMINE);
    } else {
        error::unknowNucleotide(&na);
    }
}

float Nucleotides::ReNaturationP(Nucleotides* partner, float temp) {
    temp += 273.15; // convert to kelvin
    
    return std::min(1., 
        exp(
            -(
                (this->deltaH(partner) + temp * this->deltaS(partner)) /
                (boltzmann_constant * temp)
            )
        )
    );
}

float Nucleotides::DeNaturationP(Nucleotides* partner, float temp) {
    temp += 273.15; // convert to kelvin
    
    return std::min(1., 
        exp(
            -(
                ((-this->deltaH(partner)) - temp * this->deltaS(partner)) /
                (boltzmann_constant * temp)
            )
        )
    );
}

float Nucleotides::deltaH(Nucleotides* partner) {
    float sum = 0;

    for (unsigned i = 0; i < std::min(this->length(), partner->length()); i++) {
        sum += Nucleotides::bindPower(this->get(i), partner->get(i));
    }

    return sum;
}

float Nucleotides::deltaS(Nucleotides* partner) {
    return 0.023 * std::min(this->length(), partner->length());
}

std::string Nucleotides::getStr() {
    std::string out = "";

    for (const auto& val : *(this->strand)) {
        if      (val == ADENINE)  out += 'A';
        else if (val == THYMINE)  out += 'T';
        else if (val == CYTOSINE) out += 'C';
        else if (val == GUANINE)  out += 'G';
    }

    return out;
}