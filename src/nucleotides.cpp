#include <cstddef>
#include "nucleotides.h"
#include <iostream>

Nucleotides::Nucleotides() {
    this->strand = new std::vector<NUCLEOTIDE_DT>;
}

Nucleotides::~Nucleotides() {
    delete this->strand;
}

void Nucleotides::add(NUCLEOTIDE_DT na) {
    this->strand->push_back(na);
}

void Nucleotides::addFromStr(char na) {
    if        (na == 'A' || na == 'a') {
        return this->add(ADENINE);
    } else if (na == 'C' || na == 'c') {
        return this->add(CYTOSINE);
    } else if (na == 'G' || na == 'g') {
        return this->add(GUANINE);
    } else if (na == 'T' || na == 't' || na == 'U' || na == 'u') {
        return this->add(THYMINE);
    }
}

float Nucleotides::bindPower(NUCLEOTIDE_DT n1, NUCLEOTIDE_DT n2) {

    //Klump and Ackermann 1971 data
    if ((n1 == GUANINE && n2 == CYTOSINE) || (n2 == GUANINE && n1 == CYTOSINE)) {
        return -9.0;
    } else if ((n1 == ADENINE && n2 == THYMINE) || (n2 == ADENINE && n1 == THYMINE)) {
        return -7.2;
    } else {
        return -5.4;
    }

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