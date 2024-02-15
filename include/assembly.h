#pragma once
#include "molecule.h"

namespace assembly {
    void parse(const char* filename, std::vector<Molecule*> &registers, std::vector<std::string> &dictionary);
    Molecule* parseMolecule(std::string strMol, std::vector<std::string> &dictionary);
}