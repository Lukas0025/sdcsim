#pragma once

#include "register.h"
#include "molecule.h"

namespace assembly {
    void parse(const char* filename, std::vector<Register*> &registers, std::vector<std::vector<Molecule*>> &inst, std::vector<std::string> &dictionary);
    Molecule* parseMolecule(std::string strMol, std::vector<std::string> &dictionary);
}