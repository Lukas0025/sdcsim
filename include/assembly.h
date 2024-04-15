#pragma once

#include "register.h"
#include "molecule.h"
#include "nucleotides.h"

namespace assembly {
    void parse(const char* filename, std::vector<Register*> &registers, std::vector<std::vector<Molecule*>> &inst, std::vector<std::string> &dictionary, std::vector<std::pair<std::string, std::string>> &macros, std::map<DOMAIN_DT, Nucleotides*> &nucleotides);
    Molecule* parseMolecule(std::string strMol, std::vector<std::string> &dictionary, std::map<DOMAIN_DT, Nucleotides*> *nucleotides = NULL, int line = 0, const char* linepos = "");
    std::string createAssembly(Molecule* mol, std::vector<std::string> &names);
    void replaceAll(std::string& str, const std::string& from, const std::string& to);
}