#pragma once

#include "molecule.h"

namespace output {
    void rawPrint(Molecule &molecule, char* space, std::vector<std::string> &names);
    void asciiPrint(Molecule &molecule, char* space, std::vector<std::string> &names);
    void dummyPrint(Molecule &molecule, char* space, std::vector<std::string> &names);
    void assemblyPrint(Molecule &molecule, char* space, std::vector<std::string> &names, std::vector<std::pair<std::string, std::string>> &macros);
    void svgPrint(Molecule &molecule, std::vector<std::string> &names);

    std::string getName(DOMAIN_DT val, std::vector<std::string> &names);
    std::string getNameShort(DOMAIN_DT val, std::vector<std::string> &names);
    bool isFirstBinded(Atom* a);
    bool isLastBinded(Atom* a);
    bool isLast(Atom* a);
    bool isFirst(Atom *a);
}