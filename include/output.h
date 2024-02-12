#pragma once

#include "molecule.h"

namespace output {
    void rawPrint(Molecule &molecule, char* space, std::vector<std::string> &names);
    void asciiPrint(Molecule &molecule, char* space, std::vector<std::string> &names);
}