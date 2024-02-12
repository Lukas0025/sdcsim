#include <vector>
#include <string>
#include <iostream>
#include "output.h"

#define DEBUG

namespace output {
    
    unsigned                   maxNameCache    = 0;
    std::vector<std::string>*  maxNameCachePtr = NULL;

    inline unsigned getMaxName(std::vector<std::string> &names) {
        if (maxNameCachePtr != &names) {
            maxNameCache = 1;

            for(std::string & name : names) {
                maxNameCache = maxNameCache < name.length() ? name.length() : maxNameCache;
            } 
        }

        return maxNameCache;
    }

    inline std::string emptyDomain(std::vector<std::string> &names) {
        std::string str = "";

        // +1 for complementary char
        for (unsigned i = 0; i < getMaxName(names) + 1; i++) {
            str += " ";
        }

        return str;

    }

    inline std::string getName(DOMAIN_DT val, std::vector<std::string> &names) {

        const char* complement = IS_COMPLEMENTARY(val) ? "*" : " ";

        if (IS_COMPLEMENTARY(val)) val = ~val;

        #ifdef DEBUG
            if (names.size() <= val) return "?" + std::string(complement);  
        #endif

        std::string padName = names.at(val);

        for (unsigned i = padName.length(); i < getMaxName(names); i++) {
            padName += " ";
        }

        return padName + std::string(complement);
    }

    std::string rawSubPrint(Strand* strand, char* space, std::vector<std::string> &names, unsigned bindedAt, unsigned &negativePad, int &endAt) {
        std::string str = std::string("");
        int link = -1;
        
        for (unsigned i = 0; i < strand->length(); i++) {
            auto atom = strand->getAtom(i);

            str += getName(atom->domain.get(), names) + space;

            if (atom->partner != NULL && link < 0) { // I found link to main
                link = i;
            }
        }

        endAt = strand->length() - link;
        link  = bindedAt - link;

        // pad by link
        for (int i = 0; i < link; i++) {
            str = space + str;
            str = emptyDomain(names) + str;
        }

        if (link < 0) negativePad = -link;

        return str;

    }

    void rawPrint(Molecule &molecule, char* space, std::vector<std::string> &names) {
        std::vector<std::string> displayString;

        auto mainStrand = molecule.getStrand(0);

        displayString.push_back(std::string("")); //main chain

        int endAt = 0;

        for (unsigned i = 0; i < mainStrand->length(); i++) {
            auto atom = mainStrand->getAtom(i);

            displayString[0] += getName(atom->domain.get(), names) + space;

            if (atom->partner != NULL && endAt <= 0) {
                unsigned pad = 0;

                displayString.push_back(rawSubPrint(atom->partner->strand, space, names, i, pad, endAt));

                //todo pad all existig and pad all feature
                for (int i = 0; i < pad; i++) {
                    displayString[0] = space + displayString[0];
                    displayString[0] = emptyDomain(names) + displayString[0];
                }
            }

            endAt--;
        }

        for (int i = displayString.size() - 1; i >= 0; i--) {
            std::cout << displayString[i].c_str();
            std::cout << '\n';
        }

    }

    inline bool isLast(Atom* a) {
        return a == a->strand->getAtom(a->strand->length() - 1);
    }

    void asciiPrint(Molecule &molecule, char* space, std::vector<std::string> &names) {
        std::vector<std::string> displayString;

        auto mainStrand = molecule.getStrand(0);

        displayString.push_back(std::string(" ") + space); //naming
        displayString.push_back(std::string("O") + space); //main strand
        displayString.push_back(std::string(" ") + space); //hydrogen bound
        displayString.push_back(std::string(" ") + space); //upper strand

        for (unsigned i = 0; i < mainStrand->length(); i++) {
            auto atom = mainStrand->getAtom(i);

            displayString[0] += getName(atom->domain.get(), names) + space;
            displayString[1] += "-";

            if (atom->partner != NULL) {
                displayString[2] += "|";
                displayString[3] += isLast(atom->partner) ? ">" : "-";

                for (unsigned i = 0; i < getMaxName(names); i++) {
                    displayString[2] += " ";
                    displayString[3] += " ";
                }
            } else {
                for (unsigned i = 0; i < getMaxName(names) + 1; i++) {
                    displayString[2] += " ";
                    displayString[3] += " ";
                }
            }

            for (unsigned i = 0; i < getMaxName(names); i++) {
                displayString[1] += " ";
            }

            displayString[1] += space;
            displayString[2] += space;
            displayString[3] += space;
        }

        for (int i = displayString.size() - 1; i >= 0; i--) {
            std::cout << displayString[i].c_str();
            std::cout << '\n';
        }

    }
}