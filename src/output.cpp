#include <vector>
#include <string>
#include <iostream>
#include <cstring>

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

    inline bool isLast(Atom* a) {
        return a == a->strand->getAtom(a->strand->length() - 1);
    }

    inline bool isFirst(Atom *a) {
        return a->strand->getAtom(0) == a;
    }

    inline bool isFirstBinded(Atom* a) {
        return a->strand->getAtom(0) == a                                || 
               a[-1].partner == NULL                                     ||
               a[-1].partner->partner == NULL                            ||
               a[-1].partner->partner->strand  != a->strand              ||
               a[-1].partner->strand != a->partner->strand;
    }

    inline bool isLastBinded(Atom* a) {
        return a->strand->getAtom(a->strand->length() - 1) == a          ||
               a[1].partner == NULL                                      ||
               a[1].partner->partner == NULL                             ||
               a[1].partner->strand           != a->partner->strand      ||
               a[1].partner->partner->strand  != a->strand               ||
               a[1].partner->partner == NULL;
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

    std::string rawSubPrint(Strand* strand, char* space, std::vector<std::string> &names, unsigned bindedAt, unsigned &negativePad) {
        std::string str = std::string("");
        int link = -1;
        
        for (unsigned i = 0; i < strand->length(); i++) {
            auto atom = strand->getAtom(i);

            str += getName(atom->domain.get(), names) + space;

            if (atom->partner != NULL && link < 0) { // I found link to main
                link = i;
            }
        }

        link  = bindedAt - link;

        // pad by link
        for (int i = 0; i < link; i++) {
            str = space + str;
            str = emptyDomain(names) + str;
        }

        if (link < 0) negativePad = -link;

        return str;

    }

    inline std::string buildPadding(unsigned size) {
        std::string str = "";

        for (unsigned i = 0; i < size; i++) {
            str += " ";
        }

        return str;
    }

    void rawPrint(Molecule &molecule, char* space, std::vector<std::string> &names) {
        std::vector<std::string> displayString;
        std::vector<unsigned> paddingSizes;

        auto mainStrand = molecule.getStrand(0);

        displayString.push_back(std::string("")); //main chain
        paddingSizes .push_back(0);

        int endAt = 0;

        for (unsigned i = 0; i < mainStrand->length(); i++) {
            auto atom = mainStrand->getAtom(i);

            displayString[0] += getName(atom->domain.get(), names) + space;

            if (atom->partner != NULL && atom->partner != multiAtom && isFirstBinded(atom->partner)) {
                unsigned pad = 0;

                displayString.push_back(rawSubPrint(atom->partner->strand, space, names, i, pad));

                // compute different of curent padding with new
                int diffPad = pad - paddingSizes[0];

                // pad all existig and pad all feature
                if (diffPad > 0) {
                    for (auto &padding : paddingSizes) {
                        padding += diffPad;
                    }
                }
                
                paddingSizes.push_back(pad - diffPad);
            }
        }

        //                          MAX NAME        COMPLEMENT   SPACE STR
        const unsigned onePadSize = getMaxName(names) + 1 + std::strlen(space);

        for (int i = displayString.size() - 1; i >= 0; i--) {
            std::cout << buildPadding(onePadSize * paddingSizes[i]);
            std::cout << displayString[i].c_str();
            std::cout << '\n';
        }

    }

    std::string getPosInMol(Molecule &molecule, Atom* atom) {
        for (unsigned iS = 0; iS < molecule.size(); iS++) {
            auto strand = molecule.getStrand(iS);

            for (unsigned i = 0; i < strand->length(); i++) {
                if (strand->getAtom(i) == atom) return "S" + std::to_string(iS) + "A" + std::to_string(i);
            }
        }

        return atom == multiAtom ? std::string("MORE") : std::string("FREE");

    }

    void dummyPrint(Molecule &molecule, char* space, std::vector<std::string> &names) {

        for (unsigned iS = 0; iS < molecule.size(); iS++) {
            
            std::cout << "S" << iS << " : ";

            auto strand = molecule.getStrand(iS);

            for (unsigned i = 0; i < strand->length(); i++) {

                auto atom = strand->getAtom(i);

                std::cout << getName(atom->domain.get(), names);
                std::cout << " = " << getPosInMol(molecule, atom->partner);
                std::cout << " (" << atom->partnersCount << ")";
                std::cout << " | ";


            }

            std::cout << '\n';
        }
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

            if (atom->partner != NULL && atom->partner != multiAtom) {
                const std::string atomStr = isLast(atom->partner)                                   ? ">"  : 
                                            isLastBinded(atom->partner)                             ? "/"  : 
                                            isFirstBinded(atom->partner) && !isFirst(atom->partner) ? "\\" : "-";

                displayString[2] += "|" + buildPadding(getMaxName(names));
                displayString[3] += atomStr + buildPadding(getMaxName(names));

                const auto partner = atom->partner;
                const auto first   = partner->strand->getAtom(0);
                const auto last    = partner->strand->getAtom(partner->strand->length() - 1);

                //if first binded?
                //we need display pre strend
                // \ < THIS
                //  -
                //  |
                //  -
                if (isFirstBinded(partner)) {
                    auto     overhang = partner;
                    int      pos      = displayString[3].length() - 1 - getMaxName(names) - std::strlen(space);
                    unsigned height   = 3;

                    while(overhang != first) {
                        //iterate over over hang
                        overhang = &(overhang[-1]);

                        //calculate new pos
                        height += 1;
                        pos    -= getMaxName(names) + 1 + std::strlen(space);

                        // do correction
                        if (pos < 0) {
                            //ok pad all others expect as and get same extra space
                            //we need -pos extra space for us
                            for (auto &str : displayString) {
                                str = buildPadding(-pos) + str;
                            }

                            pos = 0;
                        }

                        //push new
                        if (height >= displayString.size()) {
                            displayString.push_back(std::string(""));
                        }

                        //pad with spaces
                        if (pos >= displayString[height].length()) {
                            displayString[height] += buildPadding(pos - displayString[height].length() + 1);
                        }

                        displayString[height][pos] = displayString[height][pos] == '/' ? 'X' : '\\';
                    }
                }

                //if last binded?
                //we need display post strend
                //  / < THIS
                // -
                // |
                // -
                if (isLastBinded(partner)) {
                    auto     overhang = partner;
                    // -1 for convert len to pos
                    int      pos      = displayString[3].length() - 1 - getMaxName(names) - std::strlen(space);
                    unsigned height   = 3;

                    while(overhang != last) {
                        //iterate over over hang
                        overhang = &(overhang[1]);

                        //calculate new pos
                        height += 1;
                        pos    += getMaxName(names) + 1 + std::strlen(space);

                        //push new
                        if (height >= displayString.size()) {
                            displayString.push_back(std::string(""));
                        }

                        //pad with spaces
                        if (pos >= displayString[height].length()) {
                            displayString[height] += buildPadding(pos - displayString[height].length() + 1);
                        }

                        displayString[height][pos] = displayString[height][pos] == '\\' ? 'X' : '/';
                    }
                }

            } else {
                displayString[2] += buildPadding(getMaxName(names) + 1);
                displayString[3] += buildPadding(getMaxName(names) + 1);
            }

            displayString[1] += buildPadding(getMaxName(names)) + space;
            displayString[2] += space;
            displayString[3] += space;
        }

        for (int i = displayString.size() - 1; i >= 0; i--) {
            std::cout << displayString[i].c_str();
            std::cout << '\n';
        }

    }
}