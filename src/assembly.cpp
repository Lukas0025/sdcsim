#include <string>
#include <fstream>
#include <streambuf>
#include "molecule.h"
#include "output.h"
#include "assembly.h"
#include "error.h"
#include <regex>
#include <iostream>

#define NOT_BINDED -1

namespace assembly {
    const std::regex whitespaceRegex  ("[\t\f ]+");
    const std::regex defineRegex      ("\n *define *: *");
    const std::regex dataRegex        ("\n *data *: *");
    const std::regex domainRegex      ("\n *domains *: *");
    const std::regex instructionsRegex("\n *instructions *: *");
    const std::regex controlLabelRegex(".*:");

    void replaceAll(std::string& str, const std::string& from, const std::string& to) {
        if(from.empty()) return;

        size_t start_pos = 0;

        while((start_pos = str.find(from, start_pos)) != std::string::npos) {
            str.replace(start_pos, from.length(), to);
            start_pos += to.length();
        }
    }

    std::string createAssembly(Molecule* mol, std::vector<std::string> &names) {
        std::string code = "";

        auto mainStrand = mol->getStrand(0);

        bool isLower = false;

        for (unsigned i = 0; i < mainStrand->length(); i++) {

            auto atom = mainStrand->getAtom(i);

            if (atom->partner != NULL && atom->partner != multiAtom) {

                if (isLower) code += "}";

                if (output::isFirstBinded(atom->partner)) {
                    if (!output::isFirst(atom->partner)) {
                        auto strand = atom->partner->strand;
                        
                        code += "<";
                        
                        for (unsigned j = 0; j < strand->length(); j++) {
                            if (strand->getAtom(j) == atom->partner) break;

                            code += output::getNameShort(strand->getAtom(j)->domain.get(), names);
                        }

                        code += ">.";
                    }

                    code += "[";

                }

                code += output::getNameShort(atom->domain.get(), names);

                if (output::isLastBinded(atom->partner)) {
                    
                    code += "]";

                    if (!output::isLast(atom->partner)) {
                        auto strand = atom->partner->strand;
                        
                        std::string tmp = "";
                        
                        for (int j = strand->length() - 1; j >= 0; j--) {
                            if (strand->getAtom(j) == atom->partner) break;

                            tmp = output::getNameShort(strand->getAtom(j)->domain.get(), names) + tmp;
                        }

                        code += ".<" + tmp + ">";
                    }
                }

                isLower = false;
                
            } else {

                if (!isLower)     code += "{";

                code += output::getNameShort(atom->domain.get(), names);

                isLower                 = true;
            }


        }

        if (isLower) code += "}";

        return code;
    }

    void parseDefine(std::string asmStr, std::vector<std::pair<std::string, std::string>> &macros) {
        std::smatch defines;
        
        std::regex_search(asmStr, defines, defineRegex);
        
        // for every define
        for (unsigned i = 0; i < defines.size(); ++i) {
            // jump after label
            auto pos = defines.position(i) + defines.length(i);

            std::string left  = "";
            std::string right = "";
            bool isComment    = false;
            uint part         = 0;
            uint line         = 0;

            //ok now read char by char
            for (unsigned i = pos; i < asmStr.length(); i++) {

                if (asmStr[i] == '\n') {
                    if (left.length() > 0 && right.length() > 0) {
                        macros.push_back(
                            std::make_pair(
                                left,
                                right    
                            )
                        );

                        left  = "";
                        right = "";
                    }

                    part = 0;
                    line += 1;
                    isComment = false;
                    continue;
                }

                if (isComment) continue;

                if (asmStr[i] == ':') break; //new section 

                if (asmStr[i] == '#') {
                    isComment = true;
                    continue;
                }

                if (asmStr[i] == ' ') {
                    if (left.length() > 0) part++;

                    continue;
                }

                if (part == 0) {
                    left += asmStr[i];
                    continue;
                }

                if (part == 1) {
                    right += asmStr[i];
                    continue;
                }

                if (part > 1 && !isComment) {
                    error::assemblyUnknowControlSeq(" after define\n", line);
                }
            }
        }

    }

    DOMAIN_DT getKey(std::string name, std::vector<std::string> &dictionary) {
        
        for (DOMAIN_DT key = 0; key < dictionary.size(); key++) {
            if (name == dictionary.at(key)) return key;
        }

        dictionary.push_back(name);

        return dictionary.size() - 1;
    }

    Molecule* parseMolecule(std::string strMol, std::vector<std::string> &dictionary, std::map<DOMAIN_DT, Nucleotides*> *nucleotides, int line, const char* linepos) {
        bool up             = false;
        bool dbl            = false;
        bool down           = false;
        bool concat         = false;
        bool complement     = false;

        auto mainStrend = new Strand(nucleotides);

        // first load only main strand
        for (unsigned i = 0; i < strMol.length(); i++) {

            if (strMol[i] == '[') { // start double strand
                dbl  = true;
            } else if (strMol[i] == '{') { // start of bottom strand
                down = true;
            } else if (strMol[i] == ']') { // end of double strand
                dbl  = false;
            } else if (strMol[i] == '}') { // end of bottom strand
                down = false;
            } else {
                if (down || dbl) {
                    if (strMol[i] == '*') {
                        mainStrend->complementaryLast();
                    } else {
                        mainStrend->addDomain(
                            Domain(
                                getKey(std::string(&strMol[i], 1), dictionary)
                            )
                        );
                    }
                }
            }
        }

        // second parse all binded strands

        up         = false;
        dbl        = false;
        down       = false;
        concat     = false;
        complement = false;

        Strand*  workStrand    = NULL;
        int      startPos      = 0;
        int      bindStart     = NOT_BINDED; // not have bind start
        unsigned bindEnd       = 0;
        unsigned mainStrandPos = 0;

        strMol += "{}"; // for parse lastest strand

        auto molecule = new Molecule(mainStrend);

        for (unsigned i = 0; i < strMol.length(); i++) {

            if (strMol[i] == '<') { // start of upper strand
                if (up || down || dbl) error::assemblyUnknowControlSeq(linepos, line);

                up   = true;

                if (!concat) {
                    if (workStrand != NULL) {
                        auto id = molecule->addStrand(workStrand, startPos);
                        molecule->donePairStrands(id, bindStart - startPos, bindEnd - startPos);
                    }

                    workStrand = new Strand(nucleotides);
                    startPos   = mainStrandPos;
                    bindStart  = NOT_BINDED;
                }

                concat = false;
            } else if (strMol[i] == '[') { // start double strand
                if (up || down || dbl) error::assemblyUnknowControlSeq(linepos, line);

                if (!concat) {
                    if (workStrand != NULL) {
                        auto id = molecule->addStrand(workStrand, startPos);
                        molecule->donePairStrands(id, bindStart - startPos, bindEnd - startPos);
                    }

                    workStrand = new Strand(nucleotides);
                    startPos   = mainStrandPos;
                }

                dbl       = true;
                bindStart = mainStrandPos;

                concat = false;
            } else if (strMol[i] == '{') { // start of bottom strand
                if (up || down || dbl) error::assemblyUnknowControlSeq(linepos, line);
                down = true;

                if (concat) {
                    error::assemblyInvalidConcat(linepos, line);
                }

                if (workStrand != NULL) {
                    auto id = molecule->addStrand(workStrand, startPos);
                    molecule->donePairStrands(id, bindStart - startPos, bindEnd - startPos);
                }

                workStrand = NULL;
                concat = false;
            } else if (strMol[i] == '>') { // end of upper strand
                if (!up) error::assemblyUnknowControlSeq(linepos, line);

                up     = false;
            } else if (strMol[i] == ']') { // end of double strand
                if (!dbl) error::assemblyUnknowControlSeq(linepos, line);

                bindEnd = mainStrandPos - 1;
                dbl    = false;
            } else if (strMol[i] == '}') { // end of bottom strand
                if (!down) error::assemblyUnknowControlSeq(linepos, line);

                down   = false;
            } else if (strMol[i] == '.') { // concatenation
                if (concat) error::assemblyUnknowControlSeq(linepos, line);

                concat = true;
            } else {            
                if (up == true) {
                    if (strMol[i] == '*') {
                        if (workStrand->length() == 0) error::assemblyComplementOfNothing(linepos, line);

                        workStrand->complementaryLast();
                    } else {
                        workStrand->addDomain(
                            Domain(
                                getKey(std::string(&strMol[i], 1), dictionary)
                            )
                        );

                        //binding not started
                        if (bindStart == NOT_BINDED) {
                            startPos--;
                        }
                    }
                } else if (down == true) {
                    if (strMol[i] != '*') mainStrandPos++;
                } else if (dbl  == true) {
                    if (strMol[i] == '*') {
                        if (workStrand->length() == 0) error::assemblyComplementOfNothing(linepos, line);
                        
                        workStrand->complementaryLast();
                    } else {
                        workStrand->addDomain(
                            Domain(
                                getKey(std::string(&strMol[i], 1), dictionary)
                            )
                        );

                        workStrand->complementaryLast();

                        mainStrandPos++;
                    }
                } else {
                    error::assemblyUnknowControlSeq(linepos, line);
                }
            }
        }

        if (up || down || dbl || concat) error::assemblyUncompleteControlSeq(linepos, line);

        return molecule;
    }

    void parseDomains(std::string asmStr, std::map<DOMAIN_DT, Nucleotides*> &nucleotides, std::vector<std::string> &dictionary) {
        std::smatch domains;
        
        std::regex_search(asmStr, domains, domainRegex);
        
        // for every domains
        for (unsigned i = 0; i < domains.size(); ++i) {
            // jump after label
            auto pos = domains.position(i) + domains.length(i);

            std::string left  = "";
            std::string right = "";
            bool isComment    = false;
            uint part         = 0;
            int line          = 0;

            //ok now read char by char
            for (unsigned i = pos; i < asmStr.length(); i++) {

                if (asmStr[i] == '\n') {
                    if (left.length() > 0) {
                        left  = "";
                        right = "";
                    }

                    part = 0;
                    isComment = false;
                    line += 1;
                    continue;
                }

                if (isComment) continue;

                if (asmStr[i] == ':') break; //new section 

                if (asmStr[i] == '#') {
                    isComment = true;
                    continue;
                }

                if (asmStr[i] == ' ') {
                    if (left.length() > 0) part++;

                    continue;
                }

                if (part == 0) {
                    left += asmStr[i];
                    continue;
                }

                if (part == 1) {
                    //find in dictionary
                    auto dk = getKey(left, dictionary);

                    if (!nucleotides.count(dk)) {
                        nucleotides[dk]  = new Nucleotides();
                        nucleotides[~dk] = new Nucleotides();
                    }

                    nucleotides[dk] ->addFromStr(asmStr[i]);
                    nucleotides[~dk]->addFromStr(asmStr[i], true);
                    
                    continue;
                }
            }
        }
    }

    void parseInstructions(std::string asmStr, std::vector<std::vector<Molecule*>> &instructions, std::vector<std::string> &dictionary, std::map<DOMAIN_DT, Nucleotides*> &nucleotides) {
        std::smatch datas;
        
        std::regex_search(asmStr, datas, instructionsRegex);
        
        // for every define
        for (unsigned i = 0; i < datas.size(); ++i) {
            // jump after label
            auto pos = datas.position(i) + datas.length(i);

            std::string reading  = "";
            bool isComment       = false;
            int line             = 0;

            std::vector<Molecule*> molecules;

            //ok now read char by char
            for (unsigned i = pos; i < asmStr.length(); i++) {

                if (asmStr[i] == '\n') {
                    if (reading.length() > 0) {
                        molecules.push_back(parseMolecule(reading, dictionary, &nucleotides, line, " after instructions\n"));
                    }

                    reading = "";

                    if (molecules.size() > 0) {
                        instructions.push_back(molecules);
                        molecules.clear();
                    }

                    isComment = false;
                    line++;
                    continue;
                }

                if (isComment) continue;

                if (asmStr[i] == '#') {
                    isComment = true;
                    continue;
                }

                if (asmStr[i] == ':') break; //new section 

                if (asmStr[i] == ' ') {

                    if (reading.length() > 0) {
                        molecules.push_back(parseMolecule(reading, dictionary, &nucleotides, line, " after instructions\n"));
                        reading = "";
                    }

                    continue;
                }
                
                reading += asmStr[i];
            }
        }

    }

    void parseData(std::string asmStr, std::vector<Register*> &registers, std::vector<std::string> &dictionary, std::map<DOMAIN_DT, Nucleotides*> &nucleotides) {
        std::smatch datas;
        
        std::regex_search(asmStr, datas, dataRegex);
        
        // for every define
        for (unsigned i = 0; i < datas.size(); ++i) {
            // jump after label
            auto pos = datas.position(i) + datas.length(i);

            std::string reading  = "";
            bool isComment       = false;
            int line             = 0;

            //ok now read char by char
            for (unsigned i = pos; i < asmStr.length(); i++) {

                if (asmStr[i] == '\n') {
                    if (reading.length() > 0) {
                        registers.push_back(new Register(parseMolecule(reading, dictionary, &nucleotides, line, " after data\n")));

                        reading = "";
                    }

                    line++;
                    isComment = false;
                    continue;
                }

                if (isComment) continue;

                if (asmStr[i] == '#') {
                    isComment = true;
                    continue;
                }

                if (asmStr[i] == ':') break; //new section 

                if (asmStr[i] == ' ') {
                    continue;
                }
                
                reading += asmStr[i];
            }
        }

    }

    void parse(const char* filename, std::vector<Register*> &registers, std::vector<std::vector<Molecule*>> &inst, std::vector<std::string> &dictionary, std::vector<std::pair<std::string, std::string>> &macros, std::map<DOMAIN_DT, Nucleotides*> &nucleotides) {
        //first open file
        std::ifstream asmFile(filename);
        
        //read whole file to string -- not good for big files.
        std::string asmStr((std::istreambuf_iterator<char>(asmFile)), std::istreambuf_iterator<char>());

        //rewrite every whitespace to space + redundart
        asmStr = std::regex_replace(asmStr, whitespaceRegex, " ");

        //add new line at start
        asmStr = '\n' + asmStr + '\n'; // for working regex

        //check all control labels
        std::smatch contolLabels;
        std::string::const_iterator text_iter = asmStr.cbegin();

        bool defineD = false;
        bool dataD   = false;
        bool instD   = false;
        bool domainD = false;

        while (std::regex_search(text_iter, asmStr.cend(), contolLabels, controlLabelRegex)) {
            bool pass = false;

            if (std::regex_match("\n" + contolLabels[0].str(), defineRegex)) {
                if (defineD) error::assemblyMultipleControlLabel(contolLabels[0].str().c_str());

                defineD = true;
                pass    = true;
            }

            if (std::regex_match("\n" + contolLabels[0].str(), dataRegex)) {
                if (dataD) error::assemblyMultipleControlLabel(contolLabels[0].str().c_str());

                dataD   = true;
                pass    = true;
            }
            
            if (std::regex_match("\n" + contolLabels[0].str(), instructionsRegex)) {
                if (instD) error::assemblyMultipleControlLabel(contolLabels[0].str().c_str());

                instD   = true;
                pass    = true;
            }
            
            if (std::regex_match("\n" + contolLabels[0].str(), domainRegex)) {
                if (domainD) error::assemblyMultipleControlLabel(contolLabels[0].str().c_str());

                domainD = true;
                pass    = true;
            }

            if (!pass) {
                error::assemblyUnknowControlLabel(contolLabels[0].str().c_str());
            }

            //shift itereator
            text_iter = contolLabels.suffix().first;
        }

        //get macros
        parseDefine(asmStr, macros);

        //replace left with right in string
        for (auto &macro : macros) {
            replaceAll(asmStr, macro.first, macro.second);
        }

        // load nucleodies
        parseDomains(asmStr, nucleotides, dictionary);

        // load registers
        parseData(asmStr, registers, dictionary, nucleotides);

        // load instructions
        parseInstructions(asmStr, inst, dictionary, nucleotides);
    }
}