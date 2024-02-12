#include <string>
#include <fstream>
#include <streambuf>
#include "molecule.h"
#include <regex>
#include <iostream>

namespace assembly {
    const std::regex whitespaceRegex  ("[\t\f ]+");
    const std::regex defineRegex      ("\n *define *: *");
    const std::regex dataRegex        ("\n *data *: *");
    const std::regex instructionsRegex("[\t\f\v ]*instructions[\t\f\v ]*:[\t\f\v ]*");

    void replaceAll(std::string& str, const std::string& from, const std::string& to) {
        if(from.empty()) return;

        size_t start_pos = 0;

        while((start_pos = str.find(from, start_pos)) != std::string::npos) {
            str.replace(start_pos, from.length(), to);
            start_pos += to.length();
        }
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

    Molecule* parseMolecule(std::string strMol, std::vector<std::string> &dictionary) {
        bool up             = false;
        bool dbl            = false;
        bool down           = false;
        bool concat         = false;
        bool complement     = false;
        int  mainStrandPos  = 0;

        auto mainStrend = new Strand();

        // first load only main strand
        for (unsigned i = 0; i < strMol.length(); i++) {

            if (strMol[i] == '<') { // start of upper strand
                up   = true;
            } else if (strMol[i] == '[') { // start double strand
                dbl  = true;
            } else if (strMol[i] == '{') { // start of bottom strand
                down = true;
            } else if (strMol[i] == ']') { // end of double strand
                dbl  = false;
            } else if (strMol[i] == '}') { // end of bottom strand
                down = false;
            } else if (strMol[i] == '>') { // end of upper strand
                up   = false;
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
                } else if (up) {
                    if (mainStrend->length() == 0 && strMol[i] != '*') mainStrandPos--;
                }
            }
        }

        up         = false;
        dbl        = false;
        down       = false;
        concat     = false;
        complement = false;

        Strand* workStrand     = NULL;
        unsigned startPos      = 0;
        unsigned bindStart     = 0;
        unsigned bindEnd       = 0;

        auto molecule = new Molecule(mainStrend);

        for (unsigned i = 0; i < strMol.length(); i++) {

            if (strMol[i] == '<') { // start of upper strand
                up   = true;

                if (!concat) {
                    if (workStrand != NULL) {
                        auto id = molecule->addStrand(workStrand, startPos);
                        molecule->donePairStrands(id, bindStart, bindEnd);
                    }

                    workStrand = new Strand();
                    startPos = mainStrandPos;
                }

                concat = false;
            } else if (strMol[i] == '[') { // start double strand
                if (!concat) {
                    if (workStrand != NULL) {
                        auto id = molecule->addStrand(workStrand, startPos);
                        molecule->donePairStrands(id, bindStart - startPos, bindEnd - startPos);
                    }

                    workStrand = new Strand();
                    startPos   = mainStrandPos;
                }

                dbl       = true;
                bindStart = mainStrandPos;

                concat = false;
            } else if (strMol[i] == '{') { // start of bottom strand
                down = true;

                if (workStrand != NULL) {
                    auto id = molecule->addStrand(workStrand, startPos);
                    molecule->donePairStrands(id, bindStart - startPos, bindEnd - startPos);
                }

                workStrand = NULL;

            } else if (strMol[i] == '>') { // end of upper strand
                up   = false;
            } else if (strMol[i] == ']') { // end of double strand
                bindEnd = mainStrandPos;
                dbl     = false;
            } else if (strMol[i] == '}') { // end of bottom strand
                down = false;
            } else if (strMol[i] == '.') { // concatenation
                concat = true;
            } else {
            
                if (up == true) {
                    if (strMol[i] == '*') {
                        workStrand->complementaryLast();
                    } else {
                        workStrand->addDomain(
                            Domain(
                                getKey(std::string(&strMol[i], 1), dictionary)
                            )
                        );

                        if (mainStrandPos < 0) mainStrandPos++;
                    }
                } else if (down == true && strMol[i] != '*') {
                    mainStrandPos++;
                } else if (dbl  == true) {
                    if (strMol[i] == '*') {
                        workStrand->complementaryLast();
                    } else {
                        workStrand->addDomain(
                            Domain(
                                getKey(std::string(&strMol[i], 1), dictionary)
                            )
                        );

                        workStrand->complementaryLast();
                    }
                    mainStrandPos++;
                } else {
                    //error here
                }
            }
        }

        //add lastest
        if (workStrand != NULL) {
            auto id = molecule->addStrand(workStrand, startPos);
            molecule->donePairStrands(id, bindStart - startPos, bindEnd - startPos);
        }

        return molecule;
    }

    void parseData(std::string asmStr, std::vector<Molecule*> &registers, std::vector<std::string> &dictionary) {
        std::smatch datas;
        
        std::regex_search(asmStr, datas, dataRegex);
        
        // for every define
        for (unsigned i = 0; i < datas.size(); ++i) {
            // jump after label
            auto pos = datas.position(i) + datas.length(i);

            std::string reading  = "";
            bool isComment       = false;

            //ok now read char by char
            for (unsigned i = pos; i < asmStr.length(); i++) {

                if (asmStr[i] == '\n') {
                    if (reading.length() > 0) {
                        registers.push_back(parseMolecule(reading, dictionary));

                        reading = "";
                    }

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

    void parse(char* filename, std::vector<Molecule*> &registers, std::vector<std::string> &dictionary) {
        //first open file
        std::ifstream asmFile(filename);
        
        //read whole file to string -- not good for big files.
        std::string asmStr((std::istreambuf_iterator<char>(asmFile)), std::istreambuf_iterator<char>());

        //remove redundat spaces
        asmStr = std::regex_replace(asmStr, whitespaceRegex, " ");

        std::vector<std::pair<std::string, std::string>> macros;

        //get macros
        parseDefine(asmStr, macros);

        //replace left with right in string
        for (auto &macro : macros) {
            replaceAll(asmStr, macro.first, macro.second);
        }

        parseData(asmStr, registers, dictionary);
    }
}