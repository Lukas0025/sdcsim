/**
 * SDCSIM - strand displacement simulator on SIMD||DNA architecture
 * @file assembly.h
 * @brief Contain headers for assembly namespace
 * @author Lukáš Plevač <xpleva07@vut.cz>
 */

#pragma once

#include "register.h"
#include "molecule.h"
#include "nucleotides.h"

namespace assembly {
    /**
     * Parse SDCASM file and return structs from file using referece params
     * @param filename    path to file to parse
     * @param registers   reference to vector of registers
     * @param inst        reference to vector of vectors of pointer to Molecule as instructions
     * @param dictionary  reference to vector of strings reprezenting dictionary form string domain name to int reprezentation in registers and insts
     * @param macros      reference to vector of pairs of strings reprezenting mecros as name -> value
     * @param nucleotides reference to map from DOMAIN_DT to Nucleotides* reprezenting map of Nucleotides reprezentation for domains
     */
    void parse(const char* filename, std::vector<Register*> &registers, std::vector<std::vector<Molecule*>> &inst, std::vector<std::string> &dictionary, std::vector<std::pair<std::string, std::string>> &macros, std::map<DOMAIN_DT, Nucleotides*> &nucleotides);
    
    /**
     * Parse string reprezetation of molecule and return struct reprezentation of this molecule using reperence params
     * @param strMol      string reprezentation of molecule
     * @param dictionary  reference to vector of strings reprezenting dictionary form string domain name to int reprezentation in registers and insts (New is append)
     * @param nucleotides reference to map from DOMAIN_DT to Nucleotides* reprezenting map of Nucleotides reprezentation for domains (New is append)
     * @param line        int number of line of this string reprezentation in file use for error message when parse fail (Default 0)
     * @param linepos     string reprezenting brief of position in file use for error message when parse fail (Default "")
     * @return Molecule* of string reprezentation
     */
    Molecule* parseMolecule(std::string strMol, std::vector<std::string> &dictionary, std::map<DOMAIN_DT, Nucleotides*> *nucleotides = NULL, int line = 0, const char* linepos = "");
    
    /**
     * Create string reprezentation of molecule
     * @param mol         pointer to molecule to convert
     * @param names       reference to vecor of string reprezenting dictionary of domains
     * @return string reprezentation
     */
    std::string createAssembly(Molecule* mol, std::vector<std::string> &names);

    /**
     * Replace all strings in string with other string
     * @param str   reference to tring where we replacing (Used as output too)
     * @param from  reference to string used as search pattern
     * @param to    reference to replace pattern   
     */
    void replaceAll(std::string& str, const std::string& from, const std::string& to);
}