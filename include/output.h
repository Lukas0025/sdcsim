/**
 * SDCSIM - strand displacement simulator on SIMD||DNA architecture
 * @file output.h
 * @brief Contain headers for output namespace
 * @author Lukáš Plevač <xpleva07@vut.cz>
 */

#pragma once

#include "molecule.h"
#include "svg.h"

namespace output {
    /**
     * Print molecule in raw reprezenation
     * @param molecule molecule to print
     * @param space space char(s) between domains
     * @param names dictionary to convert in app reprezentation to string
     */
    void rawPrint(Molecule &molecule, const char* space, std::vector<std::string> &names);

    /**
     * Print molecule in ascii reprezenation
     * @param molecule molecule to print
     * @param space space char(s) between domains
     * @param names dictionary to convert in app reprezentation to string
     */
    void asciiPrint(Molecule &molecule, const char* space, std::vector<std::string> &names);

    /**
     * Print molecule in dummy reprezenation
     * @param molecule molecule to print
     * @param space space char(s) between domains
     * @param names dictionary to convert in app reprezentation to string
     */
    void dummyPrint(Molecule &molecule, const char* space, std::vector<std::string> &names);

    /**
     * Print molecule in assembly reprezenation
     * @param molecule molecule to print
     * @param space space char(s) between domains
     * @param names dictionary to convert in app reprezentation to string
     * @param macros macros to reverse decode data
     */
    void assemblyPrint(Molecule &molecule, const char* space, std::vector<std::string> &names, std::vector<std::pair<std::string, std::string>> &macros);

    /**
     * Print molecule in svg reprezenation
     * @param molecule molecule to print
     * @param names dictionary to convert in app reprezentation to string
     * @param svg object to work with svg
     */
    void svgPrint(Molecule &molecule, std::vector<std::string> &names, Svg &svg);

    /**
     * Get name of DOMAIN_DT value in string
     * @param val value to convert to str reprezenation
     * @param names dictionary to convert in app reprezentation to string
     * @return string reprezentation with complement support
     */
    std::string getName(DOMAIN_DT val, std::vector<std::string> &names);

    /**
     * Get name of DOMAIN_DT value in string without extra spacing 
     * @param val value to convert to str reprezenation
     * @param names dictionary to convert in app reprezentation to string
     * @return string reprezentation with complement support
     */
    std::string getNameShort(DOMAIN_DT val, std::vector<std::string> &names);

    /**
     * Check if domain in strand is first binded with other
     * @param a poiter to atom to check
     * @return bool if is first or not
     */
    bool isFirstBinded(Atom* a);

    /**
     * Check if domain in strand is last binded with other
     * @param a poiter to atom to check
     * @return bool if is first or not
     */
    bool isLastBinded(Atom* a);
    
    /**
     * Check if domain in strand is last
     * @param a poiter to atom to check
     * @return bool if is first or not
     */
    bool isLast(Atom* a);

    /**
     * Check if domain in strand is first
     * @param a poiter to atom to check
     * @return bool if is first or not
     */
    bool isFirst(Atom *a);
}