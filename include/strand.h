/**
 * SDCSIM - strand displacement simulator on SIMD||DNA architecture
 * @file strand.h
 * @brief Contain headers for strand class
 * @author Lukáš Plevač <xpleva07@vut.cz>
 */

#pragma once

class Strand;

#include <vector>
#include "domain.h"
#include "nucleotides.h"
#include <map>
#include <cstddef>

#define HAVE_PARTNER(X) (X->partner != NULL && X->partner != multiAtom)

/**
 * @brief stract reprezenting one atom in strand atom is structore with domain
 * poiter to partener and other support values  
 */
typedef struct atom_ {
    Domain         domain;        //< @brief domain of atom
    struct atom_*  partner;       //< @brief pointer to my partner
    Strand*        strand;        //< @brief pointer to my strand
    int            partnersCount; //< @brief number of ocupating strands (atoms)
    unsigned       offsetN;       //< @brief nucleotide offset in strand
} Atom;

extern Atom* multiAtom;           //< @brief atom to use when multiple strand ocupating one atom

/**
 * @brief Class reprezenting DNA strand 
 */
class Strand {
    public:
        /**
         * Class constructor
         * @param nucleotides map of nucleotides reprezentation
         */
        Strand(std::map<DOMAIN_DT, Nucleotides*> *nucleotides = NULL);
        ~Strand();

        /**
         * Deep copy current strand
         * @return pointer to new strand 
         */
        Strand* copy();

        /**
         * Add domian to current strand (it create new atom)
         * @param domian domain to add 
         */
        unsigned addDomain(Domain domain);

        /**
         * Full duplex bind to atoms
         * @param atom1 atom to pair
         * @param atom2 atom to pair
         */
        static void pairDomain(Atom* atom1, Atom* atom2);

        /**
         * Half duplex bind to atoms
         * @param atom1 atom to pair
         * @param atom2 atom to pair
         */
        static void halfPairDomain(Atom* atom1, Atom* atom2);

        /**
         * Do unpair of partner atom
         * @param atom1 atom to unpair with partner
         */
        static void unpairDomain(Atom* atom1);

        /**
         * Get atom on position
         * @param index index of atom
         * @return pointer to atom
         * @post method addDomin is now not supported
         */
        Atom* getAtom(unsigned index);

        /**
         * Get length of strand 
         * @return number of atoms in strand
         */
        unsigned length();

        /**
         * Set last atom as complementary 
         */
        void complementaryLast();

        /**
         * Mark current strand as deleted
         */
        void del();

        /**
         * check if delete mark is set
         * @return bool is deleted 
         */
        bool isDeleted();

        /**
         * get color of strand
         * @return char* used color
         */
        const char* getColor();

        /**
         * set color of strand
         * @param color char* of color
         */
        void setColor(const char* color);

    private:
        /**
         * Create new atom from domain
         * @param domain domain to use in atom
         * @return created atom
         */
        Atom createAtom(Domain domain);

        std::vector<Atom>* atoms; //< @brief strand as array of atoms
        
        bool readOnly;            //< @brief readOnly mark
        bool deleted;             //< @brief deleted mark

        std::map<DOMAIN_DT, Nucleotides*> *nucleotides;

        const char* color;              //< @brief color string for output namaspace

};
