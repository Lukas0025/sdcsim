/**
 * SDCSIM - strand displacement simulator on SIMD||DNA architecture
 * @file register.h
 * @brief Contain headers for register class
 * @author Lukáš Plevač <xpleva07@vut.cz>
 */

#pragma once

#include "molecule.h"

/**
 * class reprezenting registr in simd||dna architecture
 */
class Register {
    public:
        /**
         * Class constructor
         * @param init molecule reprezenting registr 
         */
        Register(Molecule *init);

        /**
         * Get molecule of registr
         * @return molecule of registr 
         */
        Molecule* get();

        /**
         * Apply instruction on registr
         * @param instruction vector of instruction molecules / strads to apply in parallel
         * @param time        time to expose registr with strads
         */
        void applyInstruction(std::vector<Molecule*> instruction, int time);

        /**
         * Enbale nucleotides level of simulation applay of istruction
         * @param nucleotides map from domain to nucleotides reprezentation
         * @param temp        temperature of solution
         * @param strands     quntity of istructions strads when applay istruction
         */
        void enableNucleotidesLevel(std::map<DOMAIN_DT, Nucleotides*> &nucleotides, float temp, int strands);

    private:
        /**
         * Do all possible bindings in molecule with other
         * @param mol molecule to bind
         */
        void doAllBinding(Molecule* mol);

        /**
         * Remove replaced strands in molecule with new
         * @param oldSize size of molecule before add instructions strands
         */
        void removeReplaced(unsigned oldSize);

        /**
         * Remove all unstable strands form molecule (1domain or less binding)
         * @param to id of strand to do this reaction
         */
        void removeUnstable(int to);

        /**
         * Remove all unbided strands from molecule when is unbinded with instruction
         * @param mol molecule with instruction
         */
        void removeUnbinded(Molecule* mol);

        /**
         * Unbind domain when multiple strands ocupating this domain
         */
        void unbindOfMultiple();

        /**
         * Free domain when not multiple strands ocupating this domain
         */
        void bindOfMultiple();
        
        /**
         * Do elution (remove of unstable) on nucleotides level
         * @param nucleotides map from domain to nucleotides reprezentation
         */
        void elution(std::map<DOMAIN_DT, Nucleotides*> &nucleotides);

        Molecule *reg; // molecule reprezenting register

        bool nucleotidesSim;
        float temp;
        int strands_count;
        std::map<DOMAIN_DT, Nucleotides*>* nucleotides;
};