/**
 * SDCSIM - strand displacement simulator on SIMD||DNA architecture
 * @file error.h
 * @brief Contain headers for error namespace
 * @author Lukáš Plevač <xpleva07@vut.cz>
 */

namespace error {
    /**
     * Error when some domain do not have nucleotides reprezentation
     */
    void nucleotideLevelNotDefine();

    /**
     * Error when some unknow control sequence in assembly file
     * @param pos brief of position
     * @param line line number of error
     */
    void assemblyUnknowControlSeq(const char *pos, int line);

    /**
     * Error when use concat operator in assembly in bad context
     * @param pos brief of position
     * @param line line number of error
     */
    void assemblyInvalidConcat(const char *pos, int line);

    /**
     * Error when some control sequence is uncoplete
     * @param pos brief of position
     * @param line line number of error
     */
    void assemblyUncompleteControlSeq(const char *pos, int line);

    /**
     * Error when try compute complement of nothing
     * @param pos brief of position
     * @param line line number of error
     */
    void assemblyComplementOfNothing(const char *pos, int line);

    /**
     * Error when use of unknow nucleotide name
     * @param nuc used unknow nucleotide
     */
    void unknowNucleotide(char* nuc);

    /**
     * Error when use of unknow control label (other that define, data, ...)
     * @param nuc used unknow nucleotide
     */
    void assemblyUnknowControlLabel(const char *label);

    /**
     * Error when use some control label more that once
     * @param label name of reused label
     */
    void assemblyMultipleControlLabel(const char *label);

    /**
     * Error unspecified in app error
     * @param name brief of error 
     */
    void inAppError(const char *name);

    /**
     * Error when try open assembly file
     * @param file path to file
     */
    void openAsmFile(const char *file);
}