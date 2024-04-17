#include <iostream>

namespace error {
    void nucleotideLevelNotDefine() {
        std::cerr << "ERROR: Not all domains are defined at the nucleotide level\n";
        exit(1);
    }

    void assemblyUnknowControlSeq(const char *pos, int line) {
        std::cerr << "ERROR: Undefined control sequence in assebly file on line ";
        std::cerr << line;
        std::cerr << pos;
        exit(1);
    }

    void assemblyInvalidConcat(const char *pos, int line) {
        std::cerr << "ERROR: Invalid use of concatenation operation in assebly file on line ";
        std::cerr << line;
        std::cerr << pos;
        exit(1);
    }

    void assemblyUncompleteControlSeq(const char *pos, int line) {
        std::cerr << "ERROR: Uncomplete or unclosed control sequence in assebly file on line ";
        std::cerr << line;
        std::cerr << pos;
        exit(1);
    }

    void assemblyComplementOfNothing(const char *pos, int line) {
        std::cerr << "ERROR: Try to compute complement of nothing in assebly file on line ";
        std::cerr << line;
        std::cerr << pos;
        exit(1);
    }

    void unknowNucleotide(char* nuc) {
        std::cerr << "ERROR: unknow nucleotide name ";
        std::cerr << nuc;
        std::cerr << " supported is only A,C,T,G,U\n";
        exit(1);
    }

    void assemblyUnknowControlLabel(const char *label) {
        std::cerr << "ERROR: Undefined control label ";
        std::cerr << label;
        std::cerr << " in assebly file\n";
        exit(1);
    }

    void assemblyMultipleControlLabel(const char *label) {
        std::cerr << "ERROR: Multiple definition of control label ";
        std::cerr << label;
        std::cerr << " in assebly file\n";
        exit(1);
    }

    void inAppError(const char *name) {
        std::cerr << "IN APP ERROR: ";
        std::cerr << name;
        std::cerr << "\n";
        exit(1);
    }
}