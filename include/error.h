namespace error {
    void nucleotideLevelNotDefine();
    void assemblyUnknowControlSeq(const char *pos, int line);
    void assemblyInvalidConcat(const char *pos, int line);
    void assemblyUncompleteControlSeq(const char *pos, int line);
    void assemblyComplementOfNothing(const char *pos, int line);
    void unknowNucleotide(char* nuc);
    void assemblyUnknowControlLabel(const char *label);
    void assemblyMultipleControlLabel(const char *label);
}