#pragma once

#include <string>
#include <map>
#include "molecule.h"

extern const char* svgColors[];
extern const unsigned svgColorsCount;

class Svg {
    public:
        Svg(int domainSize, int lineWidth, int domainHSpace, bool colorStrans, std::vector<std::string> &names);
        ~Svg();
        void        line(int x1, int x2, int y1, int y2, const char* color, int width, bool dashed = false);
        void        text(int x, int y, int size, const char* text, const char* color = "", const char* transform = "");
        void        upper(int pos, DOMAIN_DT domain);
        void        bottom(int pos, DOMAIN_DT domain);
        void        overHang(int pos, int bindDistance, DOMAIN_DT domain);
        void        upperEnd(int pos, DOMAIN_DT domain);
        std::string get();
        void        inscructionDone(const char* name);
        void        regDone();

    private:
        const char* getColorForDomain(DOMAIN_DT domain);
        int lastColor;

        std::map<DOMAIN_DT, const char *> colors;
        std::string code;
        std::string reg;
        std::string inscruction;
        
        int domainSize;
        int lineWidth;
        int domainHSpace;

        int xMax;
        int yMax;
        int yOffset;
        int absY;
        int absRegY;
        int firstReg;

        bool shortedNextUpper;
        bool colorStrans;
        bool colorDomains;

        std::vector<std::string> *names;
};