/**
 * SDCSIM - strand displacement simulator on SIMD||DNA architecture
 * @file svg.h
 * @brief Contain headers for svg class
 * @author Lukáš Plevač <xpleva07@vut.cz>
 */

#pragma once

#include <string>
#include <map>
#include "molecule.h"

// Supported color modes
#define CM_BLACK  0
#define CM_DOMAIN 1
#define CM_CHAIN  2

extern const char* svgColors[];         //< @brief supported svg colors
extern const unsigned svgColorsCount;   //< @brief supported svg colors count

/**
 * @brief class for svg generator 
 */
class Svg {
    public:
        /**
         * Class constructor
         * @param domainSize   size of one domain in pixels
         * @param lineWidth    width of strand in pixles
         * @param domainHSpace horizontal space between to binded domains in pixels
         * @param colorStrans  coloring mode from CM_BLACK, CM_DOMAIN, CM_CHAIN
         * @param names        dictionary of domains names
         */
        Svg(int domainSize, int lineWidth, int domainHSpace, uint8_t colorStrans, std::vector<std::string> &names);
        ~Svg();

        /**
         * Draw line in output svg
         * @param x1        x of fisrt point on line
         * @param x2        x of last point on line
         * @param y1        y of first point on line
         * @param y2        y of last point on line
         * @param color     color of line
         * @param width     width of line in pixels
         * @param dashed    is line dashed
         */
        void        line(int x1, int x2, int y1, int y2, const char* color, int width, bool dashed = false);

        /**
         * Draw text in output svg
         * @param x     x of center point of text 
         * @param y     y of center point of text 
         * @param text  text to draw
         * @param color color of text
         * @param transform transformating of text unsing svg transform
         */
        void        text(int x, int y, int size, const char* text, const char* color = "", const char* transform = "");

        /**
         * Draw domain in upper strand
         * @param pos position of domain by main bottom strand
         * @param domain domain to draw
         */
        void        upper(int pos, DOMAIN_DT domain);

        /**
         * Draw domain in bottom strand
         * @param pos position of domain by main bottom strand
         * @param domain domain to draw
         */
        void        bottom(int pos, DOMAIN_DT domain);

        /**
         * Draw overhang domain in upper strand
         * @param pos position of domain by main bottom strand
         * @param domain domain to draw
         */
        void        overHang(int pos, int bindDistance, DOMAIN_DT domain);

        /**
         * Draw domain what ending strand in upper strand
         * @param pos position of domain by main bottom strand
         * @param domain domain to draw
         */
        void        upperEnd(int pos, DOMAIN_DT domain);

        /**
         * Get svg coude of image
         * @return str svg
         */
        std::string get();

        /**
         * instruction is done jump to new line and name last istruction
         * @param name name of done instruction
         */
        void        inscructionDone(const char* name);

        /**
         * One register is done jump to new line 
         */
        void        regDone();

        /**
         * Set strand as currently in draw (for correct colors)
         * @param strand pointer for current strand
         */
        void setStrand(Strand *strand);

    private:
        /**
         * Get color for domain
         * @return string color
         */
        const char* getColorForDomain(DOMAIN_DT domain);

        /**
         * Last assigned color 
         */
        int lastColor;

        std::map<DOMAIN_DT, const char *> colors; //< @brief map of domains to colors
        std::string code;                         //< @brief curently generated svg code
        std::string reg;                          //< @brief curently generated register svg code
        std::string inscruction;                  //< @brief curently generated instruction svg code
        
        int domainSize;
        int lineWidth;
        int domainHSpace;

        int xMax;                                  //< @brief curently maximal X
        int yMax;                                  //< @brief curently maximal Y
        int yOffset;                               //< @brief curently offset in Y
        int absY;
        int absRegY;
        int firstReg;

        bool shortedNextUpper;
        
        uint8_t colorMode;

        Strand* currentStrand = NULL;             //< @brief currently strand in draw

        std::vector<std::string> *names;
};