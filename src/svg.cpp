#include "svg.h"
#include "output.h"

const char* svgColors[] {
    "red",
    "green",
    "blue",
    "orange",
    "purple",
    "cyan",
    "brown",
    "magenta",
    "dodgerblue",
    "orchid",
    "lawngreen",
    "lightskyblue",
    "darksalmon",
    "teal",
    "rosybrown",
    "#00bcd4",
    "#2196f3",
    "#009688",
    "#3f51b5",
    "#4caf50",
    "#8bc34a",
    "#03a9f4",
    "#cddc39",
    "#673ab7",
    "#f44336",
    "#ffeb3b",
    "#ff9800",
    "#e91e63",
    "#795548",
    "#9c27b0",
    "#607d8b"
};

const unsigned svgColorsCount = 31;

Svg::Svg(int domainSize, int lineWidth, int domainHSpace, bool colorStrans, std::vector<std::string> &names) {
    this->code             = std::string("");
    this->domainSize       = domainSize;
    this->lineWidth        = lineWidth;
    this->domainHSpace     = domainHSpace;
    this->yOffset          = lineWidth + 14;
    this->names            = &names;
    this->shortedNextUpper = true;
    this->colorStrans      = colorStrans;
    this->absRegY          = 0;
    this->firstReg         = 0;
    this->absY             = 0;
    this->xMax             = 0;
    this->yMax             = 0;
    this->lastColor        = 0;
}

Svg::~Svg() {

}


void Svg::line(int x1, int x2, int y1, int y2, const char* color, int width, bool dashed) {
    this->code += "<line x1=\"" + std::to_string(x1) + 
                     "\" y1=\"" + std::to_string(y1 - this->yOffset) +
                     "\" x2=\"" + std::to_string(x2) +
                     "\" y2=\"" + std::to_string(y2 - this->yOffset) +
                     "\" style=\"stroke:" + color + ";stroke-width:" + std::to_string(width) + 
                     (dashed ? "\" stroke-dasharray=\"" + std::to_string(width) + "," + std::to_string(width) + "\" />" : "\" />");

    this->xMax = std::max(std::max(this->xMax, x1), x2);
    this->yMax = std::min(std::min(this->yMax, y1 - width - this->yOffset ), y2 - width - this->yOffset);
}

void Svg::text(int x, int y, int size, const char* text, const char* color, const char* transform) {
    this->code += "<text x=\"" + std::to_string(x) + 
                     "\" y=\"" + std::to_string(y - this->yOffset) +
                     "\" fill=\"" + color + 
                     "\" font-size=\"" + std::to_string(size) +
                     "\" transform=\"" + transform + "\" text-anchor=\"middle\">" + text + "</text>";

    this->xMax = std::max(this->xMax, x + size);
    this->yMax = std::min(this->yMax, y - size - this->yOffset);
}

void Svg::upper(int pos, DOMAIN_DT domain) {
    auto x1 = pos * this->domainSize;
    auto x2 = pos * this->domainSize + this->domainSize;

    if (this->colorStrans && this->shortedNextUpper) this->lastColor++;

    auto color = this->getColorForDomain(domain);

    this->line(this->shortedNextUpper ? x1 + (this->domainSize - this->lineWidth) / 2 : x1, x2, -this->domainHSpace, -this->domainHSpace, color, this->lineWidth);

    //connect with downer
    this->line(x1 + this->domainSize / 2, x1 + this->domainSize / 2, 0, -this->domainHSpace, color, this->lineWidth, true);

    //write text
    //this->text(x1 + this->domainSize / 2, -this->domainHSpace - 10, 8, output::getNameShort(domain, *(this->names)).c_str(), color);

    this->shortedNextUpper = false;
}

void Svg::upperEnd(int pos, DOMAIN_DT domain) {
    auto x2 = pos * this->domainSize + this->domainSize;

    auto color = this->getColorForDomain(domain);

    //this->line(x2 - 7, x2 - 1, -this->domainHSpace + 4, -this->domainHSpace, color, this->lineWidth);
    this->line(x2 - this->domainSize / 2, x2 - 1, -this->domainHSpace - 5, -this->domainHSpace, color, this->lineWidth);

    this->shortedNextUpper = true;

    if (this->colorStrans) {
        this->lastColor++;
    }
}

void Svg::bottom(int pos, DOMAIN_DT domain) {
    auto x1 = pos * this->domainSize;
    auto x2 = pos * this->domainSize + this->domainSize;

    auto color = "black";
    if (!this->colorStrans) {
        color = this->getColorForDomain(domain);
    }

    this->line(x1, x2, 0, 0, color, this->lineWidth);

    //write text
    this->text(x1 + this->domainSize / 2, 14, 8, output::getNameShort(domain, *(this->names)).c_str(), color);
}

void Svg::overHang(int pos, int bindDistance, DOMAIN_DT domain) {
    if (bindDistance == -1) {
        this->shortedNextUpper = true;
    } else if (bindDistance == 1) {
        if (this->colorStrans && this->shortedNextUpper) this->lastColor++;
        this->shortedNextUpper = false;
    }

    bool right = false;
    if (bindDistance < 0) {
       bindDistance = -bindDistance;
       right = true;
    }

    auto x1 = pos * this->domainSize;
    auto x2 = x1 + this->domainSize;
    auto y1 = bindDistance * -this->domainHSpace;
    auto y2 = y1 - this->domainHSpace;

    if (!right) std::swap(y1, y2);

    auto color = this->getColorForDomain(domain);

    this->line(x1, x2, y1, y2, color, this->lineWidth);

    //connect with downer
    //this->line(x1 + this->domainSize / 2, x1 + this->domainSize / 2, 0, -this->domainHSpace, color, this->lineWidth, true);
    //write text
    //this->text(x1 + this->domainSize / 2, -this->domainHSpace - 10, 8, output::getNameShort(domain, *(this->names)).c_str(>

}

const char* Svg::getColorForDomain(DOMAIN_DT domain) {
    // we not need information about complementary
    domain = IS_COMPLEMENTARY(domain) ? ~domain : domain;

    if (this->colorStrans) {
        return svgColors[this->lastColor % svgColorsCount];
    }

    if (!this->colors.count(domain)) {
        this->colors.insert({domain, svgColors[this->lastColor % svgColorsCount]});
        this->lastColor++;
    }

    return this->colors[domain];
}

void Svg::inscructionDone(const char* name) {
    this->inscruction += "<text x=\"0\" y=\"" + std::to_string(this->absY + this->firstReg / 2) +
                         "\" fill=\"black\" font-size=\"15\" text-anchor=\"left\">" + name + "</text><g transform=\"translate(100," + std::to_string(this->absY) + ")\">" + this->reg + "</g>";
    this->reg          = "";
    this->absY        += this->absRegY;
    this->absRegY      = 0;
    this->firstReg     = 0;
}

void Svg::regDone() {
    //write offset
    this->reg     +=  "<g transform=\"translate(0," + std::to_string(-this->yMax + this->absRegY) + ")\">" + this->code + "</g>";
    this->absRegY += -this->yMax;
    this->firstReg = this->firstReg == 0 ? -this->yMax : this->firstReg;
    this->yMax     =  0;
    this->code     = "";
}

std::string Svg::get() {
    return "<svg height=\"" + std::to_string(this->absY) + "\" width=\"" + std::to_string(this->xMax + 100) + 
             "\" xmlns=\"http://www.w3.org/2000/svg\">" + this->inscruction + "</svg>";
}
