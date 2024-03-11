#include "svg.h"
#include "output.h"

const char* svgColors[] {
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
    "#ffc107",
    "#ff9800",
    "#ff5722",
    "#e91e63",
    "#795548",
    "#9c27b0",
    "#607d8b"
};

const unsigned svgColorsCount = 18;

Svg::Svg(int domainSize, int lineWidth, int domainHSpace, std::vector<std::string> &names) {
    this->code = std::string("");
    this->domainSize   = domainSize;
    this->lineWidth    = lineWidth;
    this->domainHSpace = domainHSpace;
    this->xMax         = 0;
    this->yMax         = 0;
    this->lastColor    = 0;
    this->yOffset      = lineWidth + 14;
    this->names        = &names;
}

Svg::~Svg() {

}


void Svg::line(int x1, int x2, int y1, int y2, const char* color, int width, bool dashed) {
    this->code += "<line x1=\"" + std::to_string(x1) + 
                     "\" y1=\"" + std::to_string(y1 - this->yOffset) +
                     "\" x2=\"" + std::to_string(x2) +
                     "\" y2=\"" + std::to_string(y2 - this->yOffset) +
                     "\" style=\"stroke:" + color + ";stroke-width:" + std::to_string(width) + 
                     (dashed ? "\" stroke-dasharray=\"6,6\" />" : "\" />");

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

    auto color = this->getColorForDomain(domain);

    this->line(x1, x2, -this->domainHSpace, -this->domainHSpace, color, this->lineWidth);

    //connect with downer
    this->line(x1 + this->domainSize / 2, x1 + this->domainSize / 2, 0, -this->domainHSpace, color, this->lineWidth, true);
}

void Svg::upperEnd(int pos, DOMAIN_DT domain) {
    auto x2 = pos * this->domainSize + this->domainSize;

    auto color = this->getColorForDomain(domain);

    this->line(x2 - 8, x2 - 2, -this->domainHSpace + 8, -this->domainHSpace, color, this->lineWidth);
    this->line(x2 - 8, x2 - 2, -this->domainHSpace - 8, -this->domainHSpace, color, this->lineWidth);
}

void Svg::bottom(int pos, DOMAIN_DT domain) {
    auto x1 = pos * this->domainSize;
    auto x2 = pos * this->domainSize + this->domainSize;

    auto color = this->getColorForDomain(domain);

    this->line(x1, x2, 0, 0, color, this->lineWidth);

    //write text
    this->text(x1 + this->domainSize / 2, 14, 12, output::getNameShort(domain, *(this->names)).c_str(), color);
}

void Svg::overHang(int pos, int bindDistance, DOMAIN_DT domain) {

}

const char* Svg::getColorForDomain(DOMAIN_DT domain) {
    // we not need information about complementary
    domain = IS_COMPLEMENTARY(domain) ? ~domain : domain;

    if (!this->colors.count(domain)) {
        this->colors.insert({domain, svgColors[this->lastColor % svgColorsCount]});
        this->lastColor++;
    } 

    return this->colors[domain];
}

std::string Svg::get() {
    return "<svg viewBox=\"0 " + std::to_string(this->yMax) + " " + std::to_string(this->xMax) + " " + std::to_string(std::abs(this->yMax)) + 
             "\" height=\"" + std::to_string(std::abs(this->yMax)) + "\" width=\"" + std::to_string(this->xMax) + 
             "\" xmlns=\"http://www.w3.org/2000/svg\">" + this->code + "</svg>";
}