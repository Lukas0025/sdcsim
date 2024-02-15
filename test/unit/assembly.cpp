#include <gtest/gtest.h>
#include "assembly.h"

#define  getStrandAtom(MOL, STRAND, ATOM)   (MOL->getStrand(STRAND)->getAtom(ATOM))
#define  getStrandDomain(MOL, STRAND, ATOM) (MOL->getStrand(STRAND)->getAtom(ATOM)->domain.get())
#define KgetStrandDomain(MOL, STRAND, ATOM) ((~MOL->getStrand(STRAND)->getAtom(ATOM)->domain).get())


TEST(assembly, leftOverHang) {
    std::vector<std::string> dictionary;
    
    auto mol = assembly::parseMolecule(std::string("<A>.[B]"), dictionary);

    EXPECT_EQ(mol->size(),       2);
    EXPECT_EQ(dictionary.size(), 2);
    
    // strands lenghts
    EXPECT_EQ(mol->getStrand(0)->length(), 1);
    EXPECT_EQ(mol->getStrand(1)->length(), 2);

    // check domains
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 0)), "B");
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 1, 1)), "B");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 1, 0)), "A");

    //check links
    EXPECT_EQ(getStrandAtom(mol, 0, 0)->partner, getStrandAtom(mol, 1, 1)); //?
    EXPECT_EQ(getStrandAtom(mol, 1, 1)->partner, getStrandAtom(mol, 0, 0)); //?

    EXPECT_TRUE(getStrandAtom(mol, 1, 0)->partner == NULL);
}

TEST(assembly, rightOverHang) {
    std::vector<std::string> dictionary;
    
    auto mol = assembly::parseMolecule(std::string("[B].<A>"), dictionary);

    EXPECT_EQ(mol->size(),       2);
    EXPECT_EQ(dictionary.size(), 2);
    
    // strands lenghts
    EXPECT_EQ(mol->getStrand(0)->length(), 1);
    EXPECT_EQ(mol->getStrand(1)->length(), 2);

    // check domains
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 0)), "B");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 1, 1)), "A");
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 1, 0)), "B");

    //check links
    EXPECT_EQ(getStrandAtom(mol, 0, 0)->partner, getStrandAtom(mol, 1, 0));
    EXPECT_EQ(getStrandAtom(mol, 1, 0)->partner, getStrandAtom(mol, 0, 0));

    EXPECT_TRUE(getStrandAtom(mol, 1, 1)->partner == NULL);
}

TEST(assembly, doubleStrand) {
    std::vector<std::string> dictionary;
    
    auto mol = assembly::parseMolecule(std::string("[AB]"), dictionary);

    EXPECT_EQ(mol->size(),       2);
    EXPECT_EQ(dictionary.size(), 2);
    
    // strands lenghts
    EXPECT_EQ(mol->getStrand(0)->length(), 2);
    EXPECT_EQ(mol->getStrand(1)->length(), 2);

    // check domains
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 0)), "A");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 1)), "B");
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 1, 0)), "A");
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 1, 1)), "B");

    //check links
    EXPECT_EQ(getStrandAtom(mol, 0, 0)->partner, getStrandAtom(mol, 1, 0));
    EXPECT_EQ(getStrandAtom(mol, 1, 0)->partner, getStrandAtom(mol, 0, 0));

    EXPECT_EQ(getStrandAtom(mol, 0, 1)->partner, getStrandAtom(mol, 1, 1));
    EXPECT_EQ(getStrandAtom(mol, 1, 1)->partner, getStrandAtom(mol, 0, 1));
}

TEST(assembly, bottomStrand) {
    std::vector<std::string> dictionary;
    
    auto mol = assembly::parseMolecule(std::string("{AB}"), dictionary);

    EXPECT_EQ(mol->size(),       1);
    EXPECT_EQ(dictionary.size(), 2);
    
    // strands lenghts
    EXPECT_EQ(mol->getStrand(0)->length(), 2);

    // check domains
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 0)), "A");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 1)), "B");

    //check links
    EXPECT_TRUE(getStrandAtom(mol, 0, 0)->partner == NULL);
    EXPECT_TRUE(getStrandAtom(mol, 0, 1)->partner == NULL);
}

TEST(assembly, bottomStrandKomplement) {
    std::vector<std::string> dictionary;
    
    auto mol = assembly::parseMolecule(std::string("{A*B}"), dictionary);

    EXPECT_EQ(mol->size(),       1);
    EXPECT_EQ(dictionary.size(), 2);
    
    // strands lenghts
    EXPECT_EQ(mol->getStrand(0)->length(), 2);

    // check domains
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 0, 0)), "A");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 1)), "B");

    //check links
    EXPECT_TRUE(getStrandAtom(mol, 0, 0)->partner == NULL);
    EXPECT_TRUE(getStrandAtom(mol, 0, 1)->partner == NULL);
}

TEST(assembly, doubleStrandKomplement) {
    std::vector<std::string> dictionary;
    
    auto mol = assembly::parseMolecule(std::string("[A*B]"), dictionary);

    EXPECT_EQ(mol->size(),       2);
    EXPECT_EQ(dictionary.size(), 2);
    
    // strands lenghts
    EXPECT_EQ(mol->getStrand(0)->length(), 2);
    EXPECT_EQ(mol->getStrand(1)->length(), 2);

    // check domains
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 0, 0)), "A");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 1)), "B");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 1, 0)), "A");
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 1, 1)), "B");

    //check links
    EXPECT_EQ(getStrandAtom(mol, 0, 0)->partner, getStrandAtom(mol, 1, 0));
    EXPECT_EQ(getStrandAtom(mol, 1, 0)->partner, getStrandAtom(mol, 0, 0));

    EXPECT_EQ(getStrandAtom(mol, 0, 1)->partner, getStrandAtom(mol, 1, 1));
    EXPECT_EQ(getStrandAtom(mol, 1, 1)->partner, getStrandAtom(mol, 0, 1));
}

TEST(assembly, OverhangKomplement) {
    std::vector<std::string> dictionary;
    
    auto mol = assembly::parseMolecule(std::string("[A*B].<A*>"), dictionary);

    EXPECT_EQ(mol->size(),       2);
    EXPECT_EQ(dictionary.size(), 2);
    
    // strands lenghts
    EXPECT_EQ(mol->getStrand(0)->length(), 2);
    EXPECT_EQ(mol->getStrand(1)->length(), 3);

    // check domains
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 0, 0)), "A");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 1)), "B");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 1, 0)), "A");
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 1, 1)), "B");

    //overhang
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 1, 2)), "A");

    //check links
    EXPECT_EQ(getStrandAtom(mol, 0, 0)->partner, getStrandAtom(mol, 1, 0));
    EXPECT_EQ(getStrandAtom(mol, 1, 0)->partner, getStrandAtom(mol, 0, 0));

    EXPECT_EQ(getStrandAtom(mol, 0, 1)->partner, getStrandAtom(mol, 1, 1));
    EXPECT_EQ(getStrandAtom(mol, 1, 1)->partner, getStrandAtom(mol, 0, 1));

    EXPECT_TRUE(getStrandAtom(mol, 1, 2)->partner == NULL);
}

TEST(assembly, LotOfKomplement) {
    std::vector<std::string> dictionary;
    
    auto mol = assembly::parseMolecule(std::string("{A**B*A***A****A******A*******}"), dictionary);

    EXPECT_EQ(mol->size(),       1);
    EXPECT_EQ(dictionary.size(), 2);
    
    // strands lenghts
    EXPECT_EQ(mol->getStrand(0)->length(), 6);

    // check domains
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 0)), "A");
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 0, 1)), "B");
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 0, 2)), "A");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 3)), "A");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 4)), "A");
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 0, 5)), "A");


    //check links
    EXPECT_TRUE(getStrandAtom(mol, 0, 0)->partner == NULL);
    EXPECT_TRUE(getStrandAtom(mol, 0, 1)->partner == NULL);
    EXPECT_TRUE(getStrandAtom(mol, 0, 2)->partner == NULL);
    EXPECT_TRUE(getStrandAtom(mol, 0, 3)->partner == NULL);
    EXPECT_TRUE(getStrandAtom(mol, 0, 4)->partner == NULL);
    EXPECT_TRUE(getStrandAtom(mol, 0, 5)->partner == NULL);
}

TEST(assembly, OverhangDual) {
    std::vector<std::string> dictionary;
    
    auto mol = assembly::parseMolecule(std::string("<B>.[A*B].<A*>"), dictionary);

    EXPECT_EQ(mol->size(),       2);
    EXPECT_EQ(dictionary.size(), 2);
    
    // strands lenghts
    EXPECT_EQ(mol->getStrand(0)->length(), 2);
    EXPECT_EQ(mol->getStrand(1)->length(), 4);

    // check domains
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 0, 0)), "A");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 1)), "B");

    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 1, 0)), "B");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 1, 1)), "A");
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 1, 2)), "B");
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 1, 3)), "A");

    //check links
    EXPECT_EQ(getStrandAtom(mol, 0, 0)->partner, getStrandAtom(mol, 1, 1));
    EXPECT_EQ(getStrandAtom(mol, 1, 1)->partner, getStrandAtom(mol, 0, 0));

    EXPECT_EQ(getStrandAtom(mol, 0, 1)->partner, getStrandAtom(mol, 1, 2));
    EXPECT_EQ(getStrandAtom(mol, 1, 2)->partner, getStrandAtom(mol, 0, 1));

    EXPECT_TRUE(getStrandAtom(mol, 1, 3)->partner == NULL);
    EXPECT_TRUE(getStrandAtom(mol, 1, 0)->partner == NULL);
}

TEST(assembly, complex) {
    std::vector<std::string> dictionary;
    
    auto mol = assembly::parseMolecule(std::string("<B>.[A*B].<A*>{AH}[G].<O>"), dictionary);

    EXPECT_EQ(mol->size(),       3);
    EXPECT_EQ(dictionary.size(), 5);
    
    // strands lenghts
    EXPECT_EQ(mol->getStrand(0)->length(), 5);
    EXPECT_EQ(mol->getStrand(1)->length(), 4);
    EXPECT_EQ(mol->getStrand(2)->length(), 2);

    // check domains
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 0, 0)), "A");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 1)), "B");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 2)), "A");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 3)), "H");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 0, 4)), "G");

    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 1, 0)), "B");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 1, 1)), "A");
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 1, 2)), "B");
    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 1, 3)), "A");

    EXPECT_EQ(dictionary.at(KgetStrandDomain(mol, 2, 0)), "G");
    EXPECT_EQ(dictionary.at( getStrandDomain(mol, 2, 1)), "O");

    //check links
    EXPECT_EQ(getStrandAtom(mol, 0, 0)->partner, getStrandAtom(mol, 1, 1));
    EXPECT_EQ(getStrandAtom(mol, 1, 1)->partner, getStrandAtom(mol, 0, 0));

    EXPECT_EQ(getStrandAtom(mol, 0, 1)->partner, getStrandAtom(mol, 1, 2));
    EXPECT_EQ(getStrandAtom(mol, 1, 2)->partner, getStrandAtom(mol, 0, 1));


    EXPECT_EQ(getStrandAtom(mol, 1, 3)->partner, getStrandAtom(mol, 0, 2));

    EXPECT_TRUE(getStrandAtom(mol, 1, 0)->partner == NULL);

    EXPECT_TRUE(getStrandAtom(mol, 0, 2)->partner == NULL);
    EXPECT_TRUE(getStrandAtom(mol, 0, 3)->partner == NULL);

    EXPECT_EQ(getStrandAtom(mol, 2, 0)->partner, getStrandAtom(mol, 0, 4));
    EXPECT_EQ(getStrandAtom(mol, 0, 4)->partner, getStrandAtom(mol, 2, 0));

    EXPECT_TRUE(getStrandAtom(mol, 2, 1)->partner == NULL);

}