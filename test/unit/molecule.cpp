#include <gtest/gtest.h>
#include "molecule.h"
#include <stdexcept>

TEST(Molecule, addStrand) {
  auto d1 = Domain(1);
  auto d2 = Domain(2);
  auto d3 = Domain(3);

  auto strand1 = Strand();
  auto strand2 = Strand();

  strand1.addDomain(d1);
  strand1.addDomain(d2);
  strand1.addDomain(d3);
 
  strand2.addDomain(~d1);
  strand2.addDomain(~d2);
  strand2.addDomain(d3);

  //simply create molecule
  auto mol = Molecule(&strand2);

  mol.addStrand(&strand1, 0);

  EXPECT_EQ(mol.getStrand(1)->getAtom(0)->partner, mol.getStrand(0)->getAtom(0));
  EXPECT_EQ(mol.getStrand(1)->getAtom(1)->partner, mol.getStrand(0)->getAtom(1));

  EXPECT_TRUE(mol.getStrand(0)->getAtom(0)->partner == NULL);
  EXPECT_TRUE(mol.getStrand(0)->getAtom(1)->partner == NULL);

  EXPECT_TRUE(mol.getStrand(1)->getAtom(2)->partner == NULL);
  EXPECT_TRUE(mol.getStrand(0)->getAtom(2)->partner == NULL);

}

TEST(Molecule, addStrandNegative) {
  auto d1 = Domain(1);
  auto d2 = Domain(2);
  auto d3 = Domain(3);

  auto strand1 = Strand();
  auto strand2 = Strand();

  strand1.addDomain(d1);
  strand1.addDomain(d2);
  strand1.addDomain(d3);
 
  strand2.addDomain(~d1);
  strand2.addDomain(~d1);
  strand2.addDomain(d3);

  //simply create molecule
  auto mol = Molecule(&strand1);

  mol.addStrand(&strand2, -1);

  EXPECT_TRUE(mol.getStrand(1)->getAtom(0)->partner == NULL);
  EXPECT_EQ(mol.getStrand(1)->getAtom(1)->partner, mol.getStrand(0)->getAtom(0));

  EXPECT_TRUE(mol.getStrand(0)->getAtom(0)->partner == NULL);
  EXPECT_TRUE(mol.getStrand(0)->getAtom(1)->partner == NULL);

  EXPECT_TRUE(mol.getStrand(1)->getAtom(2)->partner == NULL);
  EXPECT_TRUE(mol.getStrand(0)->getAtom(2)->partner == NULL);

}

TEST(Molecule, donePairStrands) {
  auto d1 = Domain(1);
  auto d2 = Domain(2);
  auto d3 = Domain(3);

  auto strand1 = Strand();
  auto strand2 = Strand();

  strand1.addDomain(d1);
  strand1.addDomain(d2);
  strand1.addDomain(d3);
 
  strand2.addDomain(~d1);
  strand2.addDomain(~d2);
  strand2.addDomain(d3);

  //simply create molecule
  auto mol = Molecule(&strand2);

  mol.addStrand(&strand1, 0);

  mol.donePairStrands(1, 0, strand1.length() - 1);

  EXPECT_EQ(mol.getStrand(1)->getAtom(0)->partner, mol.getStrand(0)->getAtom(0));
  EXPECT_EQ(mol.getStrand(1)->getAtom(1)->partner, mol.getStrand(0)->getAtom(1));

  EXPECT_EQ(mol.getStrand(0)->getAtom(0)->partner, mol.getStrand(1)->getAtom(0));
  EXPECT_EQ(mol.getStrand(0)->getAtom(1)->partner, mol.getStrand(1)->getAtom(1));

  EXPECT_TRUE(mol.getStrand(1)->getAtom(2)->partner == NULL);
  EXPECT_TRUE(mol.getStrand(0)->getAtom(2)->partner == NULL);

}

TEST(Molecule, uncompletePairStrands) {
  auto d1 = Domain(1);
  auto d2 = Domain(2);
  auto d3 = Domain(3);

  auto strand1 = Strand();
  auto strand2 = Strand();

  strand1.addDomain(d1);
  strand1.addDomain(d2);
  strand1.addDomain(d3);
 
  strand2.addDomain(~d1);
  strand2.addDomain(~d2);
  strand2.addDomain(d3);

  //simply create molecule
  auto mol = Molecule(&strand2);

  mol.addStrand(&strand1, 0);

  mol.pairUncomplete(1);

  EXPECT_EQ(mol.getStrand(1)->getAtom(0)->partner, mol.getStrand(0)->getAtom(0));
  EXPECT_EQ(mol.getStrand(1)->getAtom(1)->partner, mol.getStrand(0)->getAtom(1));

  EXPECT_EQ(mol.getStrand(0)->getAtom(0)->partner, mol.getStrand(1)->getAtom(0));
  EXPECT_EQ(mol.getStrand(0)->getAtom(1)->partner, mol.getStrand(1)->getAtom(1));

  EXPECT_TRUE(mol.getStrand(1)->getAtom(2)->partner == NULL);
  EXPECT_TRUE(mol.getStrand(0)->getAtom(2)->partner == NULL);

}

TEST(Molecule, uncompletePairTwoStrands) {
  auto d1 = Domain(1);
  auto d2 = Domain(2);
  auto d3 = Domain(3);

  auto strand1 = Strand();
  auto strand2 = Strand();
  auto strand3 = Strand();

  strand1.addDomain(d1);
  strand1.addDomain(d2);
  strand1.addDomain(d3);
 
  strand2.addDomain(~d1);
  strand2.addDomain(~d2);
  strand2.addDomain(d3);

  strand3.addDomain(d1);
  strand3.addDomain(~d2);
  strand3.addDomain(~d3);

  //simply create molecule
  auto mol = Molecule(&strand1);

  mol.addStrand(&strand2, 0);
  mol.addStrand(&strand3, 0);

  mol.pairUncomplete(1);
  mol.pairUncomplete(2);

  EXPECT_EQ(mol.getStrand(1)->getAtom(0)->partner, mol.getStrand(0)->getAtom(0));
  EXPECT_EQ(mol.getStrand(1)->getAtom(1)->partner, mol.getStrand(0)->getAtom(1));

  EXPECT_EQ(mol.getStrand(0)->getAtom(0)->partner, mol.getStrand(1)->getAtom(0));
  EXPECT_EQ(mol.getStrand(0)->getAtom(1)->partner, mol.getStrand(1)->getAtom(1));

  EXPECT_TRUE(mol.getStrand(1)->getAtom(2)->partner == NULL);

  EXPECT_EQ(mol.getStrand(0)->getAtom(2)->partner, mol.getStrand(2)->getAtom(2));

}

TEST(Molecule, addStrandNegative2) {
  auto d1 = Domain(1);
  auto d2 = Domain(2);

  auto strand1 = Strand();
  auto strand2 = Strand();

  strand1.addDomain(d1);
 
  strand2.addDomain(d2);
  strand2.addDomain(~d1);

  //simply create molecule
  auto mol = Molecule(&strand1);

  mol.addStrand(&strand2, -1);

  EXPECT_TRUE(mol.getStrand(1)->getAtom(0)->partner == NULL);
  EXPECT_EQ(mol.getStrand(1)->getAtom(1)->partner, mol.getStrand(0)->getAtom(0));

}