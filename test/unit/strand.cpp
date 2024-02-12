#include <gtest/gtest.h>
#include "strand.h"
#include <stdexcept>

TEST(Strand, addDomains) {
  auto d1 = Domain(1);
  auto d2 = Domain(2);
  auto d3 = Domain(3);

  auto strand = Strand();

  strand.addDomain(d1);
  auto s2 = strand.addDomain(d2);
  auto s3 = strand.addDomain(d3);

  EXPECT_TRUE(strand.getAtom(s2)->domain == d2);
  EXPECT_TRUE(strand.getAtom(s3)->domain == d3);
  EXPECT_TRUE(strand.length()            == 3);
  EXPECT_TRUE(strand.getAtom(2)          == strand.getAtom(s3));
}

TEST(Strand, pairStrands) {
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

  //do pairing
  strand1.pairDomain(strand1.getAtom(0), strand2.getAtom(0));
  strand1.pairDomain(strand1.getAtom(1), strand2.getAtom(1));

  EXPECT_EQ(strand1.getAtom(0)->partner, strand2.getAtom(0));
  EXPECT_EQ(strand1.getAtom(1)->partner, strand2.getAtom(1));

  EXPECT_EQ(strand2.getAtom(0)->partner, strand1.getAtom(0));
  EXPECT_EQ(strand2.getAtom(1)->partner, strand1.getAtom(1));

  EXPECT_TRUE(strand1.getAtom(2)->partner == NULL);
  EXPECT_TRUE(strand2.getAtom(2)->partner == NULL);
}

TEST(Strand, saftyReadOnlyStrand) {
  auto d1 = Domain(1);
  auto d2 = Domain(2);
  auto d3 = Domain(3);

  auto strand1 = Strand();

  strand1.addDomain(d1);
  strand1.addDomain(d2);
  strand1.addDomain(d3);

  //I Want pointer this makes strand readOnly
  auto pointer = strand1.getAtom(0);

  //Now I cant add to strand because its readonly
  try {
    strand1.addDomain(~d1);
    FAIL() << "addDomain() should throw an error\n";
  } catch (std::runtime_error& exception) {
    //ok here
  } catch (...) {
    FAIL() << "ERROR: Unexpected exception thrown: " << std::current_exception << std::endl;
  }

  // pointer cant be changed
  EXPECT_EQ(pointer, strand1.getAtom(0));
}

TEST(Strand, halfPairStrands) {
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

  //do pairing
  strand1.halfPairDomain(strand1.getAtom(0), strand2.getAtom(0));
  strand1.halfPairDomain(strand1.getAtom(1), strand2.getAtom(1));

  EXPECT_EQ(strand1.getAtom(0)->partner, strand2.getAtom(0));
  EXPECT_EQ(strand1.getAtom(1)->partner, strand2.getAtom(1));

  EXPECT_TRUE(strand2.getAtom(0)->partner == NULL);
  EXPECT_TRUE(strand2.getAtom(1)->partner == NULL);

  EXPECT_TRUE(strand1.getAtom(2)->partner == NULL);
  EXPECT_TRUE(strand2.getAtom(2)->partner == NULL);
}