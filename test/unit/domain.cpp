#include <gtest/gtest.h>
#include "domain.h"

TEST(Domain, Complementarity) {
  auto d1 = Domain(1);
  auto c1 = ~Domain(1);

  EXPECT_TRUE(~d1 ==  c1);
  EXPECT_TRUE( d1 == ~c1);
  EXPECT_TRUE( d1 == ~~d1);
  EXPECT_TRUE( d1 != c1);
}

TEST(Domain, Overwrap) {
  Domain domains[100];

  //okrajové podminky
  domains[0].set(32767);
  domains[1].set(0);

  // test uvnitř
  for (unsigned i = 2; i < 100; i++) {
    domains[i].set(i * 327);
  }

  for (unsigned i = 0; i < 100; i++) {
    for (unsigned j = 0; j < 100; j++) {
      if (i == j) continue;

      EXPECT_TRUE( domains[i] != domains[j]);
      EXPECT_TRUE(~domains[i] != domains[j]);
      EXPECT_TRUE( domains[i] != ~domains[j]);
      EXPECT_TRUE(~domains[i] != ~domains[j]);
    }
  }
}

TEST(Domain, SetGet) {

  Domain d = Domain();

  d.set(32767);

  EXPECT_TRUE(d.get() == 32767);

  d.set(0);

  EXPECT_TRUE(d.get() == 0);
}