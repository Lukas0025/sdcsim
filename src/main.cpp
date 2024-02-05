#include <iostream>
#include "domain.h"
#include "strand.h"

int main(int argc, char *argv[])
{
   auto chain1 = Strand();

   chain1.addDomain(Domain(0));
   chain1.addDomain(Domain(1));
   chain1.addDomain(Domain(2));
   chain1.addDomain(Domain(3));

   chain1.alignmentPrint();
   
   return 0;
}