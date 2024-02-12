#include <iostream>
#include "molecule.h"
#include "assembly.h"
#include "output.h"

int main(int argc, char *argv[])
{
   std::vector<Molecule*> test;
   std::vector<std::string> dictionary;


   assembly::parse("sdcasm/rule110.sdcasm", test, dictionary);

   output::asciiPrint(*test[0], "", dictionary);
   
   
   return 0;
}