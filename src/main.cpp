#include <iostream>
#include <argparse/argparse.hpp>
#include "molecule.h"
#include "assembly.h"
#include "output.h"

int main(int argc, char *argv[])
{
   // init arg parser
   argparse::ArgumentParser args("sdcsim");

   // add args
   args.add_argument("sdcasm")
       .help("File with SDC Assembly code of SDC program.")
       .required();

   args.add_argument("-o", "--output")
       .help("Define output tipe of simulation.")
       .default_value(std::string{"ascii"})
       .choices("ascii", "raw", "svg");

   // parse args
   try {
      args.parse_args(argc, argv);
   } catch (const std::exception& err) {
      std::cerr << err.what() << std::endl;
      std::cerr << args;
      std::exit(1);
   }

   auto sdcAsmFile = args.get<std::string>("sdcasm");
   auto outputType = args.get<std::string>("-o");


   std::vector<Molecule*> test;
   std::vector<std::string> dictionary;

   assembly::parse(sdcAsmFile.c_str(), test, dictionary);

   if (outputType == "ascii") {
      output::asciiPrint(*test[0], "", dictionary);
   } else if (outputType == "raw") {
      std::cout << test[0]->size() << '\n';
      output::rawPrint(*test[0], "", dictionary);
   }
   
   
   return 0;
}