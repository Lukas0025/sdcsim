#include <iostream>
#include <argparse/argparse.hpp>
#include "molecule.h"
#include "assembly.h"
#include "output.h"
#include "register.h"

inline void printRegisters(std::vector<Register*> &registers, std::vector<std::string> &dictionary, const std::string &outputType, std::vector<std::pair<std::string, std::string>> &macros) {
   if (outputType == "ascii") {
      for (auto &reg : registers) {
         output::asciiPrint(*reg->get(), "", dictionary);
         std::cout << '\n';
      }
   } else if (outputType == "raw") {
      for (auto &reg : registers) {
         output::rawPrint(*reg->get(), "", dictionary);
         std::cout << '\n';
      }
   } else if (outputType == "dummy") {
      for (auto &reg : registers) {
         output::dummyPrint(*reg->get(), "", dictionary);
         std::cout << '\n';
      }
   } else if (outputType == "assembly") {
      for (auto &reg : registers) {
         output::assemblyPrint(*reg->get(), "", dictionary, macros);
         std::cout << '\n';
      }
   } else if (outputType == "svg") {
      for (auto &reg : registers) {
         output::svgPrint(*reg->get(), dictionary);
         std::cout << '\n';
      }
   }
}

int main(int argc, char *argv[])
{
   // init arg parser
   argparse::ArgumentParser args("sdcsim");

   // add args
   args.add_argument("sdcasm")
       .help("File with SDC Assembly code of SDC program.")
       .required();

   args.add_argument("-f", "--format")
       .help("Define output fromat of simulation.")
       .default_value(std::string{"assembly"})
       .choices("ascii", "raw", "svg", "dummy", "assembly");

   args.add_argument("--silent")
       .help("Print only necessary outputs")
       .default_value(false)
       .implicit_value(true);

   // parse args
   try {
      args.parse_args(argc, argv);
   } catch (const std::exception& err) {
      std::cerr << err.what() << std::endl;
      std::cerr << args;
      std::exit(1);
   }

   auto sdcAsmFile = args.get<std::string>("sdcasm");
   auto outputType = args.get<std::string>("-f");
   bool silent     = args["--silent"] == true;


   std::vector<std::vector<Molecule*>>              instructions;
   std::vector<Register*>                           registers;
   std::vector<std::string>                         dictionary;
   std::vector<std::pair<std::string, std::string>> macros;

   //parse assembly file
   assembly::parse(sdcAsmFile.c_str(), registers, instructions, dictionary, macros);

   if (!silent) {
      std::cout << "--------------------------------\n";
      std::cout << "|       INITAL REGISTERS       |\n";
      std::cout << "--------------------------------\n";
      std::cout << '\n';
   }

   printRegisters(registers, dictionary, outputType, macros);

   // apply instruction
   for (auto &instruction : instructions) {
      for (auto &reg : registers) {
         reg->applyInstruction(instruction);
      }

      if (!silent) {
         std::cout << '\n';
         std::cout << "--------------------------------\n";
         std::cout << "|         INSTRUCTION          |\n";
         std::cout << "--------------------------------\n";
         std::cout << '\n';
      }

      printRegisters(registers, dictionary, outputType, macros);  
   }

   //print output
   //end
   
   
   return 0;
}