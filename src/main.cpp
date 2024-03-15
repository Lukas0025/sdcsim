#include <iostream>
#include <argparse/argparse.hpp>
#include "molecule.h"
#include "assembly.h"
#include "output.h"
#include "register.h"
#include "svg.h"

inline void printRegisters(std::vector<Register*> &registers, std::vector<std::string> &dictionary, const std::string &outputType, std::vector<std::pair<std::string, std::string>> &macros, Svg &svg, const char* space) {
   if (outputType == "ascii") {
      for (auto &reg : registers) {
         output::asciiPrint(*reg->get(), space, dictionary);
         std::cout << '\n';
      }
   } else if (outputType == "raw") {
      for (auto &reg : registers) {
         output::rawPrint(*reg->get(), space, dictionary);
         std::cout << '\n';
      }
   } else if (outputType == "dummy") {
      for (auto &reg : registers) {
         output::dummyPrint(*reg->get(), space, dictionary);
         std::cout << '\n';
      }
   } else if (outputType == "assembly") {
      for (auto &reg : registers) {
         output::assemblyPrint(*reg->get(), space, dictionary, macros);
         std::cout << '\n';
      }
   } else if (outputType == "svg") {
      for (auto &reg : registers) {
         output::svgPrint(*reg->get(), dictionary, svg);
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
       .default_value(std::string{"ascii"})
       .choices("ascii", "raw", "svg", "dummy", "assembly");

   args.add_argument("--silent")
       .help("Print only necessary outputs")
       .default_value(false)
       .implicit_value(true);

   args.add_argument("-d", "--decode")
       .help("Print decoded version of registers (by macros)")
       .default_value(false)
       .implicit_value(true);

   args.add_argument("-a", "--all")
       .help("Show all instruction interactions")
       .default_value(false)
       .implicit_value(true);

   args.add_argument("-c", "--color")
       .help("Select color scheme for svg output")
       .default_value(std::string{"domain"})
       .choices("domain", "chain", "black");

   args.add_argument("-b", "--break")
       .help("Break before instruction")
       .scan<'i', int>()
       .default_value(INT_MAX);

   args.add_argument("-t", "--time")
       .help("registr expose time to instruction")
       .scan<'i', int>()
       .default_value(INT_MAX);
   
   args.add_argument("-s", "--spaceing")
       .help("Space char/s between domains in ascii formats")
       .default_value("");

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
   auto spaceing   = (args.get<std::string>("-s")).c_str();
   auto breakAt    = args.get<int>("-b");
   auto time       = args.get<int>("-t");
   bool decode     = args["-d"] == true;
   bool all        = args["-a"] == true;
   bool silent     = args["--silent"] == true;
   uint8_t colorM  = CM_DOMAIN;

   if (args.get<std::string>("-c") == "chain") colorM  = CM_CHAIN;
   if (args.get<std::string>("-c") == "black") colorM  = CM_BLACK;


   std::vector<std::vector<Molecule*>>              instructions;
   std::vector<Register*>                           registers;
   std::vector<std::string>                         dictionary;
   std::vector<std::pair<std::string, std::string>> macros;

   //parse assembly file
   assembly::parse(sdcAsmFile.c_str(), registers, instructions, dictionary, macros);

   auto svg = Svg(18, 4, 12, colorM, dictionary);

   if (!silent && outputType != "svg") {
      std::cout << "---------------------------------\n";
      std::cout << "|        INITAL REGISTERS       |\n";
      std::cout << "---------------------------------\n";
      std::cout << '\n';
   }

   printRegisters(registers, dictionary, outputType, macros, svg, spaceing);
   svg.inscructionDone("inital");

   if (decode) {
      std::cout << "----------- DECODED -------------\n";

      printRegisters(registers, dictionary, "assembly", macros, svg, spaceing);
   }

   int idInstruction = 0;

   // apply instruction
   for (auto &instruction : instructions) {
      
      if (breakAt <= idInstruction) break;

      for (auto &reg : registers) {
         reg->applyInstruction(instruction);
      }

      if (all) {
         if (!silent && outputType != "svg") {
            std::cout << '\n';
            std::cout << "---------------------------------\n";
            std::cout << "|        INSTRUCTION " << std::setfill('0') << std::setw(3) << idInstruction << "        |\n";
            std::cout << "---------------------------------\n";
            std::cout << '\n';
         }

         printRegisters(registers, dictionary, outputType, macros, svg, spaceing);  
         svg.inscructionDone(("instruction " + std::to_string(idInstruction)).c_str());
      }
      
      idInstruction++;
   }

   if (!silent) {
      if (outputType != "svg") {
         std::cout << '\n';
         std::cout << "---------------------------------\n";
         std::cout << "|             FINAL             |\n";
         std::cout << "---------------------------------\n";
         std::cout << '\n';
      }

      printRegisters(registers, dictionary, outputType, macros, svg, spaceing);
      svg.inscructionDone("final");
   }

   if (decode) {
      std::cout << "----------- DECODED -------------\n";

      printRegisters(registers, dictionary, "assembly", macros, svg, spaceing);
   }

   //print output
   if (outputType == "svg") {
      std::cout << svg.get();
   }
   
   
   return 0;
}