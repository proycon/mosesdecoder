#include <iostream>
#include "TypeDef.h"
#include "CompactPT/PhrasetableCreator.h"
#include "CompactPT/CanonicalHuffman.h"

using namespace Moses;

void printHelp()
{
  std::cerr << "Usage:\n"
            "options: \n"
            "\t-in  string -- input table file name\n"
            "\t-out string -- prefix of binary table file\n"
            "\t-enc string -- Encoding type (None REnc PREnc)\n"
            "\n";
}

int main(int argc,char **argv) {
  std::cerr << "processPhraseTableHashed by Marcin Junczys-Dowmunt\n";
    
  std::string inFilePath;
  std::string outFilePath("out");
  PhrasetableCreator::Coding coding = PhrasetableCreator::None;
  
  if(1 >= argc) {
    printHelp();
    return 1;
  }
  for(int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if("-in" == arg && i+1 < argc) {
      ++i;
      inFilePath = argv[i];
    } else if("-out" == arg && i+1 < argc) {
      ++i;
      outFilePath = argv[i];
    } else if("-enc" == arg && i+1 < argc) {
      ++i;
      std::string val(argv[i]);
      if(val == "None") {
        coding = PhrasetableCreator::None;
      }
      else if(val == "REnc") {
        coding = PhrasetableCreator::REnc;
      }
      else if(val == "PREnc") {
        coding = PhrasetableCreator::PREnc;
      }
    } else {
      //somethings wrong... print help
      printHelp();
      return 1;
    }
  }
  
  size_t numScoreComponent = 5;  
  PhrasetableCreator(inFilePath, outFilePath, coding, 10, 16);
}
