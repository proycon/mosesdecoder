#include <iostream>
#include "TypeDef.h"
#include "CompactPT/PhrasetableCreator.h"

using namespace Moses;

void printHelp()
{
  std::cerr << "Usage:\n"
            "options: \n"
            "\t-in  string -- input table file name\n"
            "\t-out string -- prefix of binary table file\n"
            "\n";
}

int main(int argc,char **argv) {
  std::cerr << "processPhraseTableHashed by Marcin Junczys-Dowmunt\n";
  
  std::string inFilePath;
  std::string outFilePath("out");
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
    } else {
      //somethings wrong... print help
      printHelp();
      return 1;
    }
  }

  PhrasetableCreator(inFilePath, outFilePath, 10, 16);
  
  //size_t numScoreComponent = 5;  
  //PhraseDictionaryMemoryHashed pt(numScoreComponent, MemoryHashedText, NULL);
  //pt.LoadText(inFilePath);
  //std::cerr << "Saving to " << outFilePath << std::endl;
  //pt.SaveBinary(outFilePath);
}
