#include <iostream>
#include "TypeDef.h"
#include "PhraseDictionaryMemoryHashed.h"
#include "CompactPT/PhrasetableCreator.h"

#include <boost/date_time/posix_time/posix_time.hpp>

using namespace Moses;

void printHelp()
{
  std::cerr << "Usage:\n"
            "options: \n"
            "\t-in  string -- input table file name\n"
            "\n";
}

int main(int argc,char **argv) {
  std::cerr << "phraseTableTester by Marcin Junczys-Dowmunt\n";
   
  if(1 >= argc) {
    printHelp();
    return 1;
  }
  for(int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if("-in" == arg && i+1 < argc) {
      ++i;
      inFilePath = argv[i];
    }
    else {
      //somethings wrong... print help
      printHelp();
      return 1;
    }
  }
  
  size_t numScoreComponent = 5;  
  PhraseDictionaryMemoryHashed pt(numScoreComponent, MemoryHashedBinary, NULL);
  
  std::vector<FactorType> input;
  input.push_back(0);
  
  std::vector<FactorType> output;
  output.push_back(0);
  
  std::vector<float> weight(numScoreComponent, 0.2);
  size_t tableLimit = 20;
  float weightWP = 1;
  
  LMList languageModels;
  
  pt.Load(input, output, inFilePath, weight, tableLimit, languageModels, weightWP);
  
  //std::FILE* pFile = std::fopen(inFilePath.c_str() , "r");
  
  //BlockHashIndex hash(10, 16);
  //hash.LoadIndex(pFile);
  
  std::string line;
  while(std::getline(std::cin, line)) {
    Phrase sentence(Input,0);
    
    sentence.CreateFromString(input, line, "|||");
    std::cout << sentence << std::endl;
    
    boost::posix_time::ptime pt_start = boost::posix_time::microsec_clock::local_time();
    
    for(size_t i = 0; i < sentence.GetSize(); i++)
      for(size_t j = i; j < sentence.GetSize() && j-i < 7; j++) {
        const TargetPhraseCollection* tpc
          = pt.GetTargetPhraseCollection(sentence.GetSubString(WordsRange(i, j)));
        //if(tpc != NULL)
        //  std::cout << tpc->GetSize() << std::endl;
        //std::string sourcePhraseString = sentence.GetSubString(WordsRange(i, j)).GetStringRep(input);
        //size_t index = hash[sourcePhraseString];
        
      }
      
    boost::posix_time::ptime pt_end = boost::posix_time::microsec_clock::local_time();
    boost::posix_time::time_duration delta = pt_end - pt_start;
    int d = delta.total_milliseconds();
  
    std::cout << "Time: " << d << std::endl;
    
    pt.CleanUp();
    //hash.KeepNLastRanges(0.01, 0.2);
   }
}
