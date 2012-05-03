#ifndef PHRASEDECODER_H__
#define PHRASEDECODER_H__

#include <sstream>
#include <vector>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <fstream>
#include <string>
#include <iterator>
#include <algorithm>
#include <sys/stat.h>

#include "TypeDef.h"
#include "FactorCollection.h"
#include "Word.h"
#include "Util.h"
#include "InputFileStream.h"
#include "StaticData.h"
#include "WordsRange.h"
#include "UserMessage.h"
#include "PhraseDictionaryMemoryHashed.h"

#include "StringVector.h"
#include "CanonicalHuffman.h"
#include "ConsistantPhrases.h"
#include "TargetPhraseCollectionCache.h"

namespace Moses {

class PhraseDictionaryMemoryHashed;

class PhraseDecoder {
  private:
    
    typedef std::pair<unsigned char, unsigned char> AlignPoint;
    typedef std::pair<unsigned, unsigned> SrcTrg;
        
    enum Coding { None, REnc, PREnc } m_coding;
    
    size_t m_numScoreComponent;
    bool m_containsAlignmentInfo;
    
    boost::unordered_map<std::string, unsigned> m_sourceSymbolsMap;
    StringVector<unsigned char, unsigned, std::allocator> m_sourceSymbols;
    StringVector<unsigned char, unsigned, std::allocator> m_targetSymbols;
    
    std::vector<size_t> m_lexicalTableIndex;
    std::vector<SrcTrg> m_lexicalTable;
    
    CanonicalHuffman<unsigned>* m_symbolTree;
    CanonicalHuffman<float>* m_scoreTree;
    CanonicalHuffman<AlignPoint>* m_alignTree;
    
    PhraseDictionaryMemoryHashed& m_phraseDictionary;   
    
    // ***********************************************
    
    const std::vector<FactorType>* m_input;
    const std::vector<FactorType>* m_output;
    const PhraseDictionaryFeature* m_feature;
    const std::vector<float>* m_weight;
    float m_weightWP;
    const LMList* m_languageModels;
  
    // ***********************************************
  
  public:
    
    unsigned getSourceSymbolId(std::string& s);
    std::string getTargetSymbol(unsigned id) const;
    
    size_t getREncType(unsigned encodedSymbol);
    size_t getPREncType(unsigned encodedSymbol);
    
    unsigned getTranslation(unsigned srcIdx, size_t rank);
    
    unsigned decodeREncSymbol1(unsigned encodedSymbol);
    unsigned decodeREncSymbol2Rank(unsigned encodedSymbol);
    unsigned decodeREncSymbol2Position(unsigned encodedSymbol);
    unsigned decodeREncSymbol3(unsigned encodedSymbol);
    
    unsigned decodePREncSymbol1(unsigned encodedSymbol);
    int decodePREncSymbol2Left(unsigned encodedSymbol);
    int decodePREncSymbol2Right(unsigned encodedSymbol);
    unsigned decodePREncSymbol2Rank(unsigned encodedSymbol);
    
  public:
    
    PhraseDecoder(
      PhraseDictionaryMemoryHashed &phraseDictionary,
      const std::vector<FactorType>* &input,
      const std::vector<FactorType>* &output,
      const PhraseDictionaryFeature* feature,
      size_t numScoreComponent,
      const std::vector<float>* weight,
      float weightWP,
      const LMList* languageModels
    );
    
    ~PhraseDecoder();
     
    size_t load(std::FILE* in);
    size_t save(std::FILE* out);
    
    TargetPhraseVectorPtr decodeCollection(std::string encoded,
                                             const Phrase &sourcePhrase,
                                             TargetPhraseCollectionCache &cache);
    
    void CleanUp(bool = false);
    
};

}

#endif