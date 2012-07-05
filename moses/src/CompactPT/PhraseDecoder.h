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
  protected:
    
    friend class PhraseDictionaryMemoryHashed;
    
    typedef std::pair<unsigned char, unsigned char> AlignPoint;
    typedef std::pair<unsigned, unsigned> SrcTrg;
        
    enum Coding { None, REnc, PREnc } m_coding;
    
    size_t m_numScoreComponent;
    bool m_containsAlignmentInfo;
    size_t m_maxRank;
    size_t m_maxPhraseLength;
    
    boost::unordered_map<std::string, unsigned> m_sourceSymbolsMap;
    StringVector<unsigned char, unsigned, std::allocator> m_sourceSymbols;
    StringVector<unsigned char, unsigned, std::allocator> m_targetSymbols;
    
    std::vector<size_t> m_lexicalTableIndex;
    std::vector<SrcTrg> m_lexicalTable;
    
    CanonicalHuffman<unsigned>* m_symbolTree;
    CanonicalHuffman<AlignPoint>* m_alignTree;
    
    bool m_multipleScoreTrees;
    std::vector<CanonicalHuffman<float>*> m_scoreTrees;
    
    TargetPhraseCollectionCache m_decodingCache;
    
    PhraseDictionaryMemoryHashed& m_phraseDictionary;   
    
    // ***********************************************
    
    const std::vector<FactorType>* m_input;
    const std::vector<FactorType>* m_output;
    const PhraseDictionaryFeature* m_feature;
    const std::vector<float>* m_weight;
    float m_weightWP;
    const LMList* m_languageModels;
    
    std::string m_separator;
  
    // ***********************************************
    
    unsigned getSourceSymbolId(std::string& s);
    std::string getTargetSymbol(unsigned id) const;
    
    size_t getREncType(unsigned encodedSymbol);
    size_t getPREncType(unsigned encodedSymbol);
    
    unsigned getTranslation(unsigned srcIdx, size_t rank);
    
    size_t getMaxSourcePhraseLength();
    
    unsigned decodeREncSymbol1(unsigned encodedSymbol);
    unsigned decodeREncSymbol2Rank(unsigned encodedSymbol);
    unsigned decodeREncSymbol2Position(unsigned encodedSymbol);
    unsigned decodeREncSymbol3(unsigned encodedSymbol);
    
    unsigned decodePREncSymbol1(unsigned encodedSymbol);
    int decodePREncSymbol2Left(unsigned encodedSymbol);
    int decodePREncSymbol2Right(unsigned encodedSymbol);
    unsigned decodePREncSymbol2Rank(unsigned encodedSymbol);
    
    std::string makeSourceKey(std::string &);
    
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
    
    TargetPhraseVectorPtr createTargetPhraseCollection(const Phrase &sourcePhrase,
                                                       bool topLevel = false);
    
    TargetPhraseVectorPtr decodeCollection(TargetPhraseVectorPtr tpv,
                                           BitStream<> &encodedBitStream,
                                           const Phrase &sourcePhrase,
                                           bool topLevel);
    
    void pruneCache() {
      m_decodingCache.prune();
    }
    
};

}

#endif