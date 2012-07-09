#ifndef moses_LexicalReorderingTableCompact_h
#define moses_LexicalReorderingTableCompact_h

#include "LexicalReorderingTable.h"
#include "StaticData.h"
#include "PhraseDictionary.h"
#include "GenerationDictionary.h"
#include "TargetPhrase.h"
#include "TargetPhraseCollection.h"

#include "CompactPT/BlockHashIndex.h"
#include "CompactPT/CanonicalHuffman.h"
#include "CompactPT/StringVector.h"

namespace Moses {

class LexicalReorderingTableCompact: public LexicalReorderingTable
{
  private:
    bool m_inMemory;
    
    BlockHashIndex m_hash;
    CanonicalHuffman<float>* m_scoreTree;
    StringVector<unsigned char, unsigned long, MmapAllocator>  m_scoresMapped;
    StringVector<unsigned char, unsigned long, std::allocator> m_scoresMemory;

    std::string MakeKey(const Phrase& f, const Phrase& e, const Phrase& c) const;
    std::string MakeKey(const std::string& f, const std::string& e, const std::string& c) const;
    
  public:
    LexicalReorderingTableCompact( const std::string& filePath,
                                const std::vector<FactorType>& f_factors,
                                const std::vector<FactorType>& e_factors,
                                const std::vector<FactorType>& c_factors);
  
    LexicalReorderingTableCompact(
                                const std::vector<FactorType>& f_factors,
                                const std::vector<FactorType>& e_factors,
                                const std::vector<FactorType>& c_factors);
  
    virtual ~LexicalReorderingTableCompact();

    virtual std::vector<float> GetScore(const Phrase& f, const Phrase& e, const Phrase& c);
    void Load(const std::string& filePath);
};

}

#endif