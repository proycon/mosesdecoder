#include <deque>

#include "PhraseDecoder.h"

namespace Moses {

PhraseDecoder::PhraseDecoder(
  PhraseDictionaryMemoryHashed &phraseDictionary,
  const std::vector<FactorType>* &input,
  const std::vector<FactorType>* &output,
  const PhraseDictionaryFeature* feature,
  size_t numScoreComponent,
  const std::vector<float>* weight,
  float weightWP,
  const LMList* languageModels
)
  : m_coding(None), m_numScoreComponent(numScoreComponent),
  m_containsAlignmentInfo(true), m_symbolTree(0), m_multipleScoreTrees(false),
  m_scoreTrees(1), m_alignTree(0),
  m_phraseDictionary(phraseDictionary), m_input(input), m_output(output),
  m_feature(feature), m_weight(weight),
  m_weightWP(weightWP), m_languageModels(languageModels),
  m_separator(" ||| ")
{ }

PhraseDecoder::~PhraseDecoder() {
  if(m_symbolTree)
    delete m_symbolTree;
  
  for(size_t i = 0; i < m_scoreTrees.size(); i++)
    if(m_scoreTrees[i])
      delete m_scoreTrees[i];
  
  if(m_alignTree)
    delete m_alignTree;
}

std::string& PhraseDecoder::getSeparator() {
  return m_separator;
}

inline unsigned PhraseDecoder::getSourceSymbolId(std::string& symbol) {
  boost::unordered_map<std::string, unsigned>::iterator it
    = m_sourceSymbolsMap.find(symbol);
  if(it != m_sourceSymbolsMap.end())
    return it->second;
    
  size_t idx = m_sourceSymbols.find(symbol);
  m_sourceSymbolsMap[symbol] = idx;
  return idx;
}

inline std::string PhraseDecoder::getTargetSymbol(unsigned idx) const {
  if(idx < m_targetSymbols.size())
    return m_targetSymbols[idx];
  return std::string("##ERROR##");
}

inline size_t PhraseDecoder::getREncType(unsigned encodedSymbol) {
  return (encodedSymbol >> 30) + 1;
}

inline size_t PhraseDecoder::getPREncType(unsigned encodedSymbol) {
  return (encodedSymbol >> 31) + 1;
}

inline unsigned PhraseDecoder::getTranslation(unsigned srcIdx, size_t rank) {
  size_t srcTrgIdx = m_lexicalTableIndex[srcIdx];
  return m_lexicalTable[srcTrgIdx + rank].second;
}

size_t PhraseDecoder::getMaxSourcePhraseLength() {
  return m_maxPhraseLength;
}

inline unsigned PhraseDecoder::decodeREncSymbol1(unsigned encodedSymbol) {
  return encodedSymbol &= ~(3 << 30);
}

inline unsigned PhraseDecoder::decodeREncSymbol2Rank(unsigned encodedSymbol) {
  return encodedSymbol &= ~(255 << 24);
}

inline unsigned PhraseDecoder::decodeREncSymbol2Position(unsigned encodedSymbol) {
  encodedSymbol &= ~(3 << 30);
  encodedSymbol >>= 24;
  return encodedSymbol;
}

inline unsigned PhraseDecoder::decodeREncSymbol3(unsigned encodedSymbol) {
  return encodedSymbol &= ~(3 << 30);
}

inline unsigned PhraseDecoder::decodePREncSymbol1(unsigned encodedSymbol) {
  return encodedSymbol &= ~(1 << 31);
}

inline int PhraseDecoder::decodePREncSymbol2Left(unsigned encodedSymbol) {
  return ((encodedSymbol >> 25) & 63) - 32;
}

inline int PhraseDecoder::decodePREncSymbol2Right(unsigned encodedSymbol) {
  return ((encodedSymbol >> 19) & 63) - 32;
}

inline unsigned PhraseDecoder::decodePREncSymbol2Rank(unsigned encodedSymbol) {
  return (encodedSymbol & 524287);
}

size_t PhraseDecoder::load(std::FILE* in) {
  size_t start = std::ftell(in);
  
  std::fread(&m_coding, sizeof(m_coding), 1, in);
  std::fread(&m_numScoreComponent, sizeof(m_numScoreComponent), 1, in);
  std::fread(&m_containsAlignmentInfo, sizeof(m_containsAlignmentInfo), 1, in);
  std::fread(&m_maxPhraseLength, sizeof(m_maxPhraseLength), 1, in);
  
  if(m_coding == REnc) {
    m_sourceSymbols.load(in);
    
    size_t size;
    std::fread(&size, sizeof(size_t), 1, in);
    m_lexicalTableIndex.resize(size);
    std::fread(&m_lexicalTableIndex[0], sizeof(size_t), size, in);
    
    std::fread(&size, sizeof(size_t), 1, in);
    m_lexicalTable.resize(size);
    std::fread(&m_lexicalTable[0], sizeof(SrcTrg), size, in);
  }
  
  m_targetSymbols.load(in);
  
  m_symbolTree = new CanonicalHuffman<unsigned>(in);
  
  std::fread(&m_multipleScoreTrees, sizeof(m_multipleScoreTrees), 1, in);
  if(m_multipleScoreTrees) {
    m_scoreTrees.resize(m_numScoreComponent);
    for(size_t i = 0; i < m_numScoreComponent; i++)
      m_scoreTrees[i] = new CanonicalHuffman<float>(in);
  }
  else {
    m_scoreTrees.resize(1);
    m_scoreTrees[0] = new CanonicalHuffman<float>(in);
  }
  
  if(m_containsAlignmentInfo)
    m_alignTree = new CanonicalHuffman<AlignPoint>(in);
  
  size_t end = std::ftell(in);
  return end - start;
}
    
TargetPhraseVectorPtr PhraseDecoder::decodeCollection(
  std::string encoded, const Phrase &sourcePhrase,
  TargetPhraseCollectionCache &cache) {

  typedef std::pair<size_t, size_t> AlignPoint2;

  std::vector<int> sourceWords;
  if(m_coding == REnc) {
    for(size_t i = 0; i < sourcePhrase.GetSize(); i++) {
      std::string sourceWord
        = sourcePhrase.GetWord(i).GetString(*m_input, false);
      unsigned idx = getSourceSymbolId(sourceWord);
      sourceWords.push_back(idx);
    }
  }
  
  TargetPhraseVectorPtr tpv(new TargetPhraseVector());

  unsigned phraseStopSymbol = 0;
  AlignPoint alignStopSymbol(-1, -1);
  
  std::vector<float> scores;
  std::set<AlignPoint2> alignment;
  
  enum DecodeState { New, Symbol, Score, Alignment, Add } state = New;
  
  int node = 0;
  size_t srcSize = sourcePhrase.GetSize();
  
  TargetPhrase* targetPhrase;
  
  BitStream<> encodedBitStream(encoded);
  
  while(encodedBitStream.remainingBits()) {
     
    if(state == New) {
      // Heap allocation in threads!          
      tpv->push_back(TargetPhrase(Output));
      targetPhrase = &tpv->back();
      targetPhrase->SetSourcePhrase(&sourcePhrase);
      alignment.clear();
      scores.clear();
        
      state = Symbol;
    }
    
    if(state == Symbol) {
      unsigned symbol = m_symbolTree->nextSymbol(encodedBitStream);
      
      if(symbol == phraseStopSymbol) {
        state = Score;
      }
      else {
        if(m_coding == REnc) {
          std::string wordString;
          size_t type = getREncType(symbol);
          
          if(type == 1) {
            unsigned decodedSymbol = decodeREncSymbol1(symbol);
            wordString = getTargetSymbol(decodedSymbol);
          }
          else if (type == 2) {
            size_t rank = decodeREncSymbol2Rank(symbol);
            size_t srcPos = decodeREncSymbol2Position(symbol);
            
            if(srcPos >= sourceWords.size())
              return TargetPhraseVectorPtr();  
            
            wordString = getTargetSymbol(getTranslation(sourceWords[srcPos], rank));
            if(StaticData::Instance().UseAlignmentInfo()) {
              size_t trgPos = targetPhrase->GetSize();
              alignment.insert(AlignPoint(srcPos, trgPos));
            }
          }
          else if(type == 3) {
            size_t rank = decodeREncSymbol3(symbol);
            size_t srcPos = targetPhrase->GetSize();
            
            if(srcPos >= sourceWords.size())
              return TargetPhraseVectorPtr();  
                            
            wordString = getTargetSymbol(getTranslation(sourceWords[srcPos], rank));   
            if(StaticData::Instance().UseAlignmentInfo()) {
              size_t trgPos = srcPos;
              alignment.insert(AlignPoint(srcPos, trgPos));
            }
          }
          
          Word word;
          word.CreateFromString(Output, *m_output, wordString, false);
          targetPhrase->AddWord(word);
        }
        else if(m_coding == PREnc) {
          if(getPREncType(symbol) == 1) {
            unsigned decodedSymbol = decodePREncSymbol1(symbol);
            Word word;
            word.CreateFromString(Output, *m_output,
                                  getTargetSymbol(decodedSymbol), false);
            targetPhrase->AddWord(word);
          }
          else {
            int left = decodePREncSymbol2Left(symbol);
            int right = decodePREncSymbol2Right(symbol);
            unsigned rank = decodePREncSymbol2Rank(symbol);
            
            int srcStart = left + targetPhrase->GetSize();
            int srcEnd   = srcSize - right - 1;
            
            // false positive consistency check
            if(0 > srcStart || srcStart > srcEnd || srcEnd >= srcSize) {
              return TargetPhraseVectorPtr();  
            }
            
            TargetPhraseVectorPtr subTpv = tpv;
            if(srcEnd - srcStart + 1 != srcSize) {
              Phrase subPhrase = sourcePhrase.GetSubString(WordsRange(srcStart, srcEnd));
              subTpv = m_phraseDictionary.CreateTargetPhraseCollection(subPhrase);
            }
            
            // false positive consistency check
            if(subTpv != NULL && rank < subTpv->size()) {
              TargetPhrase& subTp = subTpv->at(rank);
              targetPhrase->Append(subTp);
              if(StaticData::Instance().UseAlignmentInfo()) {
                // @TODO: dodaj srcStart oraz targetPhrase->GetSize
                //alignment.insert(subTp->GetAlignmentInfo().begin(),
                //                subTp->GetAlignmentInfo().end());
              }
            }
            else {
              return TargetPhraseVectorPtr();
            }
          }
        }
        else {
            Word word;
            word.CreateFromString(Output, *m_output,
                                  getTargetSymbol(symbol), false);
            targetPhrase->AddWord(word);
        }
      }
    }
    else if(state == Score) {
      size_t idx = m_multipleScoreTrees ? scores.size() : 0;
      float score = m_scoreTrees[idx]->nextSymbol(encodedBitStream);
      scores.push_back(score);
      
      if(scores.size() == m_numScoreComponent) {
        targetPhrase->SetScore(m_feature, scores, *m_weight, m_weightWP, *m_languageModels);
        
        if(m_containsAlignmentInfo)
          state = Alignment;
        else
          state = Add;
      }
    }
    else if(state == Alignment) {
      AlignPoint alignPoint = m_alignTree->nextSymbol(encodedBitStream);
      if(alignPoint == alignStopSymbol) {
        state = Add;
      }
      else {
        if(StaticData::Instance().UseAlignmentInfo())  
          alignment.insert(AlignPoint2(alignPoint));
      }
    }
    
    if(state == Add) {
      if(StaticData::Instance().UseAlignmentInfo())
        targetPhrase->SetAlignmentInfo(alignment);
      
      // skip over filling bits.
      if(encodedBitStream.remainingBits() <= 8)
        break;
      
      state = New;
    }    
  }
  
  if(m_coding == PREnc)
    cache.cache(sourcePhrase, tpv);
  
  return tpv;
}

}
