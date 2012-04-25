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
: m_containsAlignmentInfo(true), m_symbolTree(0), m_scoreTree(0),
  m_alignTree(0), m_phraseDictionary(phraseDictionary), m_input(input),
  m_output(output), m_feature(feature), m_numScoreComponent(numScoreComponent),
  m_weight(weight), m_weightWP(weightWP),
  m_languageModels(languageModels)
{ }

PhraseDecoder::~PhraseDecoder() {
  if(m_symbolTree)
    delete m_symbolTree;
  
  if(m_scoreTree)
    delete m_scoreTree;
  
  if(m_alignTree)
    delete m_alignTree;
}

unsigned PhraseDecoder::getSourceSymbolId(std::string& symbol) {
  size_t idx = m_sourceSymbols.find(symbol);
  return idx;
}

std::string PhraseDecoder::getTargetSymbol(unsigned idx) const {
  if(idx < m_targetSymbols.size())
    return m_targetSymbols[idx];
  return std::string("##ERROR##");
}

size_t PhraseDecoder::getREncType(unsigned encodedSymbol) {
  return (encodedSymbol >> 30) + 1;
}

size_t PhraseDecoder::getPREncType(unsigned encodedSymbol) {
  return (encodedSymbol >> 31) + 1;
}

unsigned PhraseDecoder::getTranslation(unsigned srcIdx, size_t rank) {
  size_t srcTrgIdx = m_lexicalTableIndex[srcIdx];
  return m_lexicalTable[srcTrgIdx + rank].second;
}

unsigned PhraseDecoder::decodeREncSymbol1(unsigned encodedSymbol) {
  return encodedSymbol &= ~(1 << 30);
}

unsigned PhraseDecoder::decodeREncSymbol2(unsigned encodedSymbol) {
  return encodedSymbol &= ~(2 << 30);
}

unsigned PhraseDecoder::decodeREncSymbol3Rank(unsigned encodedSymbol) {
  return encodedSymbol &= ~(255 << 24);
}

unsigned PhraseDecoder::decodeREncSymbol3Position(unsigned encodedSymbol) {
  encodedSymbol &= ~(3 << 30);
  encodedSymbol >>= 24;
  return encodedSymbol;
}

unsigned PhraseDecoder::decodePREncSymbol1(unsigned encodedSymbol) {
  return encodedSymbol &= ~(1 << 31);
}

int PhraseDecoder::decodePREncSymbol2Left(unsigned encodedSymbol) {
  return ((encodedSymbol >> 25) & 63) - 32;
}

int PhraseDecoder::decodePREncSymbol2Right(unsigned encodedSymbol) {
  return ((encodedSymbol >> 19) & 63) - 32;
}

unsigned PhraseDecoder::decodePREncSymbol2Rank(unsigned encodedSymbol) {
  return (encodedSymbol & 524287);
}

size_t PhraseDecoder::load(std::FILE* in) {
  size_t start = std::ftell(in);
  
  m_targetSymbols.load(in);
  m_symbolTree = new Hufftree<int, unsigned>(in, true);
  
  {  
      size_t sum = 0, sumall = 0;
      for(size_t i = 0; i < m_symbolTree->size(); i++) {
          sumall += m_symbolTree->freq(i) * m_symbolTree->encode(m_symbolTree->data(i)).size();
          sum    += m_symbolTree->freq(i);
      }
      std::cerr << double(sumall)/sum << " bits per symbol " << sum << " " << sumall << " " << m_symbolTree->size() << std::endl;
  }
  
  m_scoreTree = new Hufftree<int, float>(in);
  
  m_alignTree = new Hufftree<int, AlignPoint>(in, true);
  
    {  
      size_t sum = 0, sumall = 0;
      for(size_t i = 0; i < m_alignTree->size(); i++) {
          sumall += m_alignTree->freq(i) * m_alignTree->encode(m_alignTree->data(i)).size();
          sum    += m_alignTree->freq(i);
      }
      std::cerr << double(sumall)/sum << " bits per align " << sum << " " << sumall << " " << m_alignTree->size() << std::endl;
  }

  size_t end = std::ftell(in);
  return end - start;
}

size_t PhraseDecoder::save(std::FILE* out) {
  size_t start = std::ftell(out);
  
  m_targetSymbols.save(out);
  m_symbolTree->Save(out);
  m_scoreTree->Save(out);
  m_alignTree->Save(out);
  
  size_t end = std::ftell(out);
  return end - start;
}
    
TargetPhraseVectorPtr PhraseDecoder::decodeCollection(
  std::string encoded, const Phrase &sourcePhrase,
  TargetPhraseCollectionCache &cache) {

  typedef std::pair<size_t, size_t> AlignPoint2;

  TargetPhraseVectorPtr tpv(new TargetPhraseVector());

  unsigned phraseStopSymbol = 0; //m_targetSymbolsMap["__SPECIAL_STOP_SYMBOL__"];
  AlignPoint alignStopSymbol(-1, -1);
  
  std::vector<float> scores;
  std::set<AlignPoint2> alignment;
  
  enum DecodeState { New, Symbol, Score, Alignment, Add } state = New;
  
  int node = 0;
  size_t srcSize = sourcePhrase.GetSize();
  
  TargetPhrase* targetPhrase;
  for(std::string::iterator it = encoded.begin(); it != encoded.end(); it++) {
    char byte = *it;
    char mask = 1;
    
    for(int i = 0; i < 8; i++) {
      
      if(state == New || state == Symbol) {
        node = (byte & mask) ? m_symbolTree->node(node+1) : m_symbolTree->node(node);
      }
      else if(state == Score) {
        node = (byte & mask) ? m_scoreTree->node(node+1) : m_scoreTree->node(node);
      }
      else if(state == Alignment) {
        if(m_containsAlignmentInfo)
          node = (byte & mask) ? m_alignTree->node(node+1) : m_alignTree->node(node);            
      }
      
      if(node < 0) {            
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
          size_t symbol = m_symbolTree->data(-node-1);
          if(symbol == phraseStopSymbol)
            state = Score;
          else {
            if(m_coding == REnc) {
              std::string wordString;
              
              size_t type = getREncType(symbol);
              if(type == 1) {
                unsigned decodedSymbol = decodeREncSymbol1(symbol);
                wordString = getTargetSymbol(decodedSymbol);
              }
              else if(type == 2) {
                size_t rank = decodeREncSymbol2(symbol);
                size_t srcPos = targetPhrase->GetSize();
                size_t trgPos = srcPos;
                                
                std::string sourceWord
                  = sourcePhrase.GetWord(srcPos).GetString(*m_input, false);
                
                unsigned idx = getSourceSymbolId(sourceWord);
                wordString = getTargetSymbol(getTranslation(idx, rank));
                    
                //if(StaticData::Instance().UseAlignmentInfo())
                //  alignment.insert(AlignPoint(srcPos, trgPos));
              }
              else if (type == 3) {
                size_t rank = decodeREncSymbol3Rank(symbol);
                size_t srcPos = decodeREncSymbol3Position(symbol);
                size_t trgPos = targetPhrase->GetSize();
                
                std::string sourceWord
                  = sourcePhrase.GetWord(srcPos).GetString(*m_input, false);
                    
                unsigned idx = getSourceSymbolId(sourceWord);
                wordString = getTargetSymbol(getTranslation(idx, rank));
                 
                //if(StaticData::Instance().UseAlignmentInfo())
                //  alignment.insert(AlignPoint(srcPos, trgPos));
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
                  // TODO dodaj srcStart oraz targetPhrase->GetSize
                  //alignment.insert(subTp->GetAlignmentInfo().begin(),
                  //                subTp->GetAlignmentInfo().end());
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
          float score = m_scoreTree->data(-node-1);
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
          AlignPoint alignPoint = m_alignTree->data(-node-1);
          if(alignPoint == alignStopSymbol) {
            state = Add;
          }
          else {
            //alignment.insert(AlignPoint2(alignPoint));
          }
        }
        
        if(state == Add) {
          //targetPhrase->SetAlignmentInfo(alignment);
          //phraseColl->Add(targetPhrase);
          
          if(std::distance(it, encoded.end()) == 1)
            break;
          state = New;
        }
        
        node = 0;
      }
      mask = mask << 1;
    }
  }
  cache.cache(sourcePhrase, tpv);
  return tpv;
}

}
