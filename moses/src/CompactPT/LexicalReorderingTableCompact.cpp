#include "LexicalReorderingTableCompact.h"

namespace Moses {

LexicalReorderingTableCompact::LexicalReorderingTableCompact(
  const std::string& filePath,
  const std::vector<FactorType>& f_factors,
  const std::vector<FactorType>& e_factors,
  const std::vector<FactorType>& c_factors)
  : LexicalReorderingTable(f_factors, e_factors, c_factors), m_inMemory(false), m_hash(10, 16), m_scoreTree(0)
{
  Load(filePath);
}

LexicalReorderingTableCompact::LexicalReorderingTableCompact(
  const std::vector<FactorType>& f_factors,
  const std::vector<FactorType>& e_factors,
  const std::vector<FactorType>& c_factors)
  : LexicalReorderingTable(f_factors, e_factors, c_factors), m_inMemory(false), m_hash(10, 16), m_scoreTree(0)
{}

LexicalReorderingTableCompact::~LexicalReorderingTableCompact() {
  if(m_scoreTree)
    delete m_scoreTree;
}

std::vector<float> LexicalReorderingTableCompact::GetScore(const Phrase& f,
    const Phrase& e,
    const Phrase& c)
{
  std::string key;
  Scores scores;
  size_t num_scores = 6;
  
  if(0 == c.GetSize())
    key = MakeKey(f, e, c);
  else
    for(size_t i = 0; i <= c.GetSize(); ++i)
    {
      Phrase sub_c(c.GetSubString(WordsRange(i,c.GetSize()-1)));
      key = MakeKey(f,e,sub_c);
    }
    
  size_t index = m_hash[key];
  if(m_hash.GetSize() != index)
  {
    std::string scoresString;
    if(m_inMemory)
      scoresString = m_scoresMemory[index];
    else
      scoresString = m_scoresMapped[index];
      
    BitStream<> bitStream(scoresString);
    for(size_t i = 0; i < num_scores; i++)
      scores.push_back(m_scoreTree->NextSymbol(bitStream));

    return scores;
  }

  return Scores();
}

std::string  LexicalReorderingTableCompact::MakeKey(const Phrase& f,
    const Phrase& e,
    const Phrase& c) const
{
  return MakeKey(Trim(f.GetStringRep(m_FactorsF)),
                 Trim(e.GetStringRep(m_FactorsE)),
                 Trim(c.GetStringRep(m_FactorsC)));
}

std::string  LexicalReorderingTableCompact::MakeKey(const std::string& f,
    const std::string& e,
    const std::string& c) const
{
  std::string key;
  if(!f.empty())
  {
    key += f;
  }
  if(!m_FactorsE.empty())
  {
    if(!key.empty())
    {
      key += " ||| ";
    }
    key += e;
  }
  if(!m_FactorsC.empty())
  {
    if(!key.empty())
    {
      key += " ||| ";
    }
    key += c;
  }
  key += " ||| ";
  return key;
}

void LexicalReorderingTableCompact::Load(const std::string& filePath)
{
  std::cerr << "Loading hashed version of lexical reordering model" << std::endl;
  std::string file = filePath + ".mphlexr";
  std::FILE* pFile = std::fopen(file.c_str() , "r");
  m_hash.Load(pFile);
  m_scoreTree = new CanonicalHuffman<float>(pFile);
  
  if(m_inMemory)
    m_scoresMemory.load(pFile);
  else
    m_scoresMapped.load(pFile);
  
  std::fclose(pFile);
}

}