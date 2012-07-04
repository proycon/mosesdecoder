// $Id: PhraseDictionaryMemoryHashed.cpp 3908 2011-02-28 11:41:08Z pjwilliams $
// vim:tabstop=2

/***********************************************************************
Moses - factored phrase-based language decoder
Copyright (C) 2006 University of Edinburgh

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
***********************************************************************/

#include <fstream>
#include <string>
#include <iterator>
#include <queue>
#include <algorithm>
#include <sys/stat.h>

#include "PhraseDictionaryMemoryHashed.h"
#include "FactorCollection.h"
#include "Word.h"
#include "Util.h"
#include "InputFileStream.h"
#include "StaticData.h"
#include "WordsRange.h"
#include "UserMessage.h"
#include "ThreadPool.h"

using namespace std;

namespace Moses
{
  
bool PhraseDictionaryMemoryHashed::Load(const std::vector<FactorType> &input
                                  , const std::vector<FactorType> &output
                                  , const string &filePath
                                  , const vector<float> &weight
                                  , size_t tableLimit
                                  , const LMList &languageModels
                                  , float weightWP)
{
  m_input = &input;
  m_output = &output;
  m_weight = &weight;
  m_tableLimit = tableLimit;
  m_languageModels = &languageModels; 
  m_weightWP = weightWP;

  if(m_implementation == MemoryHashedBinary)
    return LoadBinary(filePath);
  else
    return LoadBinary(filePath);
    
  return false;
}

bool PhraseDictionaryMemoryHashed::LoadBinary(std::string filePath) {
    if (FileExists(filePath + ".mph"))
        filePath += ".mph";

    m_phraseDecoder = new PhraseDecoder(*this, m_input, m_output, m_feature,
                                    m_numScoreComponent, m_weight, m_weightWP,
                                    m_languageModels);
  
    std::FILE* pFile = std::fopen(filePath.c_str() , "r");
    m_hash.LoadIndex(pFile);
    //size_t hashSize = m_hash.Load(pFile);
    //std::cerr << "Total HashIndex size: " << float(hashSize)/(1024*1024) << " M" << std::endl;

    size_t coderSize = m_phraseDecoder->load(pFile);
    std::cerr << "Total PhraseCoder size: " << float(coderSize)/(1024*1024) << " M" << std::endl;

    size_t phraseSize = m_targetPhrases.load(pFile, true);
    std::cerr << "Total TargetPhrases size: " << float(phraseSize)/(1024*1024) << " M" << std::endl;
    //std::fclose(pFile);
    
    return coderSize && phraseSize;    
    //return hashSize && coderSize && phraseSize;
}

std::string PhraseDictionaryMemoryHashed::makeSourceKey(std::string &source) {
    return source + m_phraseDecoder->getSeparator();
}

TargetPhraseVectorPtr
PhraseDictionaryMemoryHashed::CreateTargetPhraseCollection(const Phrase
                                                           &sourcePhrase) {

  //TargetPhraseVectorPtr tpv = m_decodingCache.retrieve(sourcePhrase);
  //if(tpv != NULL)
  //  return tpv;

  // Retrieve source phrase identifier
  std::string sourcePhraseString = sourcePhrase.GetStringRep(*m_input);
  size_t sourcePhraseId = m_hash[makeSourceKey(sourcePhraseString)];
  
  if(sourcePhraseId != m_hash.GetSize()) {    
    // Retrieve compressed and encoded target phrase collection
    std::string encodedPhraseCollection = m_targetPhrases[sourcePhraseId];

    // Decompress and decode target phrase collection
    TargetPhraseVectorPtr decodedPhraseColl =
      m_phraseDecoder->decodeCollection(encodedPhraseCollection, sourcePhrase, m_decodingCache);
    
    return decodedPhraseColl;
  }
  else
    return TargetPhraseVectorPtr();
}

struct CompareTargetPhrase {
  bool operator() (const TargetPhrase &a, const TargetPhrase &b) {
    return a.GetFutureScore() > b.GetFutureScore();
  }
};

const TargetPhraseCollection*
PhraseDictionaryMemoryHashed::GetTargetPhraseCollection(const Phrase &sourcePhrase) const {
  
  // There is no souch source phrase if source phrase is longer than longest
  // observed source phrase during compilation 
  if(sourcePhrase.GetSize() > m_phraseDecoder->getMaxSourcePhraseLength())
    return NULL;
  
  // Only for PREnc: Check whether phrase pair has been created previously
  TargetPhraseCollection* cachedPhraseColl
    = const_cast<PhraseDictionaryMemoryHashed*>(this)->RetrieveFromCache(sourcePhrase);
  if(cachedPhraseColl != NULL)
    return cachedPhraseColl;
  
  // Retrieve target phrase collection from phrase table
  TargetPhraseVectorPtr decodedPhraseColl
    = const_cast<PhraseDictionaryMemoryHashed*>(this)->CreateTargetPhraseCollection(sourcePhrase);
  
  if(decodedPhraseColl != NULL && decodedPhraseColl->size()) {
    TargetPhraseVectorPtr tpv(new TargetPhraseVector(*decodedPhraseColl));
    TargetPhraseCollection* phraseColl = new TargetPhraseCollection();
    
    // Score phrases and if possible apply ttable_limit
    TargetPhraseVector::iterator nth =
      (m_tableLimit == 0 || tpv->size() < m_tableLimit) ?
      tpv->end() : tpv->begin() + m_tableLimit;
    std::nth_element(tpv->begin(), nth, tpv->end(), CompareTargetPhrase());
    for(TargetPhraseVector::iterator it = tpv->begin(); it != nth; it++)
      phraseColl->Add(new TargetPhrase(*it));
    
    // Cache phrase pair for for clean-up or retrieval with PREnc
    const_cast<PhraseDictionaryMemoryHashed*>(this)->CacheForCleanup(sourcePhrase, phraseColl);
    
    return phraseColl;
  }
  else
    return NULL;
  
}

PhraseDictionaryMemoryHashed::~PhraseDictionaryMemoryHashed() {
  if(m_phraseDecoder)
    delete m_phraseDecoder;
    
  CleanUp();
}

//TO_STRING_BODY(PhraseDictionaryMemoryHashed)

TargetPhraseCollection*
PhraseDictionaryMemoryHashed::RetrieveFromCache(const Phrase &sourcePhrase) {
#ifdef WITH_THREADS
  boost::mutex::scoped_lock lock(m_sentenceMutex);
  PhraseCache &ref = m_sentenceCache[pthread_self()]; 
#else
  PhraseCache &ref = m_sentenceCache; 
#endif
  PhraseCache::iterator it = ref.find(sourcePhrase);
  if(it != ref.end())
    return it->second;
  else
    return NULL;
}

void PhraseDictionaryMemoryHashed::CacheForCleanup(const Phrase &sourcePhrase,
                                                   TargetPhraseCollection* tpc) {
#ifdef WITH_THREADS
  boost::mutex::scoped_lock lock(m_sentenceMutex);
  m_sentenceCache[pthread_self()].insert(std::make_pair(sourcePhrase, tpc));
#else
  m_sentenceCache.insert(std::make_pair(sourcePhrase, tpc));
#endif
}

void PhraseDictionaryMemoryHashed::InitializeForInput(const Moses::InputType&) {}

void
PhraseDictionaryMemoryHashed::AddEquivPhrase(const Phrase &source,
                                             const TargetPhrase &targetPhrase) { }

void PhraseDictionaryMemoryHashed::CleanUp() {
  m_hash.KeepNLastRanges(0.01, 0.2);
  m_decodingCache.prune();
  
#ifdef WITH_THREADS
  boost::mutex::scoped_lock lock(m_sentenceMutex);
  PhraseCache &ref = m_sentenceCache[pthread_self()]; 
#else
  PhraseCache &ref = m_sentenceCache; 
#endif
  
  for(PhraseCache::iterator it = ref.begin(); it != ref.end(); it++) 
      delete it->second;
      
  PhraseCache temp;
  temp.swap(ref);
}

}

