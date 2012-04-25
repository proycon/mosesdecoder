#ifndef moses_PhraseCache_h
#define moses_PhraseCache_h

#include <map>
#include <set>

#include "Phrase.h"
#include "TargetPhraseCollection.h"

class TargetPhraseCollectionCache {
  private:
    size_t m_max;
    float m_tolerance;
    
    typedef std::pair<clock_t, TargetPhraseCollection*> LastUsed;
    typedef std::map<Phrase, LastUsed> CacheMap;
    
    CacheMap m_phraseCache;
    
  public:
    
    TargetPhraseCollectionCache(size_t max, float tolerance)
    : m_max(max), m_tolerance(tolerance)
    {}
    
    void CacheTargetPhraseCollection(const Phrase &sourcePhrase,
                                     TargetPhraseCollection *tpc) {
      if(m_phraseCache.count(sourcePhrase) &&
         m_phraseCache[sourcePhrase].second != tpc) {
        std::cerr << "Alert!" << std::endl;
      }
      m_phraseCache[sourcePhrase] = std::make_pair(clock(), tpc);
    }

    TargetPhraseCollection*
    CachedTargetPhraseCollection(const Phrase &sourcePhrase) {
      if(m_phraseCache.count(sourcePhrase)) {
        LastUsed &lu = m_phraseCache[sourcePhrase];
        lu.first = clock();
        return lu.second;
      }
      return NULL;
    }

    void Prune() {
      if(m_phraseCache.size() > m_max * (1 + m_tolerance)) {
        typedef std::set<std::pair<clock_t, Phrase> > Cands;
        Cands cands; 
        for(CacheMap::iterator it = m_phraseCache.begin();
            it != m_phraseCache.end(); it++) {
          LastUsed &lu = it->second;
          cands.insert(std::make_pair(lu.first, it->first));
        }
         
        for(Cands::iterator it = cands.begin(); it != cands.end(); it++) {
          const Phrase& p = it->second;
          LastUsed lu = tCache[p];
           
          delete lu.second;
          tCache.erase(p);
          if(tCache.size() < (m_max * (1 - m_tolerance)))
            break;
        }
      } 
    }

    void CleanUp() {
      for(CacheMap::iterator it = m_phraseCache.begin(); it != m_phraseCache.end(); it++) {
        LastUsed &lu = it->second;
        delete lu.second;
      }
      m_phraseCache.clear();
    }
    
};

#endif