#ifndef moses_PhraseCache_h
#define moses_PhraseCache_h

#include <map>
#include <set>
#include <vector>

#ifdef WITH_THREADS
#ifdef BOOST_HAS_PTHREADS
#include <boost/thread/mutex.hpp>
#endif
#endif

#include <boost/shared_ptr.hpp>

#include "Phrase.h"
#include "TargetPhraseCollection.h"

namespace Moses {

typedef std::vector<TargetPhrase> TargetPhraseVector;
typedef boost::shared_ptr<TargetPhraseVector> TargetPhraseVectorPtr;

class TargetPhraseCollectionCache {
  private:
    size_t m_max;
    float m_tolerance;
    bool m_destroy;
    
    typedef std::pair<clock_t, TargetPhraseVectorPtr> LastUsed;
    typedef std::map<Phrase, LastUsed> CacheMap;
    
    CacheMap m_phraseCache;
    
    boost::mutex m_mutex;
    
  public:
    
    typedef CacheMap::iterator iterator;
    typedef CacheMap::const_iterator const_iterator;
        
    TargetPhraseCollectionCache(size_t max = 5000, float tolerance = 0.2, bool destroy = true)
    : m_max(max), m_tolerance(tolerance), m_destroy(destroy)
    {}
    
    ~TargetPhraseCollectionCache() {
      if(m_destroy)
        cleanUp();
    }
    
    void setDestroy(bool destroy) {
      m_destroy = destroy;
    }
    
    iterator begin() {
      return m_phraseCache.begin();
    }
    
    const_iterator begin() const {
      return m_phraseCache.begin();
    }
    
    iterator end() {
      return m_phraseCache.end();
    }
    
    const_iterator end() const {
      return m_phraseCache.end();
    }
    
    void cache(const Phrase &sourcePhrase, TargetPhraseVectorPtr tpc) {
      boost::mutex::scoped_lock lock(m_mutex);

      iterator it = m_phraseCache.find(sourcePhrase);
      if(it != m_phraseCache.end()) {
        it->second.first = clock();
      }
      else
        m_phraseCache[sourcePhrase] = std::make_pair(clock(), tpc);
    }

    TargetPhraseVectorPtr retrieve(const Phrase &sourcePhrase) {
      boost::mutex::scoped_lock lock(m_mutex);

      iterator it = m_phraseCache.find(sourcePhrase);
      if(it != m_phraseCache.end()) {
        LastUsed &lu = it->second;
        lu.first = clock();
        return lu.second;
      }
      else
        return TargetPhraseVectorPtr();
    }

    void prune() {
      boost::mutex::scoped_lock lock(m_mutex);

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
          m_phraseCache.erase(p);
          
          if(m_phraseCache.size() < (m_max * (1 - m_tolerance)))
            break;
        }
      }
    }

    void cleanUp() {
      boost::mutex::scoped_lock lock(m_mutex);
      m_phraseCache.clear();
    }
    
};

}

#endif