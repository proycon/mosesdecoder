#ifndef moses_TargetPhraseCollectionCache_h
#define moses_TargetPhraseCollectionCache_h

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

// Avoid using new due to locking
typedef std::vector<TargetPhrase> TargetPhraseVector;
typedef boost::shared_ptr<TargetPhraseVector> TargetPhraseVectorPtr;

class TargetPhraseCollectionCache {
  private:
    size_t m_max;
    float m_tolerance;
    
    struct LastUsed {
      clock_t m_clock;
      TargetPhraseVectorPtr m_tpv;
      size_t m_bitsLeft;
      
      LastUsed()
      : m_clock(0), m_bitsLeft(0) {}
      
      LastUsed(clock_t clock, TargetPhraseVectorPtr tpv, size_t bitsLeft = 0)
      : m_clock(clock), m_tpv(tpv), m_bitsLeft(bitsLeft) {}
    };
    
    typedef std::map<Phrase, LastUsed> CacheMap;
    
    CacheMap m_phraseCache;
    
#ifdef WITH_THREADS
    boost::mutex m_mutex;
#endif

  public:
    
    typedef CacheMap::iterator iterator;
    typedef CacheMap::const_iterator const_iterator;
        
    TargetPhraseCollectionCache(size_t max = 5000, float tolerance = 0.2)
    : m_max(max), m_tolerance(tolerance)
    {}
    
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
    
    void cache(const Phrase &sourcePhrase, TargetPhraseVectorPtr tpv, size_t bitsLeft = 0) {
#ifdef WITH_THREADS
      boost::mutex::scoped_lock lock(m_mutex);
#endif

      iterator it = m_phraseCache.find(sourcePhrase);
      if(it != m_phraseCache.end())
        it->second.m_clock = clock();
      else
        m_phraseCache[sourcePhrase] = LastUsed(clock(), tpv, bitsLeft);
    }

    std::pair<TargetPhraseVectorPtr, size_t> retrieve(const Phrase &sourcePhrase) {
#ifdef WITH_THREADS
      boost::mutex::scoped_lock lock(m_mutex);
#endif

      iterator it = m_phraseCache.find(sourcePhrase);
      if(it != m_phraseCache.end()) {
        LastUsed &lu = it->second;
        lu.m_clock = clock();
        return std::make_pair(lu.m_tpv, lu.m_bitsLeft);
      }
      else
        return std::make_pair(TargetPhraseVectorPtr(), 0);
    }

    void prune() {
#ifdef WITH_THREADS
      boost::mutex::scoped_lock lock(m_mutex);
#endif

      if(m_phraseCache.size() > m_max * (1 + m_tolerance)) {
        typedef std::set<std::pair<clock_t, Phrase> > Cands;
        Cands cands; 
        for(CacheMap::iterator it = m_phraseCache.begin();
            it != m_phraseCache.end(); it++) {
          LastUsed &lu = it->second;
          cands.insert(std::make_pair(lu.m_clock, it->first));
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
#ifdef WITH_THREADS
      boost::mutex::scoped_lock lock(m_mutex);
#endif
      m_phraseCache.clear();
    }
    
};

}

#endif