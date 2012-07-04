#ifndef moses_PhrasetableCreator_h
#define moses_PhrasetableCreator_h

#include <sstream>
#include <iostream>
#include <queue>
#include <vector>
#include <boost/unordered_map.hpp>

#include "InputFileStream.h"
#include "ThreadPool.h"
#include "UserMessage.h"
#include "Util.h"

#include "CompactPT/BlockHashIndex.h"
#include "CompactPT/ConsistantPhrases.h"
#include "CompactPT/StringVector.h"
#include "CompactPT/CanonicalHuffman.h"

// @TODO: Check speed for multithread, locking etc.
// Add quantization

namespace Moses {
    
typedef std::pair<unsigned char, unsigned char> AlignPoint;

template <typename DataType> 
class Counter {
  public:
    typedef boost::unordered_map<DataType, size_t> FreqMap;
    typedef typename FreqMap::iterator iterator;
    typedef typename FreqMap::mapped_type mapped_type;
    typedef typename FreqMap::value_type value_type; 

  private:
#ifdef WITH_THREADS    
    boost::mutex m_mutex;
#endif
    FreqMap m_freqMap;
    size_t m_maxSize;
    std::vector<DataType> m_bestVec;
    
    struct FreqSorter {
        bool operator()(const value_type& a, const value_type& b) const {
            if(a.second > b.second)
                return true;
            // Check impact on translation quality!
            if(a.second == b.second && a.first > b.first)
                return true;
            return false;
        }
    };
  
  public:
    Counter() : m_maxSize(0) {}
    
    iterator begin() {
        return m_freqMap.begin();
    }
    
    iterator end() {
        return m_freqMap.end();
    }
    
    void increase(DataType data) {
#ifdef WITH_THREADS
        boost::mutex::scoped_lock lock(m_mutex);
#endif
        m_freqMap[data]++;
    }
  
    void increaseBy(DataType data, size_t num) {
#ifdef WITH_THREADS
        boost::mutex::scoped_lock lock(m_mutex);
#endif
        m_freqMap[data] += num;
    }
    
    mapped_type operator[](DataType data) {
        return m_freqMap[data];
    }
    
    size_t size() {
#ifdef WITH_THREADS
        boost::mutex::scoped_lock lock(m_mutex);
#endif
        return m_freqMap.size();
    }
    
    void quantize(size_t maxSize) {
#ifdef WITH_THREADS
        boost::mutex::scoped_lock lock(m_mutex);
#endif
        m_maxSize = maxSize;
        std::vector<std::pair<DataType, mapped_type> > freqVec;
        freqVec.insert(freqVec.begin(), m_freqMap.begin(), m_freqMap.end());
        std::sort(freqVec.begin(), freqVec.end(), FreqSorter());
        
        for(size_t i = 0; i < freqVec.size() && i < m_maxSize; i++)
            m_bestVec.push_back(freqVec[i].first);
          
        std::sort(m_bestVec.begin(), m_bestVec.end());
        
        FreqMap t_freqMap;
        for(typename std::vector<std::pair<DataType, mapped_type> >::iterator it
            = freqVec.begin(); it != freqVec.end(); it++) {
            DataType closest = lowerBound(it->first);
            t_freqMap[closest] += it->second;
        }
        
        m_freqMap.swap(t_freqMap);   
    }
    
    void clear() {
#ifdef WITH_THREADS
        boost::mutex::scoped_lock lock(m_mutex);
#endif
        m_freqMap.clear();
    }
    
    DataType lowerBound(DataType data) {
        if(m_maxSize == 0 || m_bestVec.size() == 0)
            return data;
        else {
            typename std::vector<DataType>::iterator it
              = std::lower_bound(m_bestVec.begin(), m_bestVec.end(), data);
            if(it != m_bestVec.end())
                return *it;
            else
                return m_bestVec.back();
        }
    }
};
 
class PackedItem {
  private:
    long m_line;
    std::string m_sourcePhrase;
    std::string m_packedTargetPhrase;
    size_t m_rank;
    float m_score;
    
public:
    PackedItem(long line, std::string sourcePhrase,
               std::string packedTargetPhrase, size_t rank,
               float m_score = 0);
    
    long getLine() const;
    const std::string& getSrc() const;
    const std::string& getTrg() const;
    size_t getRank() const;
    float getScore() const;
};

static bool operator<(const PackedItem &pi1, const PackedItem &pi2);

class PhrasetableCreator {
  public:
    enum Coding { None, REnc, PREnc };
    
  private:
    std::string m_inPath;
    std::string m_outPath;
    
    std::FILE* m_outFile;
    
    size_t m_numScoreComponent;
    Coding m_coding;
    size_t m_orderBits;
    size_t m_fingerPrintBits;
    bool m_useAlignmentInfo;
    bool m_multipleScoreTrees;
    size_t m_quantize;
    size_t m_maxRank;
        
    size_t m_maxPhraseLength;
    
    static std::string m_phraseStopSymbol;
    static std::string m_separator;
    
#ifdef WITH_THREADS
    size_t m_threads;
    boost::mutex m_mutex;
#endif
    
    BlockHashIndex m_srcHash;
    BlockHashIndex m_rnkHash;
    
    std::vector<unsigned> m_ranks;
    
    typedef std::pair<unsigned, unsigned> SrcTrg;
    typedef std::pair<std::string, std::string> SrcTrgString;
    typedef std::pair<SrcTrgString, float> SrcTrgProb;
    
        
    struct SrcTrgProbSorter {
      bool operator()(const SrcTrgProb& a, const SrcTrgProb& b) const {
        
        if(a.first.first < b.first.first)
          return true;
        
        // Check impact on translation quality!
        if(a.first.first == b.first.first && a.second > b.second)
          return true;
        
        if(a.first.first == b.first.first
           && a.second == b.second
           && a.first.second < b.first.second
        )
          return true;
        
        return false;
      }
    };
    
    std::vector<size_t> m_lexicalTableIndex;
    std::vector<SrcTrg> m_lexicalTable;
    
    StringVector<unsigned char, unsigned long, MmapAllocator>
        m_encodedTargetPhrases;
        
    StringVector<unsigned char, unsigned long, MmapAllocator>
        m_compressedTargetPhrases;
    
    boost::unordered_map<std::string, unsigned> m_targetSymbolsMap;
    boost::unordered_map<std::string, unsigned> m_sourceSymbolsMap;
    
    typedef Counter<unsigned> SymbolCounter;
    typedef Counter<float> ScoreCounter;
    typedef Counter<AlignPoint> AlignCounter;
    
    typedef CanonicalHuffman<unsigned> SymbolTree;
    typedef CanonicalHuffman<float> ScoreTree;
    typedef CanonicalHuffman<AlignPoint> AlignTree;
    
    SymbolCounter m_symbolCounter;
    SymbolTree* m_symbolTree;
    
    AlignCounter m_alignCounter;
    AlignTree* m_alignTree; 
    
    std::vector<ScoreCounter*> m_scoreCounters;
    std::vector<ScoreTree*> m_scoreTrees;
    
    std::priority_queue<PackedItem> m_queue;
    long m_lastFlushedLine;
    long m_lastFlushedSourceNum;
    std::string m_lastFlushedSourcePhrase;
    std::vector<std::string> m_lastSourceRange;
    std::priority_queue<std::pair<float, size_t> > m_rankQueue;
    std::vector<std::string> m_lastCollection;
    
    void save();
    void printInfo();
    
    void addSourceSymbolId(std::string& symbol);
    unsigned getSourceSymbolId(std::string& symbol);
    
    void addTargetSymbolId(std::string& symbol);
    unsigned getTargetSymbolId(std::string& symbol);
    unsigned getOrAddTargetSymbolId(std::string& symbol);
    
    unsigned getRank(unsigned srcIdx, unsigned trgIdx);
    
    unsigned encodeREncSymbol1(unsigned symbol);
    unsigned encodeREncSymbol2(unsigned position, unsigned rank);
    unsigned encodeREncSymbol3(unsigned rank);
    
    unsigned encodePREncSymbol1(unsigned symbol);
    unsigned encodePREncSymbol2(int lOff, int rOff, unsigned rank);
    
    void encodeTargetPhraseNone(std::vector<std::string>& t,
                                std::ostream& os);
    
    void encodeTargetPhraseREnc(std::vector<std::string>& s,
                                std::vector<std::string>& t,
                                std::set<AlignPoint>& a,
                                std::ostream& os);
    
    void encodeTargetPhrasePREnc(std::vector<std::string>& s,
                                 std::vector<std::string>& t,
                                 std::set<AlignPoint>& a, size_t ownRank,
                                 std::ostream& os);
    void encodeScores(std::vector<float>& scores, std::ostream& os);
    void encodeAlignment(std::set<AlignPoint>& alignment, std::ostream& os);
    
    std::string makeSourceKey(std::string&);
    std::string makeSourceTargetKey(std::string&, std::string&);
    
    void loadLexicalTable(std::string filePath);
    
    void createRankHash();
    void encodeTargetPhrases();
    void calcHuffmanCodes();
    void compressTargetPhrases();
    
    void addRankedLine(PackedItem& pi);
    void flushRankedQueue(bool force = false);
    
    std::string encodeLine(std::vector<std::string>& tokens, size_t ownRank);
    void addEncodedLine(PackedItem& pi);
    void flushEncodedQueue(bool force = false);
    
    std::string compressEncodedCollection(std::string encodedCollection);
    void addCompressedCollection(PackedItem& pi);
    void flushCompressedQueue(bool force = false);
    
  public:
    
    PhrasetableCreator(std::string inPath,
                       std::string outPath,
                       size_t numScoreComponent = 5,
                       Coding coding = PREnc,
                       size_t orderBits = 10,
                       size_t fingerPrintBits = 16,
                       bool useAlignmentInfo = false,
                       bool multipleScoreTrees = true,
                       size_t quantize = 0,
                       size_t maxRank = 100
#ifdef WITH_THREADS
                       , size_t threads = 2
#endif
                      );
    
    friend class RankingTask;
    friend class EncodingTask;
    friend class CompressionTask;
};

class RankingTask {
  private:
#ifdef WITH_THREADS
    static boost::mutex m_mutex;
    static boost::mutex m_fileMutex;
#endif
    static size_t m_lineNum;
    InputFileStream& m_inFile;
    PhrasetableCreator& m_creator;
    
  public:
    RankingTask(InputFileStream& inFile, PhrasetableCreator& creator);
    void operator()();
};

class EncodingTask {
  private:
#ifdef WITH_THREADS
    static boost::mutex m_mutex;
    static boost::mutex m_fileMutex;
#endif
    static size_t m_lineNum;
    static size_t m_sourcePhraseNum;
    static std::string m_lastSourcePhrase;
    
    InputFileStream& m_inFile;
    PhrasetableCreator& m_creator;
    
  public:
    EncodingTask(InputFileStream& inFile, PhrasetableCreator& creator);
    void operator()();
};

class CompressionTask {
  private:
#ifdef WITH_THREADS
    static boost::mutex m_mutex;
#endif
    static size_t m_collectionNum;
    StringVector<unsigned char, unsigned long, MmapAllocator>&
        m_encodedCollections;
    PhrasetableCreator& m_creator;
    
  public:
    CompressionTask(StringVector<unsigned char, unsigned long, MmapAllocator>&
                    encodedCollections, PhrasetableCreator& creator);
    void operator()();
};


}

#endif