#include <cstdio>

#include "PhrasetableCreator.h"

namespace Moses {
    
std::string PhrasetableCreator::m_phraseStopSymbol = "__SPECIAL_STOP_SYMBOL__";
std::string PhrasetableCreator::m_separator = " ||| ";
    
PhrasetableCreator::PhrasetableCreator(std::string inPath,
                                       std::string outPath,
                                       size_t numScoreComponent,
                                       Coding coding,
                                       size_t orderBits,
                                       size_t fingerPrintBits,
                                       bool useAlignmentInfo,
                                       bool multipleScoreTrees,
                                       size_t quantize
#ifdef WITH_THREADS
                                       , size_t threads
#endif
                                       )
  : m_inPath(inPath), m_outPath(outPath),
    m_outFile(std::fopen(m_outPath.c_str(), "w")), m_numScoreComponent(numScoreComponent),
    m_coding(coding), m_orderBits(orderBits), m_fingerPrintBits(fingerPrintBits),
    m_useAlignmentInfo(useAlignmentInfo),
    m_multipleScoreTrees(multipleScoreTrees),
    m_quantize(quantize),
#ifdef WITH_THREADS
    m_threads(threads),
    m_srcHash(m_orderBits, m_fingerPrintBits, 1),
    m_rnkHash(10, 24, m_threads),
#else
    m_srcHash(m_orderBits, m_fingerPrintBits),
    m_rnkHash(m_orderBits, m_fingerPrintBits),
#endif
    m_maxPhraseLength(0),
    m_lastFlushedLine(-1), m_lastFlushedSourceNum(0),
    m_lastFlushedSourcePhrase("")
{
    
    addTargetSymbolId(m_phraseStopSymbol);
    
    m_scoreCounters.resize(m_multipleScoreTrees ? m_numScoreComponent : 1);
    for(std::vector<ScoreCounter*>::iterator it = m_scoreCounters.begin();
        it != m_scoreCounters.end(); it++)
        *it = new ScoreCounter();
    m_scoreTrees.resize(m_multipleScoreTrees ? m_numScoreComponent : 1);
    
    // 0th pass
    if(m_coding == REnc) {
        size_t found = inPath.find_last_of("/\\");
        std::string path = inPath.substr(0, found);
        loadLexicalTable(path + "/lex.f2e");
    }
    else if(m_coding == PREnc) {
        createRankHash();
    }
    
    // 1st pass
    m_srcHash.BeginSave(m_outFile);   
    encodeTargetPhrases();
    
    std::cerr << "EncodedPhrases: " << m_encodedTargetPhrases.size() << std::endl;
        
    std::cerr << "SymbolCounter: " << m_symbolCounter.size() << std::endl;
    for(std::vector<ScoreCounter*>::iterator it = m_scoreCounters.begin();
        it != m_scoreCounters.end(); it++)
        std::cerr << "ScoreCounter: " << (*it)->size() << std::endl;
    std::cerr << "AlignCounter: " << m_alignCounter.size() << std::endl;
    
    calcHuffmanCodes();
    
    // 2nd pass
    compressTargetPhrases();
    
    std::cerr << "CompressedPhrases: " << m_compressedTargetPhrases.size() << std::endl;

    save();
    std::fclose(m_outFile);
}
    
void PhrasetableCreator::save() {
    // Save type of encoding
    std::fwrite(&m_coding, sizeof(m_coding), 1, m_outFile);
    std::fwrite(&m_numScoreComponent, sizeof(m_numScoreComponent), 1, m_outFile);
    std::fwrite(&m_useAlignmentInfo, sizeof(m_useAlignmentInfo), 1, m_outFile);
    std::fwrite(&m_maxPhraseLength, sizeof(m_maxPhraseLength), 1, m_outFile);
    
    if(m_coding == REnc) {
        // Save source language symbols for REnc
        std::vector<std::string> temp1;
        temp1.resize(m_sourceSymbolsMap.size());
        for(boost::unordered_map<std::string, unsigned>::iterator it
            = m_sourceSymbolsMap.begin(); it != m_sourceSymbolsMap.end(); it++)
            temp1[it->second] = it->first;
        std::sort(temp1.begin(), temp1.end());
        StringVector<unsigned char, unsigned, std::allocator> sourceSymbols;
        for(std::vector<std::string>::iterator it = temp1.begin();
            it != temp1.end(); it++)
            sourceSymbols.push_back(*it);
        sourceSymbols.save(m_outFile);
        
        // Save lexical translation table for REnc
        size_t size = m_lexicalTableIndex.size();
        std::fwrite(&size, sizeof(size_t), 1, m_outFile);
        std::fwrite(&m_lexicalTableIndex[0], sizeof(size_t), size, m_outFile);
        size = m_lexicalTable.size();
        std::fwrite(&size, sizeof(size_t), 1, m_outFile);
        std::fwrite(&m_lexicalTable[0], sizeof(SrcTrg), size, m_outFile);
    }
    
    // Save target language symbols
    std::vector<std::string> temp2;
    temp2.resize(m_targetSymbolsMap.size());
    for(boost::unordered_map<std::string, unsigned>::iterator it
        = m_targetSymbolsMap.begin(); it != m_targetSymbolsMap.end(); it++)
        temp2[it->second] = it->first;
    StringVector<unsigned char, unsigned, std::allocator> targetSymbols;
    for(std::vector<std::string>::iterator it = temp2.begin();
        it != temp2.end(); it++)
        targetSymbols.push_back(*it);
    targetSymbols.save(m_outFile);
    
    // Save Huffman codes for target language symbols
    m_symbolTree->save(m_outFile);
    
    // Save number of Huffman code sets for scores and
    // save Huffman code sets
    std::fwrite(&m_multipleScoreTrees, sizeof(m_multipleScoreTrees), 1, m_outFile);
    size_t numScoreTrees = m_scoreTrees.size();
    for(size_t i = 0; i < numScoreTrees; i++)
        m_scoreTrees[i]->save(m_outFile);
    
    // Save Huffman codes for alignments
    if(m_useAlignmentInfo)
        m_alignTree->save(m_outFile);
    
    // Save compressed target phrase collections 
    size_t phraseSize = m_compressedTargetPhrases.save(m_outFile);
}
    
void PhrasetableCreator::loadLexicalTable(std::string filePath) {
    std::vector<SrcTrgProb> t_lexTable;
      
    std::cerr << "Reading in lexical table from " << filePath << std::endl;
    std::ifstream lexIn(filePath.c_str(), std::ifstream::in);
    std::string src, trg;
    float prob;
    
    // Reading in the translation probability lexicon
    
    while(lexIn >> trg >> src >> prob) {
        if(t_lexTable.size() % 10000 == 0)
            std::cerr << ".";
        t_lexTable.push_back(SrcTrgProb(SrcTrgString(src, trg), prob));
        addSourceSymbolId(src);
        addTargetSymbolId(trg);
    }
    std::cerr << std::endl;
    
    // Sorting lexicon by source words by lexicographical order, corresponding
    // target words by decreasing probability.
    
    std::cerr << "Read in " << t_lexTable.size() << " lexical pairs" << std::endl;
    std::sort(t_lexTable.begin(), t_lexTable.end(), SrcTrgProbSorter());
    std::cerr << "Sorted" << std::endl;
    
    // Re-assigning source word ids in lexicographical order
    
    std::vector<std::string> temp1;
    temp1.resize(m_sourceSymbolsMap.size());
    for(boost::unordered_map<std::string, unsigned>::iterator it
        = m_sourceSymbolsMap.begin(); it != m_sourceSymbolsMap.end(); it++)
        temp1[it->second] = it->first;
        
    std::sort(temp1.begin(), temp1.end());
    
    for(size_t i = 0; i < temp1.size(); i++)
        m_sourceSymbolsMap[temp1[i]] = i;
    
    // Building the lexicon based on source and target word ids
    
    std::string srcWord = "";
    size_t srcIdx = 0;
    for(std::vector<SrcTrgProb>::iterator it = t_lexTable.begin();
        it != t_lexTable.end(); it++) {
      
        // If we encounter a new source word
        if(it->first.first != srcWord) {
            srcIdx = getSourceSymbolId(it->first.first);
            
            // Store position of first translation
            if(srcIdx >= m_lexicalTableIndex.size())
                m_lexicalTableIndex.resize(srcIdx + 1);
            m_lexicalTableIndex[srcIdx] = m_lexicalTable.size();
        }
        
        // Store pair of source word and target word
        size_t trgIdx = getTargetSymbolId(it->first.second);        
        m_lexicalTable.push_back(SrcTrg(srcIdx, trgIdx));
        
        if(m_lexicalTable.size() % 10000 == 0)
            std::cerr << ".";
    
        srcWord = it->first.first;
    }
    std::cerr << "Loaded " << m_lexicalTable.size() << " lexical pairs" << std::endl;
    std::cerr << std::endl;
}

void PhrasetableCreator::createRankHash() {
    
    // @TODO: compute rank in memory per source phrase 
    
    InputFileStream inFile(m_inPath);
     
    std::string line, prevSourcePhrase = "";
    size_t phr_num = 0;
    size_t line_num = 0;
    size_t numElement = NOT_FOUND;
    
    size_t step = 1ul << m_orderBits;
    
    std::priority_queue<std::pair<float, size_t> > rankQueue;
    std::vector<std::string> sourceTargetPhrases;
    
    while(std::getline(inFile, line)) {
        if(sourceTargetPhrases.size() == step) {
            m_rnkHash.AddRange(sourceTargetPhrases);
            sourceTargetPhrases.clear();
        }
   
        std::vector<std::string> tokens;
        TokenizeMultiCharSeparator(tokens, line, m_separator);
        
        if (numElement == NOT_FOUND) {
            // init numElement
            numElement = tokens.size();
            std::cerr << numElement << std::endl;
            assert(numElement >= 3);
            // extended style: source ||| target ||| scores ||| [alignment] ||| [counts]
        }
       
        if (tokens.size() != numElement) {
            std::stringstream strme;
            strme << "Syntax error at " << m_inPath << ":" << line_num;
            UserMessage::Add(strme.str());
            abort();
        }   
        
        std::string sourcePhraseString = tokens[0];
        
        bool isLHSEmpty = (sourcePhraseString.find_first_not_of(" \t", 0)
                           == std::string::npos);
        if (isLHSEmpty) {
            TRACE_ERR( m_inPath << ":" << line_num
                      << ": pt entry contains empty target, skipping\n");
            continue;
        }
       
        if(sourcePhraseString != prevSourcePhrase && prevSourcePhrase != "") {            
            ++phr_num;
            if(phr_num % 100000 == 0)
              std::cerr << ".";
            if(phr_num % 5000000 == 0)
              std::cerr << "[" << phr_num << "]" << std::endl;
        
            m_ranks.resize(line_num + 1);
            int r = 0;
            while(!rankQueue.empty()) {
                m_ranks[rankQueue.top().second] = r++;
                rankQueue.pop();
            }
        }
        
        prevSourcePhrase = sourcePhraseString;
        
        std::string sourceTargetPhrase = makeSourceTargetKey(tokens[0], tokens[1]);
        sourceTargetPhrases.push_back(sourceTargetPhrase);
        
        std::vector<float> scores = Tokenize<float>(tokens[2]);
        rankQueue.push(std::make_pair(scores[2], line_num));
        line_num++;
    }
    
    m_rnkHash.AddRange(sourceTargetPhrases);
    
#ifdef WITH_THREADS
    m_rnkHash.WaitAll();
#endif

    m_ranks.resize(line_num + 1);
    int r = 0;
    while(!rankQueue.empty()) {
        m_ranks[rankQueue.top().second] = r++;
        rankQueue.pop();
    }

    std::cerr << std::endl;
}

inline std::string PhrasetableCreator::makeSourceKey(std::string &source) {
    return source + m_separator;
}

inline std::string PhrasetableCreator::makeSourceTargetKey(std::string &source, std::string &target) {
    return source + m_separator + target + m_separator;
}

void PhrasetableCreator::encodeTargetPhrases() {
    
    InputFileStream inFile(m_inPath);

#ifdef WITH_THREADS
    boost::thread_group threads;
    for (int i = 0; i < m_threads; ++i) {
        EncodingTask* et = new EncodingTask(inFile, *this);    
        threads.create_thread(*et);
    }
    threads.join_all();
#else
    EncodingTask* et = new EncodingTask(inFile, *this);
    (*et)();
    delete et;
#endif
    flushEncodedQueue(true);
}


void PhrasetableCreator::compressTargetPhrases() {    
#ifdef WITH_THREADS
    boost::thread_group threads;
    for (int i = 0; i < m_threads; ++i) {
        CompressionTask* ct = new CompressionTask(m_encodedTargetPhrases, *this);    
        threads.create_thread(*ct);
    }
    threads.join_all();
#else
    CompressionTask* ct = new CompressionTask(m_encodedTargetPhrases, *this);
    (*ct)();
    delete ct;
#endif
    flushCompressedQueue(true);
}

void PhrasetableCreator::calcHuffmanCodes() {
    std::cerr << "Creating Huffman codes for " << m_symbolCounter.size()
        << " symbols" << std::endl;
         
    m_symbolTree = new SymbolTree(m_symbolCounter.begin(),
                                  m_symbolCounter.end());      
    {
        size_t sum = 0, sumall = 0;
        for(SymbolCounter::iterator it = m_symbolCounter.begin();
            it != m_symbolCounter.end(); it++) {
            sumall += it->second * m_symbolTree->encode(it->first).size();
            sum    += it->second;
        }
        float bits = float(sumall)/sum;
        std::cerr << bits << " bits per encoded symbol" << std::endl;
        
        float entropy = 0;
        for(SymbolCounter::iterator it = m_symbolCounter.begin();
            it != m_symbolCounter.end(); it++) {
            entropy += -float(it->second)/float(sum)
                       * log(float(it->second)/float(sum))/log(2);
        }
        std::cerr << "Entropy for encoded symbols: " << entropy << std::endl;
        
    }
    
    
    std::vector<ScoreTree*>::iterator treeIt = m_scoreTrees.begin();
    for(std::vector<ScoreCounter*>::iterator it = m_scoreCounters.begin();
        it != m_scoreCounters.end(); it++) {
        
        if(m_quantize)
            (*it)->quantize(m_quantize);
        
        std::cerr << "Creating Huffman codes for " << (*it)->size()
            << " scores" << std::endl;
        
        *treeIt = new ScoreTree((*it)->begin(), (*it)->end());
        treeIt++;
    }
    
    if(m_useAlignmentInfo) {
        std::cerr << "Creating Huffman codes for " << m_alignCounter.size()
            << " alignment points" << std::endl;
        m_alignTree = new AlignTree(m_alignCounter.begin(), m_alignCounter.end());
        {  
            size_t sum = 0, sumall = 0;
            for(AlignCounter::iterator it = m_alignCounter.begin();
                it != m_alignCounter.end(); it++) {
                sumall += it->second * m_alignTree->encode(it->first).size();
                sum    += it->second;
            }
            
            float bits = float(sumall)/sum;
            std::cerr << bits << " bits per encoded alignment point" << std::endl;
            
            float entropy = 0;
            for(AlignCounter::iterator it = m_alignCounter.begin();
                it != m_alignCounter.end(); it++) {
                entropy += -float(it->second)/float(sum)
                           * log(float(it->second)/float(sum))/log(2);
            }
            std::cerr << "Entropy for encoded alignment point: " << entropy << std::endl;
        }
    }
}


void PhrasetableCreator::addSourceSymbolId(std::string& symbol) {
  if(m_sourceSymbolsMap.count(symbol) == 0) {
    unsigned value = m_sourceSymbolsMap.size();
    m_sourceSymbolsMap[symbol] = value;
  }
}

void PhrasetableCreator::addTargetSymbolId(std::string& symbol) {
  if(m_targetSymbolsMap.count(symbol) == 0) {
    unsigned value = m_targetSymbolsMap.size();
    m_targetSymbolsMap[symbol] = value;
  }
}

unsigned PhrasetableCreator::getSourceSymbolId(std::string& symbol) {
    boost::unordered_map<std::string, unsigned>::iterator it
        = m_sourceSymbolsMap.find(symbol);
        
    if(it != m_sourceSymbolsMap.end())   
        return it->second;
    else
        return m_sourceSymbolsMap.size();
}

unsigned PhrasetableCreator::getTargetSymbolId(std::string& symbol) {
    boost::unordered_map<std::string, unsigned>::iterator it
        = m_targetSymbolsMap.find(symbol);
        
    if(it != m_targetSymbolsMap.end())   
        return it->second;
    else
        return m_targetSymbolsMap.size();
}

unsigned PhrasetableCreator::getOrAddTargetSymbolId(std::string& symbol) {
    // @TODO: check speed
#ifdef WITH_THREADS
    boost::mutex::scoped_lock lock(m_mutex);
#endif
    boost::unordered_map<std::string, unsigned>::iterator it
        = m_targetSymbolsMap.find(symbol);
        
    if(it != m_targetSymbolsMap.end())   
        return it->second;
    else {
        unsigned value = m_targetSymbolsMap.size();
        m_targetSymbolsMap[symbol] = value;
        return value;
    }
}

unsigned PhrasetableCreator::getRank(unsigned srcIdx, unsigned trgIdx) {
  size_t srcTrgIdx = m_lexicalTableIndex[srcIdx];
  while(srcTrgIdx < m_lexicalTable.size()
    && srcIdx == m_lexicalTable[srcTrgIdx].first
    && m_lexicalTable[srcTrgIdx].second != trgIdx)
    srcTrgIdx++;
  
  if(srcTrgIdx < m_lexicalTable.size()
     && m_lexicalTable[srcTrgIdx].second == trgIdx)  
    return srcTrgIdx - m_lexicalTableIndex[srcIdx];
  else
    return m_lexicalTable.size();
}

unsigned PhrasetableCreator::encodeREncSymbol1(unsigned trgIdx) {
  assert((~(1 << 31)) > trgIdx);
  return trgIdx;
}

unsigned PhrasetableCreator::encodeREncSymbol2(unsigned pos, unsigned rank) {
  unsigned symbol = rank;
  symbol |= 1 << 30;
  symbol |= pos << 24;
  return symbol;
}

unsigned PhrasetableCreator::encodeREncSymbol3(unsigned rank) {
  unsigned symbol = rank;
  symbol |= 2 << 30;
  return symbol;
}

unsigned PhrasetableCreator::encodePREncSymbol1(unsigned trgIdx) {
  assert((~(1 << 31)) > trgIdx);
  return trgIdx;
}

unsigned PhrasetableCreator::encodePREncSymbol2(int left, int right, unsigned rank) {
  // "left" and "right" must be smaller than 2^5
  // "rank" must be smaller than 2^19  
  left  = left  + 32;
  right = right + 32;
  
  assert(64 > left);
  assert(64 > right);
  assert(524288 > rank);
  
  unsigned symbol = 0;
  symbol |=    1  << 31;
  symbol |= left  << 25;
  symbol |= right << 19;
  symbol |= rank;
  return symbol;
}

void PhrasetableCreator::encodeTargetPhraseNone(std::vector<std::string>& t,
                                                std::ostream& os)
{
    std::stringstream encodedTargetPhrase;
    int j = 0;
    while(j < t.size()) {
        unsigned targetSymbolId = getOrAddTargetSymbolId(t[j]);
        
        m_symbolCounter.increase(targetSymbolId);
        os.write((char*)&targetSymbolId, sizeof(targetSymbolId));
        j++;
    }
    
    unsigned stopSymbolId = getTargetSymbolId(m_phraseStopSymbol);
    os.write((char*)&stopSymbolId, sizeof(stopSymbolId));
    m_symbolCounter.increase(stopSymbolId);
}

void PhrasetableCreator::encodeTargetPhraseREnc(std::vector<std::string>& s,
                                                std::vector<std::string>& t,
                                                std::set<AlignPoint>& a,
                                                std::ostream& os) {
  
    std::stringstream encodedTargetPhrase;
  
    std::vector<std::vector<size_t> > a2(t.size());
    for(std::set<AlignPoint>::iterator it = a.begin(); it != a.end(); it++)
        a2[it->second].push_back(it->first);

    for(size_t i = 0; i < t.size(); i++) {
        unsigned idxTarget = getOrAddTargetSymbolId(t[i]);
        unsigned encodedSymbol = -1;
        
        unsigned bestSrcPos = s.size();
        unsigned bestDiff = s.size();
        unsigned bestRank = m_lexicalTable.size();
        unsigned badRank = m_lexicalTable.size();
        
        for(std::vector<size_t>::iterator it = a2[i].begin(); it != a2[i].end(); it++) {
            unsigned idxSource = getSourceSymbolId(s[*it]);
            size_t r = getRank(idxSource, idxTarget);
            //size_t r = badRank;
            if(r != badRank) {
                if(r < bestRank) {
                    bestRank = r;
                    bestSrcPos = *it;
                    bestDiff = abs(*it-i);
                }
                else if(r == bestRank && abs(*it-i) < bestDiff) {
                    bestSrcPos = *it;
                    bestDiff = abs(*it-i);
                }
            }
        }
      
        if(bestRank != badRank && bestSrcPos < s.size()) {
            if(bestSrcPos == i)
                encodedSymbol = encodeREncSymbol3(bestRank);
            else
                encodedSymbol = encodeREncSymbol2(bestSrcPos, bestRank);          
            a.erase(AlignPoint(bestSrcPos, i));
        }
        else {
            encodedSymbol = encodeREncSymbol1(idxTarget);
        }
      
        os.write((char*)&encodedSymbol, sizeof(encodedSymbol));
        m_symbolCounter.increase(encodedSymbol);
    }
    
    // encode?
    unsigned stopSymbolId = getTargetSymbolId(m_phraseStopSymbol);
    unsigned encodedSymbol = encodeREncSymbol1(stopSymbolId);
    os.write((char*)&encodedSymbol, sizeof(encodedSymbol));
    m_symbolCounter.increase(encodedSymbol);    
}

void PhrasetableCreator::encodeTargetPhrasePREnc(std::vector<std::string>& s,
                                                 std::vector<std::string>& t,
                                                 std::set<AlignPoint>& a,
                                                 size_t ownRank,
                                                 std::ostream& os)
{
    std::vector<unsigned> encodedSymbols(t.size());
    std::vector<unsigned> encodedSymbolsLengths(t.size(), 0);
     
    ConsistantPhrases cp(s.size(), t.size(), a.begin(), a.end());
    while(cp.size()) {
        ConsistantPhrases::Phrase p = cp.pop();
        //std::cout << "\t(" << p.i << "," << p.m << "," << p.j << "," << p.n << ")" << std::endl;
        
        std::stringstream key1;
        key1 << s[p.i];
        for(int i = p.i+1; i < p.i+p.m; i++)
            key1 << " " << s[i];
            
        std::stringstream key2;
        key2 << t[p.j];
        for(int i = p.j+1; i < p.j+p.n; i++)
            key2 << " " << t[i];
        
        int rank = -1;
        std::string key1Str = key1.str(), key2Str = key2.str();
        size_t idx = m_rnkHash[makeSourceTargetKey(key1Str, key2Str)];
        if(idx != m_rnkHash.GetSize())
            rank = m_ranks[idx];
        
        //std::cout << "\t" << key.str() << " : " << rank << std::endl;     
        
        
        if(rank >= 0) {
            if(p.m != s.size() || rank < ownRank) {   
                std::stringstream encodedSymbol;
                encodedSymbols[p.j] = encodePREncSymbol2(p.i-p.j, s.size()-(p.i+p.m), rank);
                encodedSymbolsLengths[p.j] = p.n;
               
                std::set<AlignPoint> tAlignment;
                for(std::set<AlignPoint>::iterator it = a.begin();
                    it != a.end(); it++)
                    if(it->first < p.i || it->first >= p.i + p.m
                       || it->second < p.j || it->second >= p.j + p.n)
                        tAlignment.insert(*it);
                a = tAlignment;
                cp.removeOverlap(p);  
            }
        }
        
    }
    
    std::stringstream encodedTargetPhrase;
    
    int j = 0;
    while(j < t.size()) {
        if(encodedSymbolsLengths[j] > 0) {
            unsigned encodedSymbol = encodedSymbols[j];
            m_symbolCounter.increase(encodedSymbol);
            os.write((char*)&encodedSymbol, sizeof(encodedSymbol));
            j += encodedSymbolsLengths[j];
        }
        else {
            unsigned targetSymbolId = getOrAddTargetSymbolId(t[j]);
            unsigned encodedSymbol = encodePREncSymbol1(targetSymbolId);
            m_symbolCounter.increase(encodedSymbol);
            os.write((char*)&encodedSymbol, sizeof(encodedSymbol));
            j++;
        }
    }
    
    unsigned stopSymbolId = getTargetSymbolId(m_phraseStopSymbol);
    unsigned encodedSymbol = encodePREncSymbol1(stopSymbolId);
    os.write((char*)&encodedSymbol, sizeof(encodedSymbol));
    m_symbolCounter.increase(encodedSymbol);
}

void PhrasetableCreator::encodeScores(std::vector<float>& scores, std::ostream& os) {
    size_t c = 0;
    float score;
    
    while(c < scores.size()) {
        score = scores[c];
        score = FloorScore(TransformScore(score));
        os.write((char*)&score, sizeof(score));
        m_scoreCounters[m_multipleScoreTrees ? c : 0]->increase(score);
        c++;
    }
}

void PhrasetableCreator::encodeAlignment(std::set<AlignPoint>& alignment,
                                         std::ostream& os)
{
    for(std::set<AlignPoint>::iterator it = alignment.begin();
        it != alignment.end(); it++) {
        os.write((char*)&(*it), sizeof(AlignPoint));
        m_alignCounter.increase(*it);
    }
    AlignPoint stop(-1, -1);
    os.write((char*) &stop, sizeof(AlignPoint));
    m_alignCounter.increase(stop);
}

std::string PhrasetableCreator::encodeLine(std::vector<std::string>& tokens, size_t ownRank) {        
    std::string sourcePhraseStr = tokens[0];
    std::string targetPhraseStr = tokens[1];
    std::string scoresStr = tokens[2];
    std::string alignmentStr = tokens[3];
    
    std::vector<std::string> s = Tokenize(sourcePhraseStr);
    
    // is it worth it to put a lock around that?
    size_t phraseLength = s.size();
    if(m_maxPhraseLength < phraseLength)
        m_maxPhraseLength = phraseLength;
    
    std::vector<std::string> t = Tokenize(targetPhraseStr);
    std::vector<float> scores = Tokenize<float>(scoresStr);
    
    std::set<AlignPoint> a;
    if(m_coding != None || m_useAlignmentInfo) {
        std::vector<size_t> positions = Tokenize<size_t>(alignmentStr, " \t-");
        for(size_t i = 0; i < positions.size(); i += 2) {
          a.insert(AlignPoint(positions[i], positions[i+1]));
        }
    }
    
    std::stringstream encodedTargetPhrase;
    
    if(m_coding == PREnc) {
        encodeTargetPhrasePREnc(s, t, a, ownRank, encodedTargetPhrase);
    }
    else if(m_coding == REnc) {
        encodeTargetPhraseREnc(s, t, a, encodedTargetPhrase);        
    }
    else {
        encodeTargetPhraseNone(t, encodedTargetPhrase);      
    }
    
    encodeScores(scores, encodedTargetPhrase);
    
    if(m_useAlignmentInfo)
        encodeAlignment(a, encodedTargetPhrase);
    
    return encodedTargetPhrase.str();
}

std::string PhrasetableCreator::compressEncodedCollection(std::string encodedCollection) {
  
    enum EncodeState {
        ReadSymbol, ReadScore, ReadAlignment,
        EncodeSymbol, EncodeScore, EncodeAlignment };
    EncodeState state = ReadSymbol;
  
    unsigned phraseStopSymbolId;
    if(m_coding == REnc)
        phraseStopSymbolId = encodeREncSymbol1(getTargetSymbolId(m_phraseStopSymbol));
    else if(m_coding == PREnc)
        phraseStopSymbolId = encodePREncSymbol1(getTargetSymbolId(m_phraseStopSymbol));
    else
        phraseStopSymbolId = getTargetSymbolId(m_phraseStopSymbol);
    AlignPoint alignStopSymbol(-1, -1);
  
    std::stringstream encodedStream(encodedCollection);
    encodedStream.unsetf(std::ios::skipws);
    
    std::string output;
    BitStream<> bitstream(output);

    unsigned symbol;
    float score;
    size_t currScore = 0;
    AlignPoint alignPoint;
                                       
    while(encodedStream) {
        switch(state) {
            case ReadSymbol:
                encodedStream.read((char*) &symbol, sizeof(unsigned));
                state = EncodeSymbol;
                break;
            case ReadScore:
                if(currScore == m_numScoreComponent) {
                    currScore = 0;
                    if(m_useAlignmentInfo)
                      state = ReadAlignment;
                    else
                      state = ReadSymbol;
                }
                else {
                    encodedStream.read((char*) &score, sizeof(float));
                    currScore++;
                    state = EncodeScore;
                }
                break;
            case ReadAlignment:
                encodedStream.read((char*) &alignPoint, sizeof(AlignPoint));
                state = EncodeAlignment;
                break;
            
            case EncodeSymbol:
                state = (symbol == phraseStopSymbolId) ? ReadScore : ReadSymbol;
                bitstream.putCode(m_symbolTree->encode(symbol));
                break;
            case EncodeScore: {
                    state = ReadScore;
                    size_t idx = m_multipleScoreTrees ? currScore-1 : 0;
                    if(m_quantize)
                        score = m_scoreCounters[idx]->lowerBound(score);
                    bitstream.putCode(m_scoreTrees[idx]->encode(score));
                }
                break;
            case EncodeAlignment:
                state = (alignPoint == alignStopSymbol) ? ReadSymbol : ReadAlignment;
                bitstream.putCode(m_alignTree->encode(alignPoint));
                break;
        }
    }
    
    return output;
}

void PhrasetableCreator::addEncodedLine(PackedItem& pi) {
    m_queue.push(pi);
}

void PhrasetableCreator::flushEncodedQueue(bool force) {
    while(!m_queue.empty() && m_lastFlushedLine + 1 == m_queue.top().getLine()) {
        PackedItem pi = m_queue.top();
        m_queue.pop();
        m_lastFlushedLine++;
        
        if(m_lastFlushedSourcePhrase != pi.getSrc()) {
            if(m_lastCollection.size()) {
                std::stringstream targetPhraseCollection;
                for(std::vector<std::string>::iterator it =
                    m_lastCollection.begin(); it != m_lastCollection.end(); it++)
                    targetPhraseCollection << *it;
                
                m_lastSourceRange.push_back(makeSourceKey(m_lastFlushedSourcePhrase));    
                m_encodedTargetPhrases.push_back(targetPhraseCollection.str());
                
                m_lastFlushedSourceNum++;
                if(m_lastFlushedSourceNum % 100000 == 0)
                    std::cerr << ".";
                if(m_lastFlushedSourceNum % 5000000 == 0)
                    std::cerr << "[" << m_lastFlushedSourceNum << "]" << std::endl;
                
                m_lastCollection.clear();
            }
        }
        
        if(m_lastSourceRange.size() == (1ul << m_orderBits)) {
            m_srcHash.AddRange(m_lastSourceRange);
            m_srcHash.SaveLastRange();
            m_srcHash.DropLastRange();
            m_lastSourceRange.clear();
        }
        
        m_lastFlushedSourcePhrase = pi.getSrc();
        if(m_coding == PREnc) {
            if(m_lastCollection.size() <= pi.getRank())
                m_lastCollection.resize(pi.getRank() + 1);
            m_lastCollection[pi.getRank()] = pi.getTrg();
        }
        else {
            m_lastCollection.push_back(pi.getTrg());
        }
    }
    
    if(force) {
        if(m_lastCollection.size()) {
            std::stringstream targetPhraseCollection;
            for(std::vector<std::string>::iterator it =
                m_lastCollection.begin(); it != m_lastCollection.end(); it++)
                targetPhraseCollection << *it;
                
            m_encodedTargetPhrases.push_back(targetPhraseCollection.str());
            
            m_lastCollection.clear();
        }
        
        m_srcHash.AddRange(m_lastSourceRange);
        m_lastSourceRange.clear();
        
#ifdef WITH_THREADS
        m_srcHash.WaitAll();
#endif
    
        m_srcHash.SaveLastRange();
        m_srcHash.DropLastRange();
        m_srcHash.FinalizeSave();
        
        m_lastFlushedLine = -1;
        m_lastFlushedSourceNum = -1;
    }
}

void PhrasetableCreator::addCompressedCollection(PackedItem& pi) {
    m_queue.push(pi); 
}

void PhrasetableCreator::flushCompressedQueue(bool force) {
    if(force || m_queue.size() > 10000) {
        while(!m_queue.empty() && m_lastFlushedLine + 1 == m_queue.top().getLine()) {
            PackedItem pi = m_queue.top();
            m_queue.pop();
            m_lastFlushedLine++;
                
            m_compressedTargetPhrases.push_back(pi.getTrg());
            
            if((pi.getLine()+1) % 100000 == 0)
                std::cerr << ".";
            if((pi.getLine()+1) % 5000000 == 0)
                std::cerr << "[" << (pi.getLine()+1) << "]" << std::endl;
    
        }
    }
    
    if(force)
        m_lastFlushedLine = -1;    
}

//****************************************************************************//

size_t EncodingTask::m_lineNum = 0;
#ifdef WITH_THREADS
boost::mutex EncodingTask::m_mutex;
boost::mutex EncodingTask::m_fileMutex;
#endif

EncodingTask::EncodingTask(InputFileStream& inFile, PhrasetableCreator& creator)
  : m_inFile(inFile), m_creator(creator) {}
  
void EncodingTask::operator()() {
    size_t lineNum = 0;
    bool readline = false;
    
    std::vector<std::string> lines;
    size_t max_lines = 1000;
    lines.reserve(max_lines);
    
    {
#ifdef WITH_THREADS
        boost::mutex::scoped_lock lock(m_fileMutex);
#endif
        std::string line;
        while(std::getline(m_inFile, line) && lines.size() < max_lines)
            lines.push_back(line);
        lineNum = m_lineNum;
        m_lineNum += lines.size();
    }
    
    std::vector<PackedItem> result;
    result.reserve(max_lines);
    
    while(lines.size()) {
        for(size_t i = 0; i < lines.size(); i++) {
            std::vector<std::string> tokens;
            Moses::TokenizeMultiCharSeparator(tokens, lines[i], m_creator.m_separator);
            
            size_t ownRank = 0;
            if(m_creator.m_coding == PhrasetableCreator::PREnc)
                ownRank = m_creator.m_ranks[lineNum + i];
            
            std::string encodedLine = m_creator.encodeLine(tokens, ownRank);
            
            PackedItem packedItem(lineNum + i, tokens[0], encodedLine, ownRank);
            result.push_back(packedItem);
        }
        lines.clear();
        
        {
#ifdef WITH_THREADS
            boost::mutex::scoped_lock lock(m_mutex);
#endif
            for(int i = 0; i < result.size(); i++) 
                m_creator.addEncodedLine(result[i]);
            m_creator.flushEncodedQueue();  
        }
        
        result.clear();
        lines.reserve(max_lines);
        result.reserve(max_lines);
        
#ifdef WITH_THREADS
        boost::mutex::scoped_lock lock(m_fileMutex);
#endif
        std::string line;
        while(std::getline(m_inFile, line) && lines.size() < max_lines)
            lines.push_back(line);
        lineNum = m_lineNum;
        m_lineNum += lines.size();
    }
}

//****************************************************************************//

size_t CompressionTask::m_collectionNum = 0;
#ifdef WITH_THREADS
boost::mutex CompressionTask::m_mutex;
#endif

CompressionTask::CompressionTask(StringVector<unsigned char, unsigned long,
                              MmapAllocator>& encodedCollections,
                              PhrasetableCreator& creator)
  : m_encodedCollections(encodedCollections), m_creator(creator) {}
  
void CompressionTask::operator()() {
    size_t collectionNum;
    {
#ifdef WITH_THREADS
        boost::mutex::scoped_lock lock(m_mutex);
#endif
        collectionNum = m_collectionNum;
        m_collectionNum++;
    }
    
    while(collectionNum < m_encodedCollections.size()) {
        std::string collection = m_encodedCollections[collectionNum];
        std::string compressedCollection
            = m_creator.compressEncodedCollection(collection);
    
        std::string dummy;
        PackedItem packedItem(collectionNum, dummy, compressedCollection, 0);
    
#ifdef WITH_THREADS
        boost::mutex::scoped_lock lock(m_mutex);
#endif
        m_creator.addCompressedCollection(packedItem);
        m_creator.flushCompressedQueue();
        
        collectionNum = m_collectionNum;  
        m_collectionNum++;    
    }
}

//****************************************************************************//

PackedItem::PackedItem(long line, std::string sourcePhrase,
           std::string packedTargetPhrase, size_t rank)
  : m_line(line), m_sourcePhrase(sourcePhrase),
    m_packedTargetPhrase(packedTargetPhrase), m_rank(rank) {}

long PackedItem::getLine() const { return m_line; }

const std::string& PackedItem::getSrc() const { return m_sourcePhrase; }

const std::string& PackedItem::getTrg() const { return m_packedTargetPhrase; }

size_t PackedItem::getRank() const { return m_rank; }

bool operator<(const PackedItem &pi1, const PackedItem &pi2) {
    if(pi1.getLine() < pi2.getLine())
        return false;
    return true;
}

}
