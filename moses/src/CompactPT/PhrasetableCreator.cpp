#include <cstdio>

#include "PhrasetableCreator.h"

namespace Moses {
    
size_t maxPhrases = 0;
    
std::string PhrasetableCreator::m_phraseStopSymbol = "__SPECIAL_STOP_SYMBOL__";
    
PhrasetableCreator::PhrasetableCreator(std::string inPath, std::string outPath,
                                       Coding coding, size_t orderBits,
                                       size_t fingerPrintBits)
  : m_inPath(inPath), m_outPath(outPath),
    m_outFile(std::fopen(m_outPath.c_str(), "w")), m_numScoreComponent(5),
    m_coding(coding), m_orderBits(orderBits), m_fingerPrintBits(fingerPrintBits), 
#ifdef WITH_THREADS
    m_threads(6),
    m_srcHash(m_orderBits, m_fingerPrintBits, m_threads/2),
    m_rnkHash(m_orderBits, m_fingerPrintBits, m_threads/2),
#else
    m_srcHash(m_orderBits, m_fingerPrintBits),
    m_rnkHash(m_orderBits, m_fingerPrintBits),
#endif
    m_origSymbolCount(0), m_origScoreCount(0), m_origAlignCount(0),
    m_lastFlushedLine(-1), m_lastFlushedSourceNum(0),
    m_lastFlushedSourcePhrase("")
{
    addTargetSymbolId(m_phraseStopSymbol);
    
    createHashes();
    
    std::cerr << "srcHash: " << m_srcHash.GetSize() << std::endl;
    std::cerr << "TargetSymbols: " << m_targetSymbolsMap.size() << std::endl;
    
    if(m_coding == REnc) {
        std::vector<std::string> temp1;
        temp1.resize(m_sourceSymbolsMap.size());
        for(boost::unordered_map<std::string, unsigned>::iterator it
            = m_sourceSymbolsMap.begin(); it != m_sourceSymbolsMap.end(); it++)
            temp1[it->second] = it->first;
            
        std::sort(temp1.begin(), temp1.end());
        
        for(size_t i = 0; i < temp1.size(); i++)
            m_sourceSymbolsMap[temp1[i]] = i;
        
        loadLexicalTable("/home/marcinj/Poleng/Wipo/coppa/model/lex.f2e");
    }
            
    encodeTargetPhrases();
    
    std::cerr << "EncodedPhrases: " << m_encodedTargetPhrases.size() << std::endl;
        
    std::cerr << "SymbolCounter: " << m_symbolCounter.size() << std::endl;
    std::cerr << "ScoreCounter: " << m_scoreCounter.size() << std::endl;
    std::cerr << "AlignCounter: " << m_alignCounter.size() << std::endl;
    
    calcHuffmanCodes();
    
    compressTargetPhrases();
    
    std::cerr << "CompressedPhrases: " << m_compressedTargetPhrases.size() << std::endl;

    std::fwrite(&m_coding, sizeof(m_coding), 1, m_outFile);
    if(m_coding == REnc) {
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
        
        size_t size = m_lexicalTableIndex.size();
        std::fwrite(&size, sizeof(size_t), 1, m_outFile);
        std::fwrite(&m_lexicalTableIndex[0], sizeof(size_t), size, m_outFile);
        
        size = m_lexicalTable.size();
        std::fwrite(&size, sizeof(size_t), 1, m_outFile);
        std::fwrite(&m_lexicalTable[0], sizeof(SrcTrg), size, m_outFile);
    }
    
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
    
    m_symbolTree->save(m_outFile);
    m_scoreTree->save(m_outFile);
    m_alignTree->save(m_outFile);
    size_t phraseSize = m_compressedTargetPhrases.save(m_outFile);
    std::fclose(m_outFile);
}
    
void PhrasetableCreator::loadLexicalTable(std::string filePath) {
    std::vector<SrcTrgProb> t_lexTable;
      
    std::cerr << "Reading in lexical table from " << filePath << std::endl;
    std::ifstream lexIn(filePath.c_str(), std::ifstream::in);
    std::string src, trg;
    float prob;
    
    while(lexIn >> trg >> src >> prob) {
        if(t_lexTable.size() % 10000 == 0)
            std::cerr << ".";
        t_lexTable.push_back(SrcTrgProb(SrcTrgString(src, trg), prob));
    }
    std::cerr << std::endl;
    
    std::cerr << "Read in " << t_lexTable.size() << " lexical pairs" << std::endl;
    std::sort(t_lexTable.begin(), t_lexTable.end(), SrcTrgProbSorter());
    std::cerr << "Sorted" << std::endl;
    
    std::string srcWord = "";
    size_t srcIdx = 0;
    for(typename std::vector<SrcTrgProb>::iterator it = t_lexTable.begin();
        it != t_lexTable.end(); it++) {
      
        if(it->first.first != srcWord) {
            srcIdx = getSourceSymbolId(it->first.first);
            if(srcIdx < m_sourceSymbolsMap.size()) {
                if(srcIdx >= m_lexicalTableIndex.size())
                    m_lexicalTableIndex.resize(srcIdx + 1);
                m_lexicalTableIndex[srcIdx] = m_lexicalTable.size();
            }
        }
        
        size_t trgIdx = getTargetSymbolId(it->first.second);        
        if(srcIdx < m_sourceSymbolsMap.size() && trgIdx < m_targetSymbolsMap.size()) {
            m_lexicalTable.push_back(SrcTrg(srcIdx, trgIdx));
            
            if(m_lexicalTable.size() % 10000 == 0)
                std::cerr << ".";
        }
        
        srcWord = it->first.first;
    }
    std::cerr << "Kept " << m_lexicalTable.size() << " lexical pairs" << std::endl;
    std::cerr << std::endl;
}
    
void PhrasetableCreator::createHashes() {
    
    InputFileStream inFile(m_inPath);
     
    std::string line, prevSourcePhrase = "";
    size_t phr_num = 0;
    size_t line_num = 0;
    size_t numElement = NOT_FOUND;
    
    // poprawić
    size_t step = 1ul << m_orderBits;
    
    std::vector<std::string> sourcePhrases;
    std::vector<std::string> sourceTargetPhrases;

    m_srcHash.BeginSave(m_outFile);

    while(std::getline(inFile, line)) {
        
        if(sourcePhrases.size() == step) {
            m_srcHash.AddRange(sourcePhrases);
            m_srcHash.SaveLastRange();
            m_srcHash.DropLastRange();
            sourcePhrases.clear();
        }
        
        if(m_coding == PREnc) {
            if(sourceTargetPhrases.size() == step) {
                m_rnkHash.AddRange(sourceTargetPhrases);
                sourceTargetPhrases.clear();
            }
        }
        
        std::vector<std::string> tokens = Tokenize( line , "\t" );
        
        if (numElement == NOT_FOUND) {
            // init numElement
            numElement = tokens.size();
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
            sourcePhrases.push_back(Trim(prevSourcePhrase));
            
            if(maxPhrases && phr_num > maxPhrases)
                break;
            
            ++phr_num;
            if(phr_num % 100000 == 0)
              std::cerr << ".";
            if(phr_num % 5000000 == 0)
              std::cerr << "[" << phr_num << "]" << std::endl;
        }
        
        prevSourcePhrase = sourcePhraseString;
        
        std::vector<std::string> sourceSymbols = Tokenize(tokens[0]);
        for(std::vector<std::string>::iterator it = sourceSymbols.begin();
            it != sourceSymbols.end(); it++)
            addSourceSymbolId(*it);
       
        std::vector<std::string> targetSymbols = Tokenize(tokens[1]); 
        for(std::vector<std::string>::iterator it = targetSymbols.begin();
            it != targetSymbols.end(); it++)
            addTargetSymbolId(*it);
       
        if(m_coding == PREnc) {
            std::string sourceTargetPhrase = Trim(tokens[0]) + "\t" + Trim(tokens[1]);
            sourceTargetPhrases.push_back(sourceTargetPhrase);
            m_ranks.push_back(Scan<unsigned>(tokens[5]));
        }
    }
    
    m_srcHash.AddRange(sourcePhrases);
    
    if(m_coding == PREnc)
        m_rnkHash.AddRange(sourceTargetPhrases);
    
#ifdef WITH_THREADS
    m_srcHash.WaitAll();
    
    if(m_coding == PREnc)
        m_rnkHash.WaitAll();
#endif
    
    m_srcHash.SaveLastRange();
    m_srcHash.DropLastRange();
    m_srcHash.FinalizeSave();
    
    std::cerr << std::endl;
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
        std::cerr << (float(sumall)/m_origSymbolCount)
            << " bits per original symbol" << std::endl;
        
        float entropy = 0;
        for(SymbolCounter::iterator it = m_symbolCounter.begin();
            it != m_symbolCounter.end(); it++) {
            entropy += -float(it->second)/float(sum)
                       * log(float(it->second)/float(sum))/log(2);
        }
        std::cerr << "Entropy for encoded symbols: " << entropy << std::endl;
        std::cerr << "Entropy for original symbols: "
            <<  (entropy*sum/m_origSymbolCount) << std::endl;
        
    }
    
    std::cerr << "Creating Huffman codes for " << m_scoreCounter.size()
        << " scores" << std::endl;
    m_scoreTree = new ScoreTree(m_scoreCounter.begin(), m_scoreCounter.end());
    {  
        size_t sum = 0, sumall = 0;
        for(ScoreCounter::iterator it = m_scoreCounter.begin();
            it != m_scoreCounter.end(); it++) {
            sumall += it->second * m_scoreTree->encode(it->first).size();
            sum    += it->second;
        }
        std::cerr << double(sumall)/sum << " bits per score" << std::endl;
    }

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
        std::cerr << (bits*sum/m_origAlignCount)
            << " bits per original alignment point" << std::endl;
        
        float entropy = 0;
        for(AlignCounter::iterator it = m_alignCounter.begin();
            it != m_alignCounter.end(); it++) {
            entropy += -float(it->second)/float(sum)
                       * log(float(it->second)/float(sum))/log(2);
        }
        std::cerr << "Entropy for encoded alignment point: " << entropy << std::endl;
        std::cerr << "Entropy for original alignment point: "
            <<  (entropy*sum/m_origAlignCount) << std::endl;
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

void PhrasetableCreator::encodeTargetPhraseNone(std::vector<std::string>& s,
                                                std::vector<std::string>& t,
                                                std::set<AlignPoint>& a,
                                                std::ostream& os)
{
    std::stringstream encodedTargetPhrase;
    int j = 0;
    while(j < t.size()) {
        unsigned targetSymbolId = getTargetSymbolId(t[j]);
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
        unsigned idxTarget = getTargetSymbolId(t[i]);
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
        
        std::stringstream key;
        
        key << s[p.i];
        for(int i = p.i+1; i < p.i+p.m; i++)
            key << " " << s[i];
            
        key << "\t";
        
        key << t[p.j];
        for(int i = p.j+1; i < p.j+p.n; i++)
            key << " " << t[i];
        
        int rank = -1;
        size_t idx = m_rnkHash[key.str()];
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
            unsigned targetSymbolId = getTargetSymbolId(t[j]);
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
        m_scoreCounter.increase(score);
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

std::string PhrasetableCreator::encodeLine(std::vector<std::string>& tokens) {        
    std::string sourcePhraseStr = tokens[0];
    std::string targetPhraseStr = tokens[1];
    std::string scoresStr = tokens[2];
    std::string alignmentStr = tokens[3];
    
    size_t ownRank = Scan<size_t>(tokens[5]); 
    
    std::vector<std::string> s = Tokenize(sourcePhraseStr);
    std::vector<std::string> t = Tokenize(targetPhraseStr);
    
    std::vector<float> scores = Tokenize<float>(scoresStr);
    
    std::set<AlignPoint> a;
    std::vector<size_t> positions = Tokenize<size_t>(alignmentStr, " \t-");
    for(size_t i = 0; i < positions.size(); i += 2) {
      a.insert(AlignPoint(positions[i], positions[i+1]));
    }
    
    {
#ifdef WITH_THREADS
        boost::mutex::scoped_lock lock(m_mutex);
#endif
        m_origSymbolCount += t.size();
        m_origScoreCount += scores.size();
        m_origAlignCount += a.size();
    }
    
    std::stringstream encodedTargetPhrase;
    
    if(m_coding == PREnc) {
        encodeTargetPhrasePREnc(s, t, a, ownRank, encodedTargetPhrase);
    }
    else if(m_coding == REnc) {
        encodeTargetPhraseREnc(s, t, a, encodedTargetPhrase);        
    }
    else {
        encodeTargetPhraseNone(s, t, a, encodedTargetPhrase);      
    }
    encodeScores(scores, encodedTargetPhrase);
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
    std::string result;
  
    char byte = 0;
    char mask = 1;
    unsigned int pos = 0;
  
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
                    //if(m_containsAlignmentInfo)
                      state = ReadAlignment;
                    //else
                    //  state = ReadSymbol;
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
            case EncodeScore:
            case EncodeAlignment:
                boost::dynamic_bitset<> code;
                if(state == EncodeSymbol) {
                    code = m_symbolTree->encode(symbol);
                    if(symbol == phraseStopSymbolId)
                        state = ReadScore;
                    else
                        state = ReadSymbol;
                }
                else if(state == EncodeScore) {
                    //size_t idx = (currScore-1) % scoreTypes;
                    //float closestScore = m_scoreCounts[idx].lowerBound(score);
                    //code = m_scoreTrees[idx]->encode(closestScore);
                    code = m_scoreTree->encode(score);
                    state = ReadScore;
                }
                else if(state == EncodeAlignment) {
                    code = m_alignTree->encode(alignPoint);
                    if(alignPoint == alignStopSymbol)
                        state = ReadSymbol;
                    else
                        state = ReadAlignment;
                }
                
                for(int j = code.size()-1; j >= 0; j--) {
                    if(code[j])
                        byte |= mask;
                    mask = mask << 1;
                    pos++;
                    
                    if(pos % 8 == 0) {
                        result.push_back(byte);
                        mask = 1;
                        byte = 0;
                    }
                }
                break;
        }
    }
  
    // Add last byte with remaining waste bits
    if(pos % 8 != 0)
        result.push_back(byte);
    
    return result;
}

void PhrasetableCreator::addEncodedLine(PackedItem& pi) {
    m_queue.push(pi);
    flushEncodedQueue();  
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
                    
                m_encodedTargetPhrases.push_back(targetPhraseCollection.str());
                
                m_lastFlushedSourceNum++;
                if(m_lastFlushedSourceNum % 100000 == 0)
                    std::cerr << ".";
                if(m_lastFlushedSourceNum % 5000000 == 0)
                    std::cerr << "[" << m_lastFlushedSourceNum << "]" << std::endl;
                
                m_lastCollection.clear();
            }
        }
        
        m_lastFlushedSourcePhrase = pi.getSrc();
        if(m_lastCollection.size() <= pi.getRank())
            m_lastCollection.resize(pi.getRank() + 1);
        m_lastCollection[pi.getRank()] = pi.getTrg();
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
        m_lastFlushedLine = -1;
        m_lastFlushedSourceNum = -1;
    }
}

void PhrasetableCreator::addCompressedCollection(PackedItem& pi) {
    m_queue.push(pi);
    flushCompressedQueue();  
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
#endif

EncodingTask::EncodingTask(InputFileStream& inFile, PhrasetableCreator& creator)
  : m_inFile(inFile), m_creator(creator) {}
  
void EncodingTask::operator()() {
    std::string line;
    size_t lineNum = 0;
    bool readline = false;
    
    {
#ifdef WITH_THREADS
        boost::mutex::scoped_lock lock(m_mutex);
#endif
        if(std::getline(m_inFile, line))
            readline = true;
        lineNum = m_lineNum;
        m_lineNum++;
    }
    
    while(readline) {
        std::vector<std::string> tokens = Moses::Tokenize(line , "\t");
        std::string encodedLine = m_creator.encodeLine(tokens);
    
        size_t ownRank = Scan<size_t>(tokens[5]);
        PackedItem packedItem(lineNum, tokens[0], encodedLine, ownRank);
        
        readline = false;
#ifdef WITH_THREADS
        boost::mutex::scoped_lock lock(m_mutex);
#endif
        m_creator.addEncodedLine(packedItem);
        
        //if(maxPhrases && packedItem.getSrc() > maxPhrases)
        //    return;
        
        if(std::getline(m_inFile, line))
            readline = true;
        
        lineNum = m_lineNum;  
        m_lineNum++;
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
        
        //if(maxPhrases && packedItem.getSrc() > maxPhrases)
        //    return;
        
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
