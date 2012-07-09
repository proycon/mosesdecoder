#ifndef moses_LexicalReorderingTableCreator_h
#define moses_LexicalReorderingTableCreator_h

#include "PhraseTableCreator.h"

namespace Moses {

class LexicalReorderingTableCreator {
  private:
    std::string m_inPath;
    std::string m_outPath;
    
    std::string m_separator;
    
    size_t m_threads;
    
    typedef Counter<float> ScoreCounter;
    typedef CanonicalHuffman<float> ScoreTree;  
  
    std::vector<ScoreCounter*> m_scoreCounters;
    std::vector<ScoreTree*> m_scoreTrees;
    
    BlockHashIndex m_hash;
    StringVector<unsigned char, unsigned long, MmapAllocator>  m_scores;
    
    void EncodeScores();
    void CalcHuffmanCodes();
    void CompressScores();
    void Save();
    
    std::string EncodeLine(std::vector<std::string>& tokens);
    void AddEncodedLine(PackedItem& pi);
    void FlushEncodedQueue(bool force = false);
    
    std::string CompressEncodedCollection(std::string encodedCollection);
    void AddCompressedCollection(PackedItem& pi);
    void FlushCompressedQueue(bool force = false);
    
  public:
    LexicalReorderingTableCreator(std::string inPath, std::string outPath);
    
  friend class EncodingTaskReordering;
  friend class CompressionTaskReordering;
};

class EncodingTaskReordering
{
  private:
#ifdef WITH_THREADS
    static boost::mutex m_mutex;
    static boost::mutex m_fileMutex;
#endif
    static size_t m_lineNum;
    static size_t m_sourcePhraseNum;
    static std::string m_lastSourcePhrase;
    
    InputFileStream& m_inFile;
    LexicalReorderingTableCreator& m_creator;
    
  public:
    EncodingTaskReordering(InputFileStream& inFile, LexicalReorderingTableCreator& creator);
    void operator()();
};

class CompressionTaskReordering
{
  private:
#ifdef WITH_THREADS
    static boost::mutex m_mutex;
#endif
    static size_t m_collectionNum;
    StringVector<unsigned char, unsigned long, MmapAllocator> &m_encodedScores;
    LexicalReorderingTableCreator &m_creator;
    
  public:
    CompressionTaskReordering(StringVector<unsigned char, unsigned long, MmapAllocator>&
                    m_encodedScores, LexicalReorderingTableCreator& creator);
    void operator()();
};

}

/*
 
 //std::vector<float> LexicalReorderingTableMemoryHashed::UnpackScores(std::string& scoreString) {
//  std::stringstream ss(scoreString);
//  std::vector<float> p;
//  
//  float score;
//  while(ss.read((char*) &score, sizeof(score))) {
//    p.push_back(score);
//  }
//  std::transform(p.begin(),p.end(),p.begin(),TransformScore);
//  std::transform(p.begin(),p.end(),p.begin(),FloorScore);
//  
//  return p;
//}

void  LexicalReorderingTableMemoryCompact::LoadText(const std::string& filePath)
{
  std::vector<char*> tempScores;
  std::map<float, size_t> frequencies;
  
  StringVector<unsigned char, unsigned long, MmapAllocator> phrases;
  
  std::cerr << "Reading in reordering table" << std::endl;
  
  std::string fileName = filePath;
  if(!FileExists(fileName) && FileExists(fileName+".gz")) {
    fileName += ".gz";
  }
  InputFileStream file(fileName);
  std::string line(""), key("");
  int numScores = -1;
  size_t line_num = 0;
  while(!getline(file, line).eof()) {
    ++line_num;
    if(line_num % 100000 == 0)
      std::cerr << ".";
    if(line_num % 5000000 == 0)
      std::cerr << "[" << line_num << "]" << std::endl;
      
    std::vector<std::string> tokens = TokenizeMultiCharSeparator(line, "|||");
    int t = 0 ;
    std::string f(""),e(""),c("");

    if(!m_FactorsF.empty()) {
      //there should be something for f
      f = auxClearString(tokens.at(t));
      ++t;
    }
    if(!m_FactorsE.empty()) {
      //there should be something for e
      e = auxClearString(tokens.at(t));
      ++t;
    }
    if(!m_FactorsC.empty()) {
      //there should be something for c
      c = auxClearString(tokens.at(t));
      ++t;
    }
    //last token are the probs
    std::vector<float> p = Scan<float>(Tokenize(tokens.at(t)));
    //sanity check: all lines must have equall number of probs
    if(-1 == numScores) {
      numScores = (int)p.size(); //set in first line
    }
    if((int)p.size() != numScores) {
      TRACE_ERR( "found inconsistent number of probabilities... found " << p.size() << " expected " << numScores << std::endl);
      exit(0);
    }
        
    phrases.push_back(MakeKey(f,e,c));
    std::transform(p.begin(), p.end(), p.begin(), TransformScore);
    std::transform(p.begin(), p.end(), p.begin(), FloorScore);
    
    for(std::vector<float>::iterator it = p.begin(); it != p.end(); it++)
      frequencies[*it]++;
    
    char* cstring = new char[numScores * sizeof(float)];
    std::memcpy(cstring, &p[0], numScores * sizeof(float));
    tempScores.push_back(cstring);
  }
  std::cerr << std::endl;
  
  //m_hash.Create(phrases);
  //
  //{
  //  StringVector<unsigned char, unsigned long, MmapAllocator> tPhrases;
  //  phrases.swap(tPhrases);
  //}

  // TODO!
  //std::cerr << "Creating Huffman compression tree for " << frequencies.size() << " symbols" << std::endl;
  //m_tree = new Hufftree<int, float>(frequencies.begin(), frequencies.end());
  //
  //double freq_sum = 0, len_sum = 0;
  //for(std::map<float, size_t>::iterator it = frequencies.begin(); it != frequencies.end(); it++) {
  //  len_sum  += it->second * m_tree->encode(it->first).size();
  //  freq_sum += it->second;
  //}
  //std::cerr << "Average no. of bits per score: " << (len_sum/freq_sum) << std::endl;
  //std::cerr << "Compressing target phrases" << std::endl;
  //for(size_t i = 0; i < m_hash.GetSize(); i++) {
  //  if((i+1) % 100000 == 0)
  //    std::cerr << ".";
  //  if((i+1) % 5000000 == 0)
  //    std::cerr << "[" << i+1 << "]" << std::endl;
  //    
  //  std::vector<float> p(numScores, 0);
  //  char* cstring = tempScores[i];
  //  std::memcpy(&p[0], cstring, numScores * sizeof(float));
  //  delete[] cstring;
  //  
  //  std::string compressedScores = m_tree->encode(p.begin(), p.end());
  //  m_scores.push_back(compressedScores);
  //}
  std::cerr << std::endl;
}

 */

#endif