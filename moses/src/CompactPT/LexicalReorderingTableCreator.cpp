#include "LexicalReorderingTableCreator.h"

namespace Moses {

LexicalReorderingTableCreator::LexicalReorderingTableCreator(
    std::string inPath, std::string outPath)
    : m_inPath(inPath), m_outPath(outPath), m_separator(" ||| "),
      m_threads(8), m_hash(10, 16)
{
  EncodeScores();
  CompressScores();
  Save();
}

void LexicalReorderingTableCreator::EncodeScores()
{
  InputFileStream inFile(m_inPath);

#ifdef WITH_THREADS
  boost::thread_group threads;
  for (size_t i = 0; i < m_threads; ++i)
  {
    EncodingTaskReordering* et = new EncodingTaskReordering(inFile, *this);    
    threads.create_thread(*et);
  }
  threads.join_all();
#else
  EncodingTaskReordering* et = new EncodingTaskReordering(inFile, *this);
  (*et)();
  delete et;
#endif
  FlushEncodedQueue(true);
}

void LexicalReorderingTableCreator::CompressScores()
{
#ifdef WITH_THREADS
  boost::thread_group threads;
  for (size_t i = 0; i < m_threads; ++i) {
    CompressionTaskReordering* ct = new CompressionTaskReordering(m_scores, *this);    
    threads.create_thread(*ct);
  }
  threads.join_all();
#else
  CompressionTaskReordering* ct = new CompressionTaskReordering(m_scores, *this);
  (*ct)();
  delete ct;
#endif
  FlushCompressedQueue(true);
}

void LexicalReorderingTableCreator::Save()
{
  
}

std::string LexicalReorderingTableCreator::EncodeLine(std::vector<std::string>& tokens)
{
  
}

void LexicalReorderingTableCreator::AddEncodedLine(PackedItem& pi)
{
  
}

void LexicalReorderingTableCreator::FlushEncodedQueue(bool force) {
  
}

std::string LexicalReorderingTableCreator::CompressEncodedCollection(std::string encodedCollection) {
  
}

void LexicalReorderingTableCreator::AddCompressedCollection(PackedItem& pi) {
  
}

void LexicalReorderingTableCreator::FlushCompressedQueue(bool force)
{  

}

//****************************************************************************//

size_t EncodingTaskReordering::m_lineNum = 0;
#ifdef WITH_THREADS
boost::mutex EncodingTaskReordering::m_mutex;
boost::mutex EncodingTaskReordering::m_fileMutex;
#endif

EncodingTaskReordering::EncodingTaskReordering(InputFileStream& inFile, LexicalReorderingTableCreator& creator)
  : m_inFile(inFile), m_creator(creator) {}
  
void EncodingTaskReordering::operator()()
{
  size_t lineNum = 0;
  
  std::vector<std::string> lines;
  size_t max_lines = 1000;
  lines.reserve(max_lines);
  
  {
#ifdef WITH_THREADS
    boost::mutex::scoped_lock lock(m_fileMutex);
#endif
    std::string line;
    while(lines.size() < max_lines && std::getline(m_inFile, line))
        lines.push_back(line);
    lineNum = m_lineNum;
    m_lineNum += lines.size();
  }
  
  std::vector<PackedItem> result;
  result.reserve(max_lines);
  
  while(lines.size())
  {
    for(size_t i = 0; i < lines.size(); i++)
    {
      std::vector<std::string> tokens;
      Moses::TokenizeMultiCharSeparator(tokens, lines[i], m_creator.m_separator);
      
      std::string encodedLine = m_creator.EncodeLine(tokens);
      
      PackedItem packedItem(lineNum + i, tokens[0], encodedLine, i);
      result.push_back(packedItem);
    }
    lines.clear();
    
    {
#ifdef WITH_THREADS
      boost::mutex::scoped_lock lock(m_mutex);
#endif
      for(size_t i = 0; i < result.size(); i++) 
        m_creator.AddEncodedLine(result[i]);
      m_creator.FlushEncodedQueue();  
    }
    
    result.clear();
    lines.reserve(max_lines);
    result.reserve(max_lines);
    
#ifdef WITH_THREADS
    boost::mutex::scoped_lock lock(m_fileMutex);
#endif
    std::string line;
    while(lines.size() < max_lines && std::getline(m_inFile, line))
      lines.push_back(line);
    lineNum = m_lineNum;
    m_lineNum += lines.size();
  }
}

//****************************************************************************//

size_t CompressionTaskReordering::m_collectionNum = 0;
#ifdef WITH_THREADS
boost::mutex CompressionTaskReordering::m_mutex;
#endif

CompressionTaskReordering::CompressionTaskReordering(StringVector<unsigned char, unsigned long,
                              MmapAllocator>& encodedScores,
                              LexicalReorderingTableCreator& creator)
  : m_encodedScores(encodedScores), m_creator(creator) {}
  
void CompressionTaskReordering::operator()()
{
  size_t collectionNum;
  {
#ifdef WITH_THREADS
    boost::mutex::scoped_lock lock(m_mutex);
#endif
    collectionNum = m_collectionNum;
    m_collectionNum++;
  }
  
  while(collectionNum < m_encodedScores.size())
  {
    std::string collection = m_encodedScores[collectionNum];
    std::string compressedCollection
        = m_creator.CompressEncodedCollection(collection);

    std::string dummy;
    PackedItem packedItem(collectionNum, dummy, compressedCollection, 0);

#ifdef WITH_THREADS
    boost::mutex::scoped_lock lock(m_mutex);
#endif
    m_creator.AddCompressedCollection(packedItem);
    m_creator.FlushCompressedQueue();
    
    collectionNum = m_collectionNum;  
    m_collectionNum++;    
  }
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
}
