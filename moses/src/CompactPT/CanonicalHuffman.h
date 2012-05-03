#ifndef moses_CanonicalHuffman_h
#define moses_CanonicalHuffman_h

#include <string>

#include <boost/dynamic_bitset.hpp>

#include "Hufftree.h"

namespace Moses {

template <class Container = std::string>
class BitStream {
  private:
    Container& m_data;
    bool m_reverse;
    
    typename Container::const_iterator m_iterator;
    typename Container::value_type m_currentValue;
    
    size_t m_bitPos;
    size_t m_valueBits;
    
    typename Container::value_type m_mask;
    
  public:
    BitStream(Container &data, bool reverse = false)
    : m_data(data), m_iterator(m_data.begin()), m_reverse(reverse),
      m_valueBits(sizeof(typename Container::value_type) * 8),
      m_mask(m_reverse ? 1 : 1 << m_valueBits-1), m_bitPos(0) { }
    
    size_t remainingBits() {
        if(m_data.size() * m_valueBits < m_bitPos)
            return 0;
        return m_data.size() * m_valueBits - m_bitPos;
    }
    
    bool getNext() {
        if(m_bitPos % m_valueBits == 0) {
            if(m_iterator != m_data.end()) {
                m_currentValue = *m_iterator++;
            }
        }
        else {
            m_currentValue = m_reverse ? m_currentValue >> 1 : m_currentValue << 1;
        } 
        
        m_bitPos++;
        return (m_currentValue & m_mask);
    }
    
    void reset() {
        m_iterator = m_data.begin();
        m_bitPos = 0;
    }
};

template<typename PosType, typename DataType> class Hufftree;

template <typename Data, typename Code = size_t>
class CanonicalHuffman {
  private:
    std::vector<Data> m_symbols;
    
    std::vector<Code> m_firstCodes;
    std::vector<size_t> m_lengthIndex;
    
    typedef boost::unordered_map<Data, boost::dynamic_bitset<> > EncodeMap;
    EncodeMap m_encodeMap;
    
    template <class HT>
    void init(HT& huffTree) {
      std::cerr << "Canonical" << std::endl;
        
      std::vector<size_t> numLength;
      for(typename HT::encodemap::iterator it = huffTree.m_encoding.begin();
        it != huffTree.m_encoding.end(); it++) {
        size_t length = it->second.size();
        if(numLength.size() <= length)
          numLength.resize(length+1, 0);
        numLength[length]++;
      }
      
      m_symbols.resize(huffTree.m_encoding.size());
      m_lengthIndex.resize(numLength.size());
      m_lengthIndex[0] = 0; 
      for(size_t l = 1; l < numLength.size(); l++) {
        m_lengthIndex[l] = m_lengthIndex[l-1] + numLength[l-1];
        std::cerr << l << " " << numLength[l] << " " << m_lengthIndex[l] << std::endl;
      }
      
      size_t maxLength = numLength.size() - 1;
      m_firstCodes.resize(maxLength + 1, 0);
      for(size_t l = maxLength - 1; l > 0; l--)
        m_firstCodes[l] = (m_firstCodes[l+1] + numLength[l+1])/2;
      
      for(size_t l = 1; l < m_firstCodes.size(); l++) {
        if(numLength[l]) {
          boost::dynamic_bitset<> x(l, m_firstCodes[l]);
          std::cerr << l << " " << x << std::endl;
        }
      }
      
      std::vector<size_t> nextCode = m_firstCodes;
      
      for(typename HT::encodemap::iterator it = huffTree.m_encoding.begin();
        it != huffTree.m_encoding.end(); it++) {
        
        Data data = it->first;
        size_t length = it->second.size();
        
        size_t pos = m_lengthIndex[length]
                     + (nextCode[length] - m_firstCodes[length]);
        m_symbols[pos] = data;
        
        nextCode[length] = nextCode[length] + 1;
      }
    }
    
  public:
    CanonicalHuffman(std::FILE* pFile, bool forEncoding = false) {
      load(pFile);
      
      if(forEncoding)
        createCodeMap();
    }
    
    template <class Iterator>
    CanonicalHuffman(Iterator begin, Iterator end, bool forEncoding = true) {
      Hufftree<int, Data> temp(begin, end);
      init(temp);

      if(forEncoding)
        createCodeMap();
    }
    
    template <class HT>
    CanonicalHuffman(HT& huffTree, bool forEncoding = true) {
      init(huffTree);

      if(forEncoding)
        createCodeMap();
    }
    
    void createCodeMap() {
      size_t maxLength = m_lengthIndex.size() - 1;
      for(size_t l = 1; l <= maxLength; l++) {
        Code code = m_firstCodes[l];
        
        size_t num = ((l+1 < m_lengthIndex.size()) ? m_lengthIndex[l+1]
                      : m_symbols.size()) - m_lengthIndex[l];
        
        for(size_t i = 0; i < num; i++) {
          Data data = m_symbols[m_lengthIndex[l] + i];  
          boost::dynamic_bitset<> bitCode(l, code);
          m_encodeMap[data] = bitCode;
          
          //if(l <= 10)
          //  std::cerr << l << " " << m_symbols[m_lengthIndex[l] + i] << " " << bitCode << std::endl;
          //
          code++;
        }
      }
    }
    
    boost::dynamic_bitset<>& encode(Data data) {
      return m_encodeMap[data];
    }
    
    template <class BitStream>
    Data nextSymbol(BitStream& bitStream) {
      if(bitStream.remainingBits()) {
        Code code = bitStream.getNext();
        size_t length = 1;
        
        while(code < m_firstCodes[length]) {
          code = 2 * code + bitStream.getNext();
          length++;
        }
        
        size_t symbolIndex = m_lengthIndex[length]
                             + (code - m_firstCodes[length]);
                               
        return m_symbols[symbolIndex];
      }   
      return Data();
    }
    
    size_t load(std::FILE* pFile) {
      size_t start = std::ftell(pFile);
      
      size_t size;
      std::fread(&size, sizeof(size_t), 1, pFile);
      m_symbols.resize(size);
      std::fread(&m_symbols[0], sizeof(Data), size, pFile);
      
      std::fread(&size, sizeof(size_t), 1, pFile);
      m_firstCodes.resize(size);
      std::fread(&m_firstCodes[0], sizeof(Code), size, pFile);
      
      std::fread(&size, sizeof(size_t), 1, pFile);
      m_lengthIndex.resize(size);
      std::fread(&m_lengthIndex[0], sizeof(size_t), size, pFile);
      
      return std::ftell(pFile) - start;
    }
    
    size_t save(std::FILE* pFile) {
      size_t start = std::ftell(pFile);
      
      size_t size = m_symbols.size();
      std::fwrite(&size, sizeof(size_t), 1, pFile);
      std::fwrite(&m_symbols[0], sizeof(Data), size, pFile);
      
      size = m_firstCodes.size();
      std::fwrite(&size, sizeof(size_t), 1, pFile);
      std::fwrite(&m_firstCodes[0], sizeof(Code), size, pFile);
      
      size = m_lengthIndex.size();
      std::fwrite(&size, sizeof(size_t), 1, pFile);
      std::fwrite(&m_lengthIndex[0], sizeof(size_t), size, pFile);
      
      return std::ftell(pFile) - start;
    }
    
};

}

#endif