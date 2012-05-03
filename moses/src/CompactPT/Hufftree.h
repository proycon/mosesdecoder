/* Copyright (c) 2011 the authors listed at the following URL, and/or
the authors of referenced articles or incorporated external code:
http://en.literateprograms.org/Huffman_coding_(C_Plus_Plus)?action=history&offset=20090129100015

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Retrieved from: http://en.literateprograms.org/Huffman_coding_(C_Plus_Plus)?oldid=16057
*/

#ifndef HUFFTREE_H__
#define HUFFTREE_H__

#include <vector>
#include <map>

#include <queue>
#include <deque>
#include <iostream>
#include <cstdio>
#include <cassert>

#include <boost/unordered_map.hpp>

#include "CanonicalHuffman.h"

namespace Moses {

template<typename PosType, typename DataType> class Hufftree {
  public:
    template<typename InputIterator>
    Hufftree(InputIterator begin, InputIterator end, bool = true);
    Hufftree(std::FILE*, bool = false);
    ~Hufftree() {}
  
    void Save(std::FILE*);
  
    std::vector<bool> encode(DataType const& value) {
      return m_encoding[value];
    }
  
    template<typename InputIterator>
    std::string encode(InputIterator begin, InputIterator end);
    
    template<typename InputIterator>
    std::string encodeWithLength(InputIterator begin, InputIterator end);
  
    template<typename InputIterator, typename OutputIterator>
    void decode(DataType length, InputIterator begin, InputIterator end, OutputIterator iter);
  
    template<typename InputIterator, typename OutputIterator>
    void decodeWithLength(InputIterator begin, InputIterator end, OutputIterator iter);
    
    inline PosType node(size_t i) {
      return m_nodes[i];
    }

    inline DataType data(size_t i) {
      return m_data[i];
    }
    
    inline size_t freq(size_t i) {
      return m_frequencies[i];
    }
    
    inline size_t size() {
      return m_data.size();
    }
  
  private:
    typedef boost::unordered_map<DataType, std::vector<bool> > encodemap;
    
    void fill(PosType index, std::vector<bool>& prefix);
    
    std::vector<PosType> m_nodes;
    std::vector<DataType> m_data;
    std::vector<size_t> m_frequencies;
    
    bool m_forEncoding;
    encodemap m_encoding;
    
    class Node;
    class NodeOrder;
    
    template <typename T1, typename T2> friend class CanonicalHuffman;
    //friend class CanonicalDecoder<DataType>;
};

template<typename PosType, typename DataType>
struct Hufftree<PosType, DataType>::Node {
  size_t frequency;
  Node* leftChild;
  
  union {
    Node* rightChild; // if leftChild != 0
    DataType* data;  // if leftChild == 0
  };

  Node(size_t f, DataType d):
    frequency(f),
    leftChild(0),
    data(new DataType(d))
  {}

  Node(Node* left, Node* right):
    frequency(left->frequency + right->frequency),
    leftChild(left),
    rightChild(right)
  {}

  ~Node() {
    if (leftChild) {
      delete leftChild;
      delete rightChild;
    }
    else {
      delete data;
    }
  }
};

template<typename PosType, typename DataType>
struct Hufftree<PosType, DataType>::NodeOrder {
  bool operator()(Node* a, Node* b) {
    if (b->frequency < a->frequency)
      return true;
    if (a->frequency < b->frequency)
      return false;

    if (!a->leftChild && b->leftChild)
      return true;
    if (a->leftChild && !b->leftChild)
      return false;

    if (a->leftChild && b->leftChild) {
      if ((*this)(a->leftChild, b->leftChild))
        return true;
      if ((*this)(b->leftChild, a->leftChild))
        return false;
      return (*this)(a->rightChild, b->rightChild);
    }

    return *(a->data) < *(b->data);
  }
};

template<typename PosType, typename DataType>
template<typename InputIterator>
Hufftree<PosType, DataType>::Hufftree(InputIterator begin, InputIterator end, bool forEncoding)
 : m_forEncoding(forEncoding){
  std::priority_queue<Node*, std::vector<Node*>, NodeOrder> pqueue;

  if(std::distance(begin, end) == 1) {
    std::vector<bool> zero(1, false);
    m_encoding[begin->first] = zero;
    m_nodes.resize(1, -1);
    m_data.push_back(begin->first);
    m_frequencies.push_back(begin->second);
    return;
  }

  while (begin != end) {
    Node* dataNode = new Node(begin->second, begin->first);
    pqueue.push(dataNode);
    ++begin;
  }

  Node* tree = 0;
  while (!pqueue.empty()) {
    Node* top = pqueue.top();
    pqueue.pop();
    if (pqueue.empty()) {
      tree = top;
    }
    else {
      Node* top2 = pqueue.top();
      pqueue.pop();
      pqueue.push(new Node(top, top2));
    }
  }

  std::deque<Node*> queue;
  queue.push_back(tree);
  
  PosType node_no = 0;
  PosType symb_no = 0;
  
  while (!queue.empty()) {
    Node* actual = queue.front();
    queue.pop_front();
    
    if (actual->leftChild) {
      queue.push_back(actual->leftChild);
      queue.push_back(actual->rightChild);
      
      if (node_no > 0)
        m_nodes.push_back(node_no);
        
      node_no += 2;
    }
    else {
      symb_no--;
      m_nodes.push_back(symb_no);
      m_data.push_back(*actual->data);
      m_frequencies.push_back(actual->frequency);
    }
  }

  delete tree;

  if(m_forEncoding) {
    std::vector<bool> bitvec;
    fill(0, bitvec);
  }
}

template<typename PosType, typename DataType>
Hufftree<PosType, DataType>::Hufftree(std::FILE* in, bool forEncoding)
  : m_forEncoding(forEncoding) {
  std::priority_queue<Node*, std::vector<Node*>, NodeOrder> pqueue;

  size_t length1;
  fread(&length1, sizeof(size_t), 1, in);
  m_nodes.resize(length1);
  fread(&(m_nodes[0]), sizeof(PosType), length1, in);
  
  size_t length2;
  fread(&length2, sizeof(size_t), 1, in);
  m_data.resize(length2);
  fread(&(m_data[0]), sizeof(DataType), length2, in);
  
  size_t length3;
  fread(&length3, sizeof(size_t), 1, in);
  m_frequencies.resize(length3);
  fread(&(m_frequencies[0]), sizeof(size_t), length3, in);
 
  if(m_forEncoding) {
    std::vector<bool> bitvec;
    fill(0, bitvec);
  }
}

template<typename PosType, typename DataType>
void Hufftree<PosType, DataType>::fill(PosType index, std::vector<bool>& prefix) {
  if(index >= 0) {
    prefix.push_back(0);
    fill(m_nodes[index], prefix);
    prefix.back() = 1;
    fill(m_nodes[index+1], prefix);
    prefix.pop_back();
  }
  else
    m_encoding[m_data[-index-1]] = prefix;
}


template<typename PosType, typename DataType>
void Hufftree<PosType, DataType>::Save(std::FILE* out) {
  size_t length1 = m_nodes.size();
  fwrite(&length1, sizeof(size_t), 1, out);
  fwrite(&(m_nodes[0]), sizeof(PosType), length1, out);
  
  size_t length2 = m_data.size();
  fwrite(&length2, sizeof(size_t), 1, out);
  fwrite(&(m_data[0]), sizeof(DataType), length2, out);
  
  size_t length3 = m_frequencies.size();
  fwrite(&length3, sizeof(size_t), 1, out);
  fwrite(&(m_frequencies[0]), sizeof(size_t), length3, out);
}

template<typename PosType, typename DataType>
template<typename InputIterator>
std::string Hufftree<PosType, DataType>::encode(InputIterator begin, InputIterator end) {
  assert(m_forEncoding);
  
  std::string result;
  
  char byte = 0;
  char mask = 1;
  unsigned int pos = 0;
  
  while(begin != end) {    
    std::vector<bool> code = m_encoding[*begin];
    begin++;
    
    for(size_t j = 0; j < code.size(); j++) {
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
  }

  if(pos % 8 != 0)
    result.push_back(byte);
  
  return result;
}

template<typename PosType, typename DataType>
template<typename InputIterator>
std::string Hufftree<PosType, DataType>::encodeWithLength(InputIterator begin, InputIterator end) {
  assert(m_forEncoding);

  std::string result;
  
  char byte = 0;
  char mask = 1;
  unsigned int pos = 0;
  
  size_t lengthNum = std::distance(begin, end);
  DataType length = lengthNum;
  std::vector<bool> code = m_encoding[length];
  
  for(size_t j = 0; j < code.size(); j++) {
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

  while(begin != end) {    
    std::vector<bool> code = m_encoding[*begin];
    begin++;
    
    for(size_t j = 0; j < code.size(); j++) {
      
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
  }
  
  if(pos % 8 != 0)
    result.push_back(byte);
  
  return result;
}

template<typename PosType, typename DataType>
template<typename InputIterator, typename OutputIterator>
void Hufftree<PosType, DataType>::decode(DataType length, InputIterator begin, InputIterator end, OutputIterator iter) {
  size_t count = 0;
  
  PosType node = 0;
  for(InputIterator it = begin; it != end; it++) {
    char byte = *it;
    char mask = 1;
    
    for(int i = 0; i < 8 && count < length; i++) {
      node = (byte & mask) ? m_nodes[node+1] : m_nodes[node];
      if(node < 0) {
        *iter++ = m_data[-node-1];
        count++;
        node = 0;
      }
      mask = mask << 1;
    }
  
    if(count >= length)
      return;
  }
}

template<typename PosType, typename DataType>
template<typename InputIterator, typename OutputIterator>
void Hufftree<PosType, DataType>::decodeWithLength(InputIterator begin, InputIterator end, OutputIterator iter) {
  bool has_length = 0;
  DataType length = 0;
  DataType count  = 0;
  
  PosType node = 0;
  for(InputIterator it = begin; it != end; it++) {
    char byte = *it;
    char mask = 1;
    
    int i = 0;
    for(; i < 8 && !has_length; i++) {
      node = (byte & mask) ? m_nodes[node+1] : m_nodes[node];
      if(node < 0) {
        length = m_data[-node-1];
        node = 0;
        has_length = true;
      }
      mask = mask << 1;
    }
  
    if(!has_length)
      continue;
    
    for(; i < 8; i++) {
      node = (byte & mask) ? m_nodes[node+1] : m_nodes[node];
      if(node < 0) {
        *iter++ = m_data[-node-1];
        node = 0;
        
        count++;
        if(count >= length)
          return;
      }
      mask = mask << 1;
    }
  }
}

}

#endif
