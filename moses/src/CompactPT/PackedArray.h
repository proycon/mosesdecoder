#ifndef PACKEDARRAY_H__
#define PACKEDARRAY_H__

#include <vector>
#include <cmath>
#include <cstring>
#include <cstdio>

namespace Moses {

template <typename T = size_t, typename D = unsigned char>
class PackedArray {
  protected:
    static size_t m_bits;
    static size_t m_dataBits;
    
    size_t m_size;
    size_t m_storageSize;
    D* m_storage;
    
  public:
    PackedArray() {
        m_size = 0;
        m_storageSize = 0;
        m_storage = new D[0];
    }
    
    PackedArray(size_t bits, size_t size) {
        m_bits = bits;
        m_size = size;
        
        m_storageSize = ceil(float(m_bits * size) / float(m_dataBits));
        m_storage = new D[m_storageSize];
    }
    
    PackedArray(const PackedArray<T, D> &c) {
        m_size = c.m_size;
        
        m_storageSize = c.m_storageSize;
        m_storage = new D[m_storageSize];
        
        std::memcpy(m_storage, c.m_storage, m_storageSize * sizeof(D));
    }
    
    ~PackedArray() {
        delete [] m_storage;
        m_size = 0;
        m_storageSize = 0;
        m_storage = 0;
    }
    
    T get(size_t i) const {
        T out = 0;
        
        size_t bitstart = (i * m_bits);
        size_t bitpos = bitstart;
        
        size_t zero = ((1ul << (m_bits)) - 1);
        
        while(bitpos - bitstart < m_bits) {
            size_t pos = bitpos / m_dataBits;
            size_t off = bitpos % m_dataBits;
               
            out |= (T(m_storage[pos]) << (bitpos - bitstart)) >> off;

            bitpos += (m_dataBits - off);
        }
        
        out &= zero;
        return out;
    }
    
    void set(size_t i, T v) {
        size_t bitstart = (i * m_bits);
        size_t bitpos = bitstart;
        
        while(bitpos - bitstart < m_bits) {
            size_t pos = bitpos / m_dataBits;
            size_t off = bitpos % m_dataBits;
            
            size_t rest = m_bits - (bitpos - bitstart);
            D zero = ~((1ul << (rest + off)) - 1) | ((1ul << off) - 1);
            
            m_storage[pos] &= zero;
            m_storage[pos] |= v << off;
            v = v >> (m_dataBits - off);
            bitpos += (m_dataBits - off);
        }
    }
    
    virtual D*& getStorage() {
        return m_storage;
    }
    
    virtual size_t getStorageSize() const {
        return m_storageSize;
    }
    
    virtual size_t size() const {
        return m_size;
    }
    
    virtual size_t load(std::FILE* in) {
        size_t a1 = std::ftell(in);
        
        std::fread(&m_size, sizeof(m_size), 1, in);
        std::fread(&m_storageSize, sizeof(m_storageSize), 1, in);
        delete [] m_storage;
        m_storage = new D[m_storageSize];
        std::fread(m_storage, sizeof(D), m_storageSize, in);
        
        size_t a2 = std::ftell(in);
        return a2 - a1;
    }
    
    virtual size_t save(std::FILE* out) {
        size_t a1 = std::ftell(out);
        
        std::fwrite(&m_size, sizeof(m_size), 1, out);
        std::fwrite(&m_storageSize, sizeof(m_storageSize), 1, out);
        std::fwrite(m_storage, sizeof(D), m_storageSize, out);
        
        size_t a2 = std::ftell(out);
        return a2 - a1;
    }
    
};

template <typename T, typename D>
size_t PackedArray<T, D>::m_bits = 26;

template <typename T, typename D>
size_t PackedArray<T, D>::m_dataBits = sizeof(D)*8;

template <typename T = size_t, typename D = unsigned char>
class PairedPackedArray : public PackedArray<T,D> {
  private:
    static size_t m_bits1;
    static size_t m_bits2;
    
  public:
    PairedPackedArray() : PackedArray<T,D>(m_bits1 + m_bits2, 0) {}
    
    PairedPackedArray(size_t bits1, size_t bits2, size_t size)
    : PackedArray<T, D>(bits1 + bits2, size) {
        m_bits1 = bits1;
        m_bits2 = bits2;
    }
    
    void set(size_t i, T a, T b) {
        T c = 0;
        c = a | (b << m_bits1);
        PackedArray<T,D>::set(i, c);
    }
    
    void set(size_t i, std::pair<T,T> p) {
        T c = 0;
        c = p.second | (p.first << m_bits1);
        PackedArray<T, D>::set(i, c);
    }
    
    std::pair<T, T> get(size_t i) {
        T v = PackedArray<T, D>::get(i);
        T a = v & ((1 << m_bits1) - 1);
        T b = v >> m_bits1;
        return std::pair<T, T>(a, b);
    }
};

template <typename T, typename D>
size_t PairedPackedArray<T, D>::m_bits1 = 10;

template <typename T, typename D>
size_t PairedPackedArray<T, D>::m_bits2 = 16;

}

#endif