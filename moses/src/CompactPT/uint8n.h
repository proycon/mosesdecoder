#ifndef UINT8N_H__
#define UINT8N_H__

#include <iostream>

template<size_t n>
class uint8n {
  private:
    unsigned char bytes[n];
    
  public:
    uint8n() {};
    
    uint8n(long a) {
        for(int i = 1; i <= n; i++) {
            bytes[n-i] = a;
            a = a >> 8;
        }
    }
    
    uint8n& operator=(uint8n& a) {
        return a;
    }
    
    uint8n& operator=(long a) {
        for(int i = 1; i <= n; i++) {
            bytes[n-i] = a;
            a >>= 8;
        }
        return *this;
    }
    
    bool operator==(const uint8n& a) {
        for(int i = 1; i <= n; i++)
            if(bytes[i] != a.bytes[i])
                return false;
        return true;
    }
    
    operator long() const {
        long a = bytes[0];
        for(int i = 1; i < n; i++) {
            a = a << 8;
            a |= bytes[i];
        }
        return a;
    }
    
    std::istream& operator<<(std::istream &is) {
      for(int i = 0; i < n; i++)
          is.get((char&) bytes[i]);
      return is;  
    }
    
    std::ostream& operator>>(std::ostream &os) {
      for(int i = 0; i < n; i++)
          os.put((char&) bytes[i]);
      return os;  
    }
};

typedef uint8n<1> uint8;
typedef uint8n<2> uint16;
typedef uint8n<3> uint24;
typedef uint8n<4> uint32;

#endif
