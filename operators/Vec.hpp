#ifndef _Vec
#define _Vec

#include "Mondrian_base.hpp"


class Vec{
public:

  Vec(){};

public:

  

public:

  virtual string str() const=0;

};


inline ostream& operator<<(ostream& stream, const Vec& x){
  stream<<x.str(); return stream;}


#endif
