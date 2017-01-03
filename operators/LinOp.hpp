#ifndef _LinOp
#define _LinOp

#include "Vec.hpp"


class LinOp{
public:

  virtual ~LinOp(){};

public:

  virtual string str() const=0;

};


inline ostream& operator<<(ostream& stream, const LinOp& x) {stream<<x.str(); return stream;}


#endif


