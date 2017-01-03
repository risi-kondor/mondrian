#ifndef _Detachable
#define _Detachable

#include "Detached.hpp"

namespace Mondrian{

class Detachable{
public:

  //virtual Detachable shallow()=0;
  virtual void detach()=0;

};

}


#endif
