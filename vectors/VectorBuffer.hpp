#ifndef _VectorBuffer
#define _VectorBuffer

#include "Mondrian_base.hpp"
#include "Flushed.hpp"


namespace Mondrian{

  class VectorBuffer{
  public:
    
    vector<ivpair> vec;
    mutex buf_mx;
    mutable mutex access_mx;

    

  public:
    
    

    
  };

  class VectorBuffers{
  public:
    
    
  };

}



#endif
