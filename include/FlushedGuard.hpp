#ifndef _FlushedGuard
#define _FlushedGuard

namespace Mondrian{

  template<class TYPE>
  class FlushedGuard{
  public:
    
    lock_guard<mutex> lock;
    
    FlushedGuard(TYPE& _owner): lock(_owner.flushedmx){_owner.flush();}
    
  };

}

#endif
