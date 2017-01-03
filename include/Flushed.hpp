#ifndef _Flushed
#define _Flushed

namespace Mondrian{

  template<class TYPE> 
  class Locked: public TYPE{
  public:
    
    TYPE& obj;

    Locked(const TYPE& _obj) :obj(_obj){
      obj.access_mx.lock();
    }
    
    ~Locked(){obj.access_mx.unlock();}
    
  };


  template<class TYPE> 
  class Flushed: public TYPE{
  public:
    
    TYPE& obj;
    
    Flushed(const TYPE& _obj) :obj(_obj){
      obj.access_mx.lock();
      obj.flush_while_locked();
    }
    
    ~Flushed(){obj.access_mx.unlock();}
    
  };

  template<class TYPE> 
  Flushed<TYPE> flush(const TYPE& _obj){return _obj;}

  

}


#endif
