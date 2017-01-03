#ifndef _Detached
#define _Detached

namespace Mondrian{

template<class CLASS>
class Detached: public CLASS{
public:

  //using CLASS::CLASS;

  Detached(){}

  Detached(const CLASS& obj): CLASS(obj.shallow()){}

  Detached(const Detached<CLASS>& x): CLASS(x.shallow()){}

  Detached<CLASS>& operator=(const Detached<CLASS>& x){
    *this=x.shallow(); return *this;}

  Detached<CLASS> copy() const {return *this;}

  ~Detached(){CLASS::detach();}

  Detached<CLASS>& operator=(const CLASS& x){
    CLASS::assign(x); return* this;}

  // need to stop upcasts!
  // what about moving?

  //template<class CLASS2>
  //Detached(const CLASS2& x);

};

}

#endif


