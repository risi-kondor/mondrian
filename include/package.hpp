#ifndef _package
#define _package

namespace Mondrian{

template<class TYPE1, class TYPE2>
class package{
public:

  TYPE1 obj1;
  TYPE2 obj2;

public:

  package(TYPE1& _obj1, TYPE2& _obj2): obj1(std::move(_obj1)), obj2(std::move(_obj2)){}
  //package(const TYPE1& _obj1, const TYPE2& _obj2): obj1(_obj1), obj2(_obj2){}

  package(TYPE1&& _obj1, TYPE2&& _obj2): obj1(std::move(_obj1)), obj2(std::move(_obj2)){}

public:

  package(package&& x): obj1(std::move(x.obj1)), obj2(std::move(x.obj2)){}

  package operator=(package&& x){
    return package(std::move(x.obj1),std::move(x.obj2));}


public:

  TYPE1 first(){return std::move(obj1);}
  TYPE2 second(){return std::move(obj2);}

  package& operator>>(TYPE1& target1){target1=std::move(obj1); return *this;}
  package& operator>>(TYPE2& target2){target2=std::move(obj2); return *this;}

  operator TYPE1(){return std::move(obj1);}
  operator TYPE2(){return std::move(obj2);}

  /*
  class{
  public:
    operator TYPE1(){return std::move(obj1);}
  } first;

  class{
  public:
    operator TYPE1(){return std::move(obj2);}
  } second;
  */

};

} // namespace Mondrian

#endif
