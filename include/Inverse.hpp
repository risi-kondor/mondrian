#ifndef _Inverse
#define _Inverse

namespace Mondrian{

template<class OBJ>
class Inverse{
public:

  Inverse(const OBJ& _obj): obj(_obj){};
  ~Inverse(){} // what to put here?

  OBJ& obj;

public:

  OBJ& operator~() const {return obj;}
  //OBJ& transp() const {return obj;}


};


}

#endif
