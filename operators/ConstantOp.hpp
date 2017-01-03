#ifndef _ConstantOp
#define _ConstantOp

#include "LinOp.hpp"


class ConstantOp: public LinOp{
public:

  ConstantOp(const SCALAR _val=1):val(_val){}

public:

  template<class TYPE> 
  TYPE operator*(const TYPE& x) const {return x*val;} // for TYPE=SCALAR what happens?
  template<class TYPE> 
  TYPE dot(const TYPE& x) const {return x*val;} 

  ConstantOp operator+(const ConstantOp& x) const {return ConstantOp(val+x.val);}
  ConstantOp operator-(const ConstantOp& x) const {return ConstantOp(val-x.val);}

public:

  string str() const{
    ostringstream oss; oss<<"ConstantOp("<<val<<")"; return oss.str();}


public:

  SCALAR val;

};


ConstantOp operator*(const SCALAR x, const ConstantOp v){
  return ConstantOp(x*v.val);}
ConstantOp operator*(const ConstantOp v, const SCALAR x){ 
  return ConstantOp(x*v.val);}



#endif
