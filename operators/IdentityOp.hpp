#ifndef _IdentityOp
#define _IdentityOp

#include "LinOp.hpp"
#include "ConstantOp.hpp"


class IdentityOp: public LinOp{
public:

  template<class TYPE> 
  TYPE operator*(const TYPE& x) const {return x;}
  template<class TYPE> 
  TYPE dot(const TYPE& x) const {return x;} 

  ConstantOp operator+(const ConstantOp& x) const {return ConstantOp(x.val+1);}
  ConstantOp operator-(const ConstantOp& x) const {return ConstantOp(1-x.val);}

public:

  string str() const{
    ostringstream oss; oss<<"IdentityOp()"; return oss.str();}

};


ConstantOp operator*(const FIELD x, const IdentityOp v){
  return ConstantOp(x);}
ConstantOp operator*(const IdentityOp v, const FIELD x){ 
  return ConstantOp(x);}



#endif
