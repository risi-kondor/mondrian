#ifndef _ConstantVec
#define _ConstantVec

#include "Vec.hpp"


class ConstantVec: public Vec{
public:

  ConstantVec(const SCALAR _val=0):val(_val){}


public:

  SCALAR dot(const ConstantVec& x) const {return val*x.val;}


public:
  
  string str() const{
    ostringstream oss; oss<<"ConstantVec("<<val<<")"; return oss.str();}


public:

  SCALAR val;

};


ConstantVec operator*(const SCALAR x, const ConstantVec v){
  return ConstantVec(x*v.val);}
ConstantVec operator*(const ConstantVec v, const SCALAR x){
  return ConstantVec(x*v.val);}


#endif



  //ConstantVec operator*(SCALAR x) const {return ConstantVec(x*val);}
