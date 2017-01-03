#ifndef _ProductOp
#define _ProductOp

#include "LinOp.hpp"
#include <vector>


template<class OP>
class ProductOp: public LinOp{
public:


  ~ProductOp() {for(auto p: op) delete p;}

public:

  template<class VEC> VEC operator*(const VEC& x) const;  

  
public:

  string str() const;

public:

  vector<OP*> op;

};


template<class OP>
template<class VEC>
VEC ProductOp<OP>::operator*(const VEC& x) const {
  if(op.size()==0) return x; 
  VEC y=(*op[0])*x;
  for(int i=1; i<op.size(); i++) y=(*op[i])*y;
  return y;
}


template<class OP>
string ProductOp<OP>::str() const{
  ostringstream oss;
  oss<<"ProductOp";
  return oss.str();
}


#endif
