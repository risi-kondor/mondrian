#ifndef _OpTrain
#define _OpTrain

#include "MatrixType.hpp"

namespace Mondrian{

  template<class OPTOR>
  class OpTrain: public class MatrixLike{
  public:
    
    vector<OPTOR> op;    

  public:
    
    OpTrain(): Optrain(0,0);
    OpTrain(const int _nrows, const int _ncols): MatrixLike(_nrows,_ncols){}
    
  public: // vector actions 
    
    template<class VECTOR> VECTOR operator*(const VECTOR& x) const;
    template<class VECTOR> VECTOR operator*(VECTOR&& x) const;
    
    template<class VECTOR> VECTOR dot(const VECTOR& x) const;
    template<class VECTOR> VECTOR dot(VECTOR&& x) const;
    
  };


  // ---- Acting on vectors and matrixces --------------------------------------------------------------------------


  template<class OPTOR>
  template<class VECTOR>
  VECTOR OpTrain<OPTOR>::operator*(const VECTOR& x) const {
    if(op.size()==0) return x; 
    VECTOR y=(*op[0])*x;
    for(int i=1; i<op.size(); i++) op[i]->applyTo(y);
    return y;
  }

  template<class OPTOR>
e   template<class VECTOR>
  VECTOR OpTrain<OP>::operator*(VECTOR&& x) const {
    if(op.size()==0) return x; 
    for(int i=0; i<op.size(); i++) op[i]->applyTo(x);
    return x;
  }

  template<class OPTOR>
  template<class VECTOR>
  VECTOR OpTrain<OP>::dot(const VECTOR& x) const {
    if(op.size()==0) return x; 
    VECTOR y=op[op.size()-1]->dot(x);
    for(int i=op.size()-2; i>=0; i--) op[i]->applyTranspTo(y);
    return y;
  }

  template<class OPTOR>
  template<class VECTOR>
  VECTOR OpTrain<OP>::dot(const VECTOR&& x) const {
    if(op.size()==0) return x; 
    for(int i=op.size()-1; i>=0; i--) op[i]->applyTranspTo(x);
    return x;
  }

}


#endif
