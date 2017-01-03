#ifndef _MatrixSum
#define _MatrixSum

#include "Matrix.hpp"

namespace Mondrian{

  template<class MATRIX>
  class MatrixSum: public Matrix{
  public:
    
    using Matrix::Matrix;

    vector<MATRIX*> matrix;

    
  public: // element access
    
    SCALAR operator()(const int i, const int j) const{ 
      assert(i<nrows); assert(j<ncols); 
      SCALAR t=0; for(auto& p:matrix) t+=p.array[j*nrows+i];
      return t;
    } 
    

  public: // vector actions 

    template<class VECTOR> 
    VECTOR operator*(const VECTOR& x) const{
      if(matrix.size==0) return VECTOR::Zero(x.n);
      VECTOR r=(*matrix[0])*x;
      for(int i=1; i<matrix.size(); i++) r+=(*matrix[i])*x;
      return r;
    }
    
    template<class VECTOR> 
    VECTOR dot(const VECTOR& x) const{
      if(matrix.size==0) return VECTOR::Zero(x.n);
      VECTOR r=matrix[0]->dot(x);
      for(int i=1; i<matrix.size(); i++) r+=matrix[i]->dot(x);
      return r;
    }
    
    




  };

}

#endif
