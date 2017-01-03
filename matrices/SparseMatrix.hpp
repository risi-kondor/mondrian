#ifndef _SparseMatrix
#define _SparseMatrix

#include "Matrix.hpp"

namespace Mondrian{

class SparseMatrix: public Matrix{
public:

  using Matrix::Matrix;


public: // attributes

  bool isSparseFormat() const {return true;}


public: // element access

  template<class VECTOR>
  VECTOR row(const int i) const {assert(i<nrows);
    VECTOR x=VECTOR::Zero(ncols); 
    for(int j=0; j<ncols; j++){
      SCALAR v=(*this)(i,j);
      if(v!=0) x(j)=v;} 
    return x;
  }

  template<class VECTOR>
  VECTOR column(const int j) const {assert(j<ncols);
    VECTOR x=VECTOR::Zero(nrows); 
    for(int i=0; i<nrows; i++){
      SCALAR v=(*this)(i,j);
      if(v!=0) x(i)=v;} 
    return x;
  }

  template<class VECTOR>
  VECTOR diag() const {int t=min(nrows,ncols);
    VECTOR x=VECTOR::Zero(t); 
    for(int i=0; i<t; i++){
      SCALAR v=(*this)(i,i);
      if(v!=0) x(i)=v;} 
    return x;
  }





};

} // namespace Mondrian

#endif
