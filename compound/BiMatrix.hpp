#ifndef _BiMatrix
#define _BiMatrix

#include "Matrix.hpp"

namespace Mondrian{

template<class MATRIX1, class MATRIX2>
class BiMatrix: public Matrix{
public:

  bool type=0;

  MATRIX1 M1;
  MATRIX2 M2;

public:


};

}

#endif
