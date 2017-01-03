#ifndef _DiagBlockedMatrix
#define _DiagBlockedMatrix

#include "Matrix.hpp"
#include "BlockedArray.hpp"


namespace Mondrian{

template<class MATRIX>
class DiagBlockedMatrix: public BlockedMatrix<MATRIX>{
public:



public:


};

template<class BLOCK>
DiagBlockedMatrix<BLOCK> operator*(const SCALAR c, const DiagBlockedMatrix<BLOCK>& M){return M*c;}

}

#endif
