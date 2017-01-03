#ifndef _HODLRmatrix
#define _HODLRmatrix

#include "BlockedMatrix.hpp"

namespace Mondrian{

template<class MATRIX>
class HODLRmatrix: BlockedMatrix< BiMatrix< HODLRmatrix<MATRIX>, OuterProduct<MATRIX> > >{
public:



public:


};


template<class MATRIX>
HODLRmatrix<MATRIX> operator*(const SCALAR c, const HODLRmatrix<MATRIX>& M){return M*c;}

}

#endif
