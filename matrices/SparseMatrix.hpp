/* ---------------------------------------------------------------------------
 
  Mondrian: open source multithreaded matrix library 

  Copyright (C) 2017 Imre Risi Kondor

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 3
  of the License, or (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.
 
--------------------------------------------------------------------------- */


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
