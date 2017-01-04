/* -----------------------------------------------------------------------------
 
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
 
----------------------------------------------------------------------------- */


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
