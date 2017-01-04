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
