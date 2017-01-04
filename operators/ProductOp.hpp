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
