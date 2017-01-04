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


#ifndef _IdentityOp
#define _IdentityOp

#include "LinOp.hpp"
#include "ConstantOp.hpp"


class IdentityOp: public LinOp{
public:

  template<class TYPE> 
  TYPE operator*(const TYPE& x) const {return x;}
  template<class TYPE> 
  TYPE dot(const TYPE& x) const {return x;} 

  ConstantOp operator+(const ConstantOp& x) const {return ConstantOp(x.val+1);}
  ConstantOp operator-(const ConstantOp& x) const {return ConstantOp(1-x.val);}

public:

  string str() const{
    ostringstream oss; oss<<"IdentityOp()"; return oss.str();}

};


ConstantOp operator*(const FIELD x, const IdentityOp v){
  return ConstantOp(x);}
ConstantOp operator*(const IdentityOp v, const FIELD x){ 
  return ConstantOp(x);}



#endif
