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


#ifndef _ConstantOp
#define _ConstantOp

#include "LinOp.hpp"


class ConstantOp: public LinOp{
public:

  ConstantOp(const SCALAR _val=1):val(_val){}

public:

  template<class TYPE> 
  TYPE operator*(const TYPE& x) const {return x*val;} // for TYPE=SCALAR what happens?
  template<class TYPE> 
  TYPE dot(const TYPE& x) const {return x*val;} 

  ConstantOp operator+(const ConstantOp& x) const {return ConstantOp(val+x.val);}
  ConstantOp operator-(const ConstantOp& x) const {return ConstantOp(val-x.val);}

public:

  string str() const{
    ostringstream oss; oss<<"ConstantOp("<<val<<")"; return oss.str();}


public:

  SCALAR val;

};


ConstantOp operator*(const SCALAR x, const ConstantOp v){
  return ConstantOp(x*v.val);}
ConstantOp operator*(const ConstantOp v, const SCALAR x){ 
  return ConstantOp(x*v.val);}



#endif
