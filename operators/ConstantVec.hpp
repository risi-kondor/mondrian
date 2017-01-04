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


#ifndef _ConstantVec
#define _ConstantVec

#include "Vec.hpp"


class ConstantVec: public Vec{
public:

  ConstantVec(const SCALAR _val=0):val(_val){}


public:

  SCALAR dot(const ConstantVec& x) const {return val*x.val;}


public:
  
  string str() const{
    ostringstream oss; oss<<"ConstantVec("<<val<<")"; return oss.str();}


public:

  SCALAR val;

};


ConstantVec operator*(const SCALAR x, const ConstantVec v){
  return ConstantVec(x*v.val);}
ConstantVec operator*(const ConstantVec v, const SCALAR x){
  return ConstantVec(x*v.val);}


#endif



  //ConstantVec operator*(SCALAR x) const {return ConstantVec(x*val);}
