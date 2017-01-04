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


#ifndef _SparseVector
#define _SparseVector

#include "Vector.hpp"

namespace Mondrian{

class SparseVector: public Vector{
public:

  SparseVector(const int _n): Vector(_n){}


public: // attributes

  virtual bool isSparse() const {return true;}


public: // sparse vector methods

  virtual SCALAR* findptr(const int i)=0;
  virtual void insert(const int i, const SCALAR value)=0;
  virtual void append(const int i, const SCALAR value)=0;
  //virtual void zero(const int i)=0;
  virtual void sort() const {}
  virtual void tidy()=0;



};

} // Mondrian

#endif
