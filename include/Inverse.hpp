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


#ifndef _Inverse
#define _Inverse

namespace Mondrian{

template<class OBJ>
class Inverse{
public:

  Inverse(const OBJ& _obj): obj(_obj){};
  ~Inverse(){} // what to put here?

  OBJ& obj;

public:

  OBJ& operator~() const {return obj;}
  //OBJ& transp() const {return obj;}


};


}

#endif
