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


#ifndef _Detached
#define _Detached

namespace Mondrian{

template<class CLASS>
class Detached: public CLASS{
public:

  //using CLASS::CLASS;

  Detached(){}

  Detached(const CLASS& obj): CLASS(obj.shallow()){}

  Detached(const Detached<CLASS>& x): CLASS(x.shallow()){}

  Detached<CLASS>& operator=(const Detached<CLASS>& x){
    *this=x.shallow(); return *this;}

  Detached<CLASS> copy() const {return *this;}

  ~Detached(){CLASS::detach();}

  Detached<CLASS>& operator=(const CLASS& x){
    CLASS::assign(x); return* this;}

  // need to stop upcasts!
  // what about moving?

  //template<class CLASS2>
  //Detached(const CLASS2& x);

};

}

#endif


