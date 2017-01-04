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


#ifndef _package
#define _package

namespace Mondrian{

template<class TYPE1, class TYPE2>
class package{
public:

  TYPE1 obj1;
  TYPE2 obj2;

public:

  package(TYPE1& _obj1, TYPE2& _obj2): obj1(std::move(_obj1)), obj2(std::move(_obj2)){}
  //package(const TYPE1& _obj1, const TYPE2& _obj2): obj1(_obj1), obj2(_obj2){}

  package(TYPE1&& _obj1, TYPE2&& _obj2): obj1(std::move(_obj1)), obj2(std::move(_obj2)){}

public:

  package(package&& x): obj1(std::move(x.obj1)), obj2(std::move(x.obj2)){}

  package operator=(package&& x){
    return package(std::move(x.obj1),std::move(x.obj2));}


public:

  TYPE1 first(){return std::move(obj1);}
  TYPE2 second(){return std::move(obj2);}

  package& operator>>(TYPE1& target1){target1=std::move(obj1); return *this;}
  package& operator>>(TYPE2& target2){target2=std::move(obj2); return *this;}

  operator TYPE1(){return std::move(obj1);}
  operator TYPE2(){return std::move(obj2);}

  /*
  class{
  public:
    operator TYPE1(){return std::move(obj1);}
  } first;

  class{
  public:
    operator TYPE1(){return std::move(obj2);}
  } second;
  */

};

} // namespace Mondrian

#endif
