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


#ifndef _Flushed
#define _Flushed

namespace Mondrian{

  template<class TYPE> 
  class Locked: public TYPE{
  public:
    
    TYPE& obj;

    Locked(const TYPE& _obj) :obj(_obj){
      obj.access_mx.lock();
    }
    
    ~Locked(){obj.access_mx.unlock();}
    
  };


  template<class TYPE> 
  class Flushed: public TYPE{
  public:
    
    TYPE& obj;
    
    Flushed(const TYPE& _obj) :obj(_obj){
      obj.access_mx.lock();
      obj.flush_while_locked();
    }
    
    ~Flushed(){obj.access_mx.unlock();}
    
  };

  template<class TYPE> 
  Flushed<TYPE> flush(const TYPE& _obj){return _obj;}

  

}


#endif
