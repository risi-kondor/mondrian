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


#ifndef _Cluster
#define _Cluster

#include "Mondrian_base.hpp"

namespace Mondrian{

template<class OBJ>
class Cluster{ 
public:


  Cluster(){}
  Cluster(const int _k){vector<OBJ>::resize(_k);}

  vector<OBJ> v;
  
public:

  INDEX operator()(const int i){return v[i];}
  INDEX& operator[](const int i){return v[i];}
  INDEX operator[](const int i) const {return v[i];}

  void push_back(const OBJ x){v.push_back(x);}
  int size() const{return v.size();}

  string str() const{
    ostringstream result; 
    result<<"("; 
    if(v.size()>1) for(int i=0; i<v.size()-1; i++) {result<<v[i];result<<",";}
    if(v.size()>0) result<<v[v.size()-1];
    result<<")"; 
    return result.str();
  }


};

}

#endif
