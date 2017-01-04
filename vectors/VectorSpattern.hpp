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


#ifndef _VectorSpattern
#define _VectorSpattern

#include <set>

#include "Mondrian_base.hpp"
#include "IndexMap.hpp"

namespace Mondrian{


class VectorSpattern{
public:

  set<INDEX> filled;

public:

  VectorSpattern(){}

  VectorSpattern(const initializer_list<INDEX> list){
    for(INDEX v:list) filled.insert(v);
  }

  template<class VECTOR>
  VectorSpattern(const VECTOR& v){
    v.for_each_filled([this](const int i, const SCALAR v){filled.insert(i);});
  }

  template<class VECTOR>
  VectorSpattern(const vector<VECTOR*>& V){
    if (V.size()==0) return;
    V[0]->for_each_filled([this](const int i, const SCALAR v){cout<<i<<endl; filled.insert(i);});
    for(int j=1; j<V.size(); j++)
      V[j]->for_each_filled([this](const int i, const SCALAR v){filled.insert(i);});
  }

  
public: // conversions

  VectorSpattern(const IndexMap& map){
    for(int i=0; i<map.nsource; i++) filled.insert(map(i));
  }

  operator IndexMap(){
    IndexMap map(filled.size());
    std::copy(filled.begin(),filled.end(),map.forward);
    return map;
  }


public: // element access 

  bool isFilled(const int i) const {
    return filled.find(i)!=filled.end();}
    
  
public: //iterators

  void for_each_filled(std::function<void(const INDEX)> lambda){
    for(auto p:filled) lambda(p);}


public: // methods

  template<class VECTOR>
  void add(const VECTOR& v){
    v.for_each_filled([this](const int i, const SCALAR v){filled.insert(i);});
  }

  template<class VECTOR>
  void add(const VECTOR& v, const INDEX limit){
    v.for_each_filled([this,limit](const int i, const SCALAR v){if(i<=limit) filled.insert(i);});
  }


public:

  string str() const{
    ostringstream result; 
    result<<"Spattern( "; 
    for(auto p: filled) result<<p<<" "; 
    result<<")"; 
    return result.str();
  }



};

inline ostream& operator<<(ostream& stream, const VectorSpattern& v){stream<<v.str(); return stream;}



}

#endif
