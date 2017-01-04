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


#ifndef _IndexBiMap
#define _IndexBiMap

#include "IndexMap.hpp"

namespace Mondrian{


class IndexBiMap: public IndexMap{
public:

  int ndest;
  INDEX* backward;


public:

  IndexBiMap(const int _nsource, const int _ndest=-1): IndexMap(_nsource){
    if(_ndest>=0) ndest=_ndest; else ndest=_nsource;
    backward=new INDEX[ndest];
  }

  IndexBiMap(const initializer_list<INDEX> list): IndexMap(list){
    int ndest=0; for(auto j:list) if(j>ndest) ndest=j;
    backward=new INDEX[ndest];
    std::fill(backward,backward+ndest,-1);
    int i=0; for(auto j:list) {forward[i]=j; backward[j]=i; i++;}
  }


public: // copying

  IndexBiMap(const IndexBiMap& x): IndexMap(x), ndest(x.ndest){
    backward=new INDEX[ndest]; 
    std::copy(x.backward,x.backward+ndest,backward);
  }
  
  IndexBiMap(IndexBiMap&& x): IndexMap(std::move(x)), ndest(x.ndest){
    backward=snatchptr(x.backward);
  }
  
  IndexBiMap& operator=(const IndexBiMap& x){
    nsource=x.nsource; ndest=x.ndest;
    delete forward; forward=new INDEX[nsource]; std::copy(x.forward,x.forward+nsource,forward);
    delete backward; backward=new INDEX[ndest]; std::copy(x.backward,x.backward+ndest,backward);
    return *this;
  }
  
  IndexBiMap& operator=(IndexBiMap&& x) {
    nsource=x.nsource; x.nsource=0; 
    ndest=x.ndest; x.ndest=0; 
    move_over(forward,x.forward); 
    move_over(backward,x.backward); 
    return *this;
  }

  ~IndexBiMap(){delete[] backward;}


public: // named constructors

  static IndexBiMap Identity(const int _nsource){
    IndexBiMap r(_nsource); for(int i=0; i<_nsource; i++) {r.forward[i]=i; r.backward[i]=i;} return r;}


public: // conversions

  IndexBiMap(const IndexMap& x, const int _ndest=-1): IndexMap(x){
    if(_ndest!=-1) ndest=_ndest; else ndest=max();
    backward=new INDEX[ndest]; 
    fixBackwardMap();
  }

  IndexBiMap(const Inverse<IndexMap>& x, const int ndest=-1): IndexMap(0){
    *this=x.obj.inv(ndest);
  }


public: // element access

  INDEX forw(const int i) const {assert(i<nsource); return forward[i];}

  INDEX backw(const int i) const {assert(i<nsource); return backward[i];}

  IndexBiMap& swap(const int i, const int j){
    int t=forward[i]; forward[i]=forward[j]; forward[j]=t; 
    backward[forward[i]]=i; backward[forward[j]]=j;
    return *this;
  }


public: // other methods

  void sort(){
    IndexMap::sort();
    fixBackwardMap();
  }

  void fixBackwardMap(){
    std::fill(backward,backward+ndest,-1);
    for(int i=0; i<nsource; i++) 
      if(forward[i]!=-1) backward[forward[i]]=i;
  }

  IndexBiMap inv(){
    IndexBiMap x(ndest,nsource);
    for(int i=0; i<ndest; i++) x.forward[i]=backward[i];
    for(int i=0; i<nsource; i++) x.backward[i]=forward[i];
    return x;
  }

};

}

#endif
